# PEAS (Parametric EIT analysis software): Software to analyze measurements with electrical impedance tomography.
# Copyright (C) 2025 Claas Händel
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import json
import csv
from pathlib import Path
import hashlib

# for some base items
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
from scipy.spatial import ConvexHull

#for preprocessing
from scipy.signal import butter, sosfiltfilt, correlate

from utils.analysisutils import FunctionType
from utils.globalsettings import GlobalSettings

class BaseItem():
    def __init__(self, name, **kwargs):
        self.name = name
        self.kwargs = kwargs

    def generate(self, data, interval):
        return np.array([])

    @staticmethod
    def retrieve_and_sort_flags(interval):
        """
        Retrieve all flags with vertical lines and sort by time.
        """
        flags = []
        for flag in interval.flags:
            if flag.function.function_type == FunctionType.LINE and 'x' in flag.function.params:
                flag.time = flag.function.params['x']
                flags.append(flag)
        flags.sort(key=lambda flag: flag.time)
        return flags

    @staticmethod
    def filter_flags_by_type(flags, flag_types):
        return [flag for flag in flags if flag.flag_type in flag_types]

    @staticmethod
    def filter_flags_by_type_alternatingly(flags, alternating_flag_types):
        filtered_flags = []
        type_index = 0
        i = 0
        while i < len(flags):
            while i < len(flags) and flags[i].flag_type not in alternating_flag_types[type_index]:
                i += 1
            if i < len(flags):
                filtered_flags.append(flags[i])
                type_index = (type_index + 1) % len(alternating_flag_types)
                i += 1
        return filtered_flags

    @staticmethod
    def find_nearest_indices(time_list, times):
        """
        Find the nearest indices in time_list for each time in times.
        """
        return [np.argmin(np.abs(np.array(time_list) - time)) for time in times]

    @staticmethod
    def get_indices_for_flags(time_list, *flags):
        times = [flag.time for flag in flags]
        return BaseItem.find_nearest_indices(time_list, times)

    @staticmethod
    def resample_3d_array(time_list, array_3d, new_time_list):
        """
        Resample a 3D array using Cubic Spline interpolation along dimension 0.
        A separate spline is calculated for every element in dimensions 1 and 2.
        The resampled array is created using the new_time_list.
        """
        _, dim2, dim3 = array_3d.shape

        resampled_array = np.empty((len(new_time_list), dim2, dim3))

        for i in range(dim2):
            for j in range(dim3):
                data = array_3d[:, i, j]
                if np.any(~np.isfinite(data)):
                    resampled_array[:, i, j] = np.nan
                else:
                    cs = CubicSpline(time_list, data)
                    resampled_array[:, i, j] = cs(new_time_list)
        return resampled_array

    @staticmethod
    def initialize_base_item(base_item_name, **kwargs):
        if base_item_name in AVAILABLE_BASE_ITEMS:
            base_item_class = AVAILABLE_BASE_ITEMS[base_item_name]
            return base_item_class(base_item_name, **kwargs)
        else:
            raise ValueError(f"Invalid base item name: {base_item_name}")

class TidalImageItem(BaseItem):
    def __init__(self, name, inspiration=True, **kwargs):
        super().__init__(name, **kwargs)
        self.inspiration = inspiration

    def generate(self, data, interval):
        np_array, time_list = data

        flags = BaseItem.retrieve_and_sort_flags(interval)
        flags = BaseItem.filter_flags_by_type_alternatingly(flags, [GlobalSettings.EOE_FLAG_TYPES, GlobalSettings.EOI_FLAG_TYPES] if self.inspiration else [GlobalSettings.EOI_FLAG_TYPES, GlobalSettings.EOE_FLAG_TYPES])

        if len(flags)==0:
            raise RuntimeError("Empty flags list")

        result_slices = []
        for i in range(0, len(flags)-1, 2):
            index_start, index_end = BaseItem.get_indices_for_flags(time_list, flags[i], flags[i+1])
            result_slice = np_array[index_end] - np_array[index_start]
            result_slice *= 1 if self.inspiration else -1
            result_slices.append(result_slice)
        return np.stack(result_slices)

class RespiredTimeItem(BaseItem):
    def __init__(self, name, t=1, inspiration=False, **kwargs):
        super().__init__(name, **kwargs)
        self.t = t
        self.inspiration = inspiration

    def generate(self, data, interval):
        np_array, time_list = data

        flags = BaseItem.retrieve_and_sort_flags(interval)
        flags = BaseItem.filter_flags_by_type_alternatingly(flags, [GlobalSettings.EOE_FLAG_TYPES, GlobalSettings.EOI_FLAG_TYPES] if self.inspiration else [GlobalSettings.EOI_FLAG_TYPES, GlobalSettings.EOE_FLAG_TYPES])

        if len(flags)==0:
            raise RuntimeError("Empty flags list")

        result_slices = []

        for i in range(1, len(flags), 2):
            time_start = flags[i - 1].time
            time_t = time_start + self.t # TODO: sum might be out of bounds
            index_start, index_t = BaseItem.find_nearest_indices(time_list, [time_start, time_t])

            slice_start = np_array[index_start]
            slice_t = np_array[index_t]

            if self.inspiration:
                result_slice = slice_t - slice_start
            else:
                result_slice = slice_start - slice_t

            result_slices.append(result_slice)
        return np.stack(result_slices)


class RespiredFractionItem(BaseItem):
    """
    time_mode: return times if True, return flows if False
    """
    def __init__(self, name, fraction, inspiration=False, time_mode=True, sf=100, **kwargs):
        super().__init__(name, **kwargs)
        self.fraction = fraction
        self.inspiration = inspiration
        self.time_mode = time_mode
        self.sf = sf

    def generate(self, data, interval):
        np_array, time_list = data

        flags = BaseItem.retrieve_and_sort_flags(interval)
        flags = BaseItem.filter_flags_by_type_alternatingly(flags, [GlobalSettings.EOE_FLAG_TYPES, GlobalSettings.EOI_FLAG_TYPES] if self.inspiration else [GlobalSettings.EOI_FLAG_TYPES, GlobalSettings.EOE_FLAG_TYPES])

        if len(flags)==0:
            raise RuntimeError("Empty flags list")

        result_slices = []

        for i in range(1, len(flags), 2):
            time_start = flags[i-1].time
            time_end = flags[i].time
            index_start, index_end = BaseItem.find_nearest_indices(time_list, [time_start, time_end])

            breath_array = np_array[index_start:index_end]
            resampled_time_list = np.arange(0, time_end-time_start + 1/self.sf, 1/self.sf)

            resampled_breath_array = BaseItem.resample_3d_array(np.linspace(0, time_end-time_start, len(breath_array)).tolist(), breath_array, resampled_time_list.tolist())

            threshold_array = resampled_breath_array[0] + self.fraction*(resampled_breath_array[-1]-resampled_breath_array[0])
            indices = np.argmax(resampled_breath_array>threshold_array if self.inspiration else resampled_breath_array<threshold_array, axis=0)

            if self.time_mode:
                result_slice = resampled_time_list[indices]
            else:
                diff_array = np.diff(resampled_breath_array, axis=0)
                valid_indices = np.minimum(indices, len(diff_array)-1)
                xx, yy = np.meshgrid(np.arange(valid_indices.shape[0]), np.arange(valid_indices.shape[1]), indexing='ij')
                result_slice = diff_array[valid_indices, xx, yy] * self.sf

            result_slices.append(result_slice)
        return np.stack(result_slices)


class MeanFlowItem(BaseItem):
    """
    This class handles the computation of mean flow between two fractions.
    """
    def __init__(self, name, fraction, inspiration=False, sf=100, **kwargs):
        super().__init__(name, **kwargs)
        if not isinstance(fraction, list) or len(fraction) != 2:
            raise ValueError("fraction must be a list with two elements")
        self.fraction = fraction
        self.inspiration = inspiration
        self.sf = sf

    def generate(self, data, interval):
        np_array, time_list = data

        flags = BaseItem.retrieve_and_sort_flags(interval)
        flags = BaseItem.filter_flags_by_type_alternatingly(flags, [GlobalSettings.EOE_FLAG_TYPES, GlobalSettings.EOI_FLAG_TYPES] if self.inspiration else [GlobalSettings.EOI_FLAG_TYPES, GlobalSettings.EOE_FLAG_TYPES])

        if len(flags) == 0:
            raise RuntimeError("Empty flags list")

        result_slices = []

        for i in range(1, len(flags), 2):
            time_start = flags[i-1].time
            time_end = flags[i].time
            index_start, index_end = BaseItem.find_nearest_indices(time_list, [time_start, time_end])

            breath_array = np_array[index_start:index_end]
            resampled_time_list = np.arange(0, time_end-time_start + 1/self.sf, 1/self.sf)

            resampled_breath_array = BaseItem.resample_3d_array(np.linspace(0, time_end-time_start, len(breath_array)).tolist(), breath_array, resampled_time_list.tolist())

            first_threshold = resampled_breath_array[0] + self.fraction[0] * (resampled_breath_array[-1] - resampled_breath_array[0])
            second_threshold = resampled_breath_array[0] + self.fraction[1] * (resampled_breath_array[-1] - resampled_breath_array[0])

            first_indices = np.argmax(resampled_breath_array > first_threshold if self.inspiration else resampled_breath_array < first_threshold, axis=0)
            second_indices = np.argmax(resampled_breath_array > second_threshold if self.inspiration else resampled_breath_array < second_threshold, axis=0)

            first_indices = np.minimum(first_indices, len(resampled_time_list) - 1)
            second_indices = np.minimum(second_indices, len(resampled_time_list) - 1)

            xx, yy = np.meshgrid(np.arange(first_indices.shape[0]), np.arange(first_indices.shape[1]), indexing='ij')

            first_times = resampled_time_list[first_indices]
            second_times = resampled_time_list[second_indices]

            flow_values = (resampled_breath_array[second_indices, xx, yy] - resampled_breath_array[first_indices, xx, yy]) / (second_times - first_times)

            result_slices.append(flow_values)
        return np.stack(result_slices)


class FlowVolumeLoopConcavityItem(BaseItem):
    def __init__(self, name, sf=100, **kwargs):
        super().__init__(name, **kwargs)
        self.sf = sf

    def generate(self, data, interval):

        def polygon_area(x, y):
            return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

        np_array, time_list = data

        flags = BaseItem.retrieve_and_sort_flags(interval)
        flags = BaseItem.filter_flags_by_type_alternatingly(flags, [GlobalSettings.EOI_FLAG_TYPES, GlobalSettings.EOE_FLAG_TYPES])

        if len(flags)==0:
            raise RuntimeError("Empty flags list")

        result_slices = []

        for i in range(1, len(flags), 2):
            time_eoi = flags[i-1].time
            time_eoe = flags[i].time
            index_eoi, index_eoe = BaseItem.find_nearest_indices(time_list, [time_eoi, time_eoe])

            breath_array = np_array[index_eoi:index_eoe]
            resampled_time_list = np.arange(0, time_eoe-time_eoi + 1/self.sf, 1/self.sf)

            resampled_breath_array = BaseItem.resample_3d_array(np.linspace(0, time_eoe-time_eoi, len(breath_array)).tolist(), breath_array, resampled_time_list.tolist())

            _, dim1, dim2 = resampled_breath_array.shape

            result = np.zeros((dim1, dim2))

            for i in range(dim1):
                for j in range(dim2):
                    xi = resampled_breath_array[:, i, j]
                    
                    if np.all(xi == xi[0]) or np.any(np.isnan(xi)):
                        result[i, j] = np.nan
                        continue
                    
                    area_polygon = polygon_area(xi, resampled_time_list)
                    
                    if area_polygon == 0:
                        result[i, j] = np.nan
                        continue
                    
                    points = np.column_stack((xi, resampled_time_list))
                    hull = ConvexHull(points)
                    area_convex_hull = hull.volume
                    
                    result[i, j] = 1 - (area_polygon / area_convex_hull)

            result_slices.append(result)
        return np.stack(result_slices)


class RegressionItem(BaseItem):
    def __init__(self, name, func_str, variables, valid_range=(.25, .75), inspiration=False, **kwargs):
        self.p0 = kwargs.pop('p0', None)
        self.bounds = kwargs.pop('bounds', (-np.inf, np.inf))
        super().__init__(name, **kwargs)
        self.func_str = func_str
        self.variables = variables
        self.valid_range = valid_range
        self.inspiration = inspiration

    def parse_function(self, func_str, vars):
        """
        Safely parses a function string into a callable function.

        :param func_str: String representation of the function.
        :param vars: List of variable names in the function.
        :return: Callable function.
        """
        eval_globals = {"__builtins__": {}, 'np': np}
        return eval(f"lambda {', '.join(vars)}: {func_str}", eval_globals)

    def create_valid_mask(self, z_data, valid_range, inspiration):
        """
        Creates a mask based on thresholds within a given range and depending on the inspiration flag.

        Parameters:
        z_data (numpy array): The data array of shape (t, x, y).
        valid_range (tuple): A tuple (first, second) defining the valid range.
        inspiration (bool): A boolean flag that changes the comparison for threshold.

        Returns:
        numpy array: A mask of shape (t, x, y) based on the threshold conditions.
        """

        valid_range = sorted(valid_range)

        first_threshold = z_data[0] + valid_range[0] * (z_data[-1] - z_data[0])
        second_threshold = z_data[0] + valid_range[1] * (z_data[-1] - z_data[0])

        mask = np.zeros_like(z_data, dtype=bool)

        for x in range(z_data.shape[1]):
            for y in range(z_data.shape[2]):

                if inspiration:
                    first_index = np.argmax(z_data[:, x, y] > first_threshold[x, y])
                else:
                    first_index = np.argmax(z_data[:, x, y] < first_threshold[x, y])

                if inspiration:
                    second_index = np.argmax(z_data[:, x, y] > second_threshold[x, y])
                else:
                    second_index = np.argmax(z_data[:, x, y] < second_threshold[x, y])

                mask[first_index:second_index, x, y] = True
        return mask

    def fit_curve_to_masked_data(self, z_data, time, valid_mask, fit_func):
        """
        Performs a curve fit on each x/y element of z_data using the corresponding time data.
        Only data points where valid_mask is True are considered for the curve fit.
        Outputs all results of the fit (including R^2) in a numpy array.

        Parameters:
        z_data (numpy array): The data array of shape (t, x, y).
        time (numpy array): The time array of shape (t).
        valid_mask (numpy array): A boolean mask of shape (t, x, y) indicating valid data points.

        Returns:
        numpy array: An array containing the fit parameters and R^2 for each (x, y) point.
        """

        n_results = len(self.variables)

        fit_results = np.full((z_data.shape[1], z_data.shape[2], n_results), np.nan)

        for x in range(z_data.shape[1]):
            for y in range(z_data.shape[2]):
                masked_data = z_data[valid_mask[:, x, y], x, y]
                masked_time = time[valid_mask[:, x, y]]

                if len(masked_data) > 1:
                    try:
                        params, _ = curve_fit(fit_func, masked_time, masked_data, p0=self.p0, bounds=self.bounds)
                        residuals = masked_data - fit_func(masked_time, *params)
                        ss_res = np.sum(residuals**2)
                        ss_tot = np.sum((masked_data - np.mean(masked_data))**2)
                        r_squared = 1 - (ss_res / ss_tot)

                        fit_results[x, y, :len(self.variables)-1] = params
                        fit_results[x, y, len(self.variables)-1] = r_squared
                    except (RuntimeError, TypeError):
                        continue
                else:
                    continue
        return fit_results

    def generate(self, data, interval):
        np_array, time_list = data

        flags = BaseItem.retrieve_and_sort_flags(interval)
        flags = BaseItem.filter_flags_by_type_alternatingly(flags, [GlobalSettings.EOE_FLAG_TYPES, GlobalSettings.EOI_FLAG_TYPES] if self.inspiration else [GlobalSettings.EOI_FLAG_TYPES, GlobalSettings.EOE_FLAG_TYPES])

        if len(flags)==0:
            raise RuntimeError("Empty flags list")

        result_slices = []

        parsed_func = self.parse_function(self.func_str, self.variables)

        for i in range(1, len(flags), 2):
            index_start, index_end = BaseItem.get_indices_for_flags(time_list, flags[i-1], flags[i])

            breath_array = np_array[index_start:index_end]
            time_array = time_list[index_start:index_end]-time_list[index_start]

            valid_mask = self.create_valid_mask(breath_array, self.valid_range, self.inspiration)
            
            # perform fits
            result_slice = self.fit_curve_to_masked_data(breath_array, time_array, valid_mask, parsed_func)

            result_slices.append(result_slice)

        return np.stack(result_slices)

class PeakFlowItem(BaseItem):
    def __init__(self, name, inspiration=False, time_mode=False, sf=100, **kwargs):
        super().__init__(name, **kwargs)
        self.inspiration = inspiration
        self.time_mode = time_mode
        self.sf = sf

    def generate(self, data, interval):
        np_array, time_list = data

        flags = BaseItem.retrieve_and_sort_flags(interval)
        flags = BaseItem.filter_flags_by_type_alternatingly(flags, [GlobalSettings.EOE_FLAG_TYPES, GlobalSettings.EOI_FLAG_TYPES] if self.inspiration else [GlobalSettings.EOI_FLAG_TYPES, GlobalSettings.EOE_FLAG_TYPES])

        if len(flags)==0:
            raise RuntimeError("Empty flags list")

        result_slices = []

        for i in range(1, len(flags), 2):
            time_start = flags[i-1].time
            time_end = flags[i].time
            index_start, index_end = BaseItem.find_nearest_indices(time_list, [time_start, time_end])

            breath_array = np_array[index_start:index_end]
            time_list = np.linspace(0, time_end-time_start, len(breath_array)).tolist()
            resampled_time_list = np.arange(0, time_end-time_start + 1/self.sf, 1/self.sf)

            resampled_breath_array = BaseItem.resample_3d_array(time_list, breath_array, resampled_time_list.tolist())

            diff_array = np.diff(resampled_breath_array, axis=0)

            if self.time_mode:
                if self.inspiration:
                    result_slice = resampled_time_list[np.argmax(diff_array, axis=0)] 
                else:
                    result_slice = resampled_time_list[np.argmin(diff_array, axis=0)]
            else:
                if self.inspiration:
                    result_slice = np.max(diff_array, axis=0)*self.sf
                else:
                    result_slice = np.min(diff_array, axis=0)*self.sf

            result_slices.append(result_slice)

        return np.stack(result_slices)

class AbsoluteImpedanceItem(BaseItem):
    def __init__(self, name, inspiration=False, **kwargs):
        super().__init__(name, **kwargs)
        self.inspiration = inspiration

    def generate(self, data, interval):
        np_array, time_list = data

        flags = BaseItem.retrieve_and_sort_flags(interval)
        flags = BaseItem.filter_flags_by_type(flags, GlobalSettings.EOI_FLAG_TYPES if self.inspiration else GlobalSettings.EOE_FLAG_TYPES)

        if len(flags)==0:
            raise RuntimeError("Empty flags list")

        result_slices = []

        for i in range(len(flags)):
            time = flags[i].time
            index = BaseItem.find_nearest_indices(time_list, [time])[0]

            result_slice = np_array[index]
            result_slices.append(result_slice)

        return np.stack(result_slices)

class BreathTimesItem(BaseItem):
    def __init__(self, name, **kwargs):
        super().__init__(name, **kwargs)

    def generate(self, data, interval):
        flags = BaseItem.retrieve_and_sort_flags(interval)
        flags = BaseItem.filter_flags_by_type_alternatingly(flags, [GlobalSettings.EOE_FLAG_TYPES, GlobalSettings.EOI_FLAG_TYPES])

        if len(flags)==0:
            raise RuntimeError("Empty flags list")

        times = np.array([flag.time for flag in flags])
        differences = np.diff(times)
        differences = np.append(differences, np.nan)

        final_shape = (len(flags) // 2, 1, 1, 3)
        reshaped_differences = np.zeros(final_shape)
        for i in range(final_shape[0]):
            even_diff = differences[2 * i]
            odd_diff = differences[2 * i + 1]
            reshaped_differences[i, 0, 0, 0] = even_diff
            reshaped_differences[i, 0, 0, 1] = odd_diff
            reshaped_differences[i, 0, 0, 2] = even_diff + odd_diff

        return reshaped_differences

class IntervalDataItem(BaseItem):
    def __init__(self, name, attr_name):
        super().__init__(name)
        self.attr_name = attr_name

    def generate(self, data, interval):
        return getattr(interval, self.attr_name, None)

class DetectorDataItem(BaseItem):
    def __init__(self, name, attr_name):
        super().__init__(name)
        self.attr_name = attr_name

    def generate(self, data, interval):
        if hasattr(interval.detector, self.attr_name): # is attribute of detector class?
            return getattr(interval.detector, self.attr_name)
        detector_setting = interval.detector.settings.get(self.attr_name, None)
        if detector_setting is None:
            return np.nan
        else:
            return detector_setting['value']

class PassthroughItem(BaseItem):
    def __init__(self, name, **kwargs):
        super().__init__(name, **kwargs)

    def generate(self, data, interval):
        np_array, _ = data
        return np_array

class ROIItem(BaseItem):
    def __init__(self, name, **kwargs):
        super().__init__(name, **kwargs)

    def generate(self, data, interval):
        mask = np.array([[
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0],
            [0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0],
            [0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0],
            [0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0],
            [0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0],
            [0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0],
            [0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
            [0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
            [0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0],
            [0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0],
            [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
            [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
            [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
            [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
            [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
            [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
            [0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0],
            [0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
            [0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0],
            [0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        ]], dtype=bool)
        return mask

AVAILABLE_BASE_ITEMS = {
    "tidal_image": TidalImageItem,
    "expiratory_tidal_image": lambda name, **kwargs: TidalImageItem(name, inspiration=False, **kwargs),
    "expired_volume": RespiredTimeItem,
    "inspired_volume": lambda name, **kwargs: RespiredTimeItem(name, inspiration=True, **kwargs),
    "time_to_expire": RespiredFractionItem,
    "time_to_inspire": lambda name, **kwargs: RespiredFractionItem(name, inspiration=True, **kwargs),
    "flow_when_expired": lambda name, **kwargs: RespiredFractionItem(name, time_mode=False, **kwargs),
    "flow_when_inspired": lambda name, **kwargs: RespiredFractionItem(name, inspiration=True, time_mode=False, **kwargs),
    "mean_expiratory_flow": MeanFlowItem,
    "peak_expiratory_flow": PeakFlowItem,
    "peak_expiratory_flow_time": lambda name, **kwargs: PeakFlowItem(name, time_mode=True, **kwargs),
    "expiratory_concavity": FlowVolumeLoopConcavityItem,
    "end_expiratory_lung_impedance": AbsoluteImpedanceItem,
    "end_inspiratory_lung_impedance": lambda name, **kwargs: AbsoluteImpedanceItem(name, inspiration=True, **kwargs),
    "expiratory_time_constant": lambda name, **kwargs: RegressionItem(name, 'z0 * np.exp(-t/tau) + c', ['t', 'z0', 'tau', 'c'], bounds=(0, np.inf), **kwargs),
    "custom_fit": RegressionItem,
    "breath_times": BreathTimesItem,
    "interval_data": lambda name, attr_name: IntervalDataItem(name, attr_name),
    "detector_data": lambda name, attr_name: DetectorDataItem(name, attr_name),
    "passthrough": PassthroughItem,
    "thorax_roi": ROIItem
}


class Operation():
    def __init__(self, name: str, function, **kwargs):
        self.name = name
        self.function = function
        self.kwargs = kwargs

    def apply(self, *input_data):
        with np.errstate(divide='ignore'):
            if all(not np.any(np.isfinite(item)) for item in input_data):
                # return data unchanged if all inputs empty
                # if some are empty, let the operation handle NaNs
                return input_data[0]
            return self.function(*input_data, **self.kwargs)

    @staticmethod
    def initialize_operation(operation_name, **kwargs):
        if operation_name in AVAILABLE_OPERATIONS:
            operation_class = AVAILABLE_OPERATIONS[operation_name]
            return operation_class(**kwargs)
        else:
            raise ValueError(f"Invalid operation name: {operation_name}")

def mean_over_time(input_data, **kwargs):
    original_shape = input_data.shape
    result = np.nanmean(input_data, axis=0, **kwargs)
    return reshape_with_new_axis(result, original_shape, axis=0)

def median_over_time(input_data, **kwargs):
    original_shape = input_data.shape
    result = np.nanmedian(input_data, axis=0, **kwargs)
    return reshape_with_new_axis(result, original_shape, axis=0)

def mean_over_image(input_data, **kwargs):
    original_shape = input_data.shape
    result = np.nanmean(input_data, axis=(1, 2), **kwargs)
    return reshape_with_new_axis(result, original_shape, axis=(1, 2))

def sum_over_image(input_data, **kwargs):
    original_shape = input_data.shape
    result = np.nansum(input_data, axis=(1, 2), **kwargs)

    non_finite_mask = ~np.any(np.isfinite(input_data), axis=(1, 2))
    result[non_finite_mask] = np.nan

    return reshape_with_new_axis(result, original_shape, axis=(1, 2))

def percentile_over_image(input_data, percentile, **kwargs):
    original_shape = input_data.shape
    result = np.nanpercentile(input_data, percentile, axis=(1, 2), **kwargs)
    return reshape_with_new_axis(result, original_shape, axis=(1, 2))

def threshold_operation(input_data, threshold, invert):
    if invert:
        result = input_data < threshold
    else:
        result = input_data > threshold
    return result[..., np.newaxis] if result.ndim == 2 else result

def percentile_operation(input_data, percentile_value):
    threshold_value = np.nanpercentile(input_data, percentile_value)
    result = input_data > threshold_value
    return result[..., np.newaxis] if result.ndim == 2 else result

def cumulative_threshold(input_data, threshold):
    if input_data.ndim != 3 or input_data.shape[0] != 1:
        raise ValueError("Input data must be a 3D array (time, space, space) with the first dimension of length 1.")

    flattened = input_data.reshape(-1)
    flattened = np.nan_to_num(flattened, nan=0)

    sorted_indices = np.argsort(flattened)[::-1]
    sorted_values = flattened[sorted_indices]
    cumulative_sum = np.nancumsum(sorted_values)

    n_smallest = np.searchsorted(cumulative_sum, threshold, side='left')

    output_mask = np.zeros_like(flattened, dtype=bool)
    output_mask[sorted_indices[:n_smallest]] = True

    return output_mask.reshape(input_data.shape)

def coefficient_of_variation(input_data, **kwargs):
    if input_data.ndim != 3:
        raise ValueError("Input data must be a 3D array (time, space, space).")

    axes = kwargs.get("axis", tuple(range(1, input_data.ndim))) # default: apply frame-wise
    result = np.nanstd(input_data, axis=axes) / np.abs(np.nanmean(input_data, axis=axes))
    return result[np.newaxis, ...] if result.ndim == 0 else result

def global_inhomogeneity_index(input_data, **kwargs):
    if input_data.ndim != 3:
        raise ValueError("Input data must be a 3D array (time, space, space).")

    axes = kwargs.get("axis", tuple(range(1, input_data.ndim))) # default: apply frame-wise
    total_weight = np.nansum(input_data, axis=axes, keepdims=True)
    median_di = np.nanmedian(input_data, axis=axes, keepdims=True)
    abs_diff = np.abs(input_data - median_di)
    sum_abs_diff = np.nansum(abs_diff, axis=axes, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        gi = np.where(
            total_weight != 0,
            sum_abs_diff / total_weight,
            0
        )
    return gi

def apply_mask(np_array, bool_mask):
    if np_array.shape[1:] != bool_mask.shape[1:]:
        raise ValueError(f"Dimensions of np_array ({np_array.shape}) and bool_mask ({bool_mask.shape}) do not match.")

    expanded_mask = bool_mask == False
    expanded_mask = np.broadcast_to(expanded_mask, np_array.shape)
    result = np_array.copy()
    result[expanded_mask] = np.nan

    return result

def normalize_sum(input_data):
    if input_data.ndim != 3 or input_data.shape[0] != 1:
        raise ValueError("Input data must be a 3D array (time, space, space) with the first dimension of length 1.")
    total_sum = np.nansum(input_data)
    result = input_data / total_sum if total_sum != 0 else input_data
    return result[..., np.newaxis] if result.ndim == 2 else result

def normalize_max(input_data):
    if input_data.ndim != 3 or input_data.shape[0] != 1:
        raise ValueError("Input data must be a 3D array (time, space, space) with the first dimension of length 1.")
    max_value = np.nanmax(input_data)
    result = input_data / max_value if max_value != 0 else input_data
    return result[..., np.newaxis] if result.ndim == 2 else result

def centroid(input_data):
    if input_data.ndim != 3:
        raise ValueError("Input data must be a 3D array (time, space, space).")

    t, h, w = input_data.shape
    out = np.zeros((t, 1, 1, 2), dtype=float)

    axes = (1, 2)

    total_weight = np.nansum(input_data, axis=axes, keepdims=True)

    y_coords, x_coords = np.meshgrid(
        np.arange(h, dtype=float),
        np.arange(w, dtype=float),
        indexing="ij"
    )

    y_center = np.nansum(input_data * y_coords[None, :, :], axis=axes, keepdims=True) / total_weight
    x_center = np.nansum(input_data * x_coords[None, :, :], axis=axes, keepdims=True) / total_weight

    denom_y = (h - 1) if h > 1 else np.nan
    denom_x = (w - 1) if w > 1 else np.nan

    y_fraction = y_center / denom_y
    x_fraction = x_center / denom_x

    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.stack([y_fraction, x_fraction], axis=-1)
        out = np.where(total_weight[..., None] != 0, out, 0)

    return out

def slice_last_operation(input_data, index):
    """ use index to slice last dimension of data
    """
    if index < -input_data.shape[input_data.ndim-1] or index >= input_data.shape[input_data.ndim-1]:
        raise ValueError("Index out of range for the last dimension of the input data.")
    return input_data[..., index]

def slice_first_operation(input_data, index):
    """ use index to slice first dimension of data
    """
    if index < -input_data.shape[0] or index >= input_data.shape[0]:
        raise ValueError("Index out of range for the first dimension of the input data.")
    return input_data[index:index+1 or None] # None is used if index=-1

def custom_stack_operation(*arrays, **kwargs):
    # Determine the target shape for broadcasting by finding the first valid input shape
    target_shape = None
    for array in arrays:
        if array.shape != (1, 1, 1):
            target_shape = array.shape
            break
    
    if target_shape is None:
        target_shape = (1, 1, 1)

    broadcasted_arrays = []
    for array in arrays:
        if array.shape == (1, 1, 1) and np.isnan(array).all():
            broadcasted_arrays.append(np.full(target_shape, np.nan))
        else:
            broadcasted_arrays.append(array)

    return np.concatenate(broadcasted_arrays, axis=0, **kwargs)

def reshape_with_new_axis(result, original_shape, axis):
    if axis is None:
        return result
    axis = (axis,) if not isinstance(axis, tuple) else axis
    reduced_shape = [original_shape[i] for i in range(len(original_shape)) if i not in axis]
    for ax in sorted(axis):
        reduced_shape.insert(ax, 1)

    return result.reshape(reduced_shape)

def min_operation(input_data, **kwargs):
    """ Get scalar minimum of all input_data
    """
    original_shape = input_data.shape
    result = np.nanmin(input_data, **kwargs)
    return reshape_with_new_axis(result, original_shape, kwargs.get('axis'))

def max_operation(input_data, **kwargs):
    """ Get scalar maximum of all input_data
    """
    original_shape = input_data.shape
    result = np.nanmax(input_data, **kwargs)
    return reshape_with_new_axis(result, original_shape, kwargs.get('axis'))

def minimum_operation(input_data, x2):
    """ Get element-wise minimum when comparing with x2
    """
    return np.minimum(input_data, x2)

def maximum_operation(input_data, x2):
    """ Get element-wise maximum when comparing with x2
    """
    return np.maximum(input_data, x2)

def costa_approach_operation(input_array):
    max_values = np.nanmax(input_array, axis=0)
    difference_normalized = (max_values - input_array) / np.nansum(max_values)

    before_max = np.full_like(input_array, np.nan)
    after_max = np.full_like(input_array, np.nan)

    spatial_shape = input_array.shape[1:]
    it = np.ndindex(spatial_shape)

    for idx in it:
        try:
            max_idx = np.nanargmax(input_array[:, idx[0], *idx[1:]])
        except ValueError:
            continue  # all-NaN time series at this position
        for n in range(input_array.shape[0]):
            if n <= max_idx:
                before_max[(n,) + idx] = difference_normalized[(n,) + idx]
            if n >= max_idx:
                after_max[(n,) + idx] = difference_normalized[(n,) + idx]

    return np.stack((before_max, after_max), axis=-1)

AVAILABLE_OPERATIONS =  {
    "mean_over_time": lambda: Operation("mean_over_time", mean_over_time),
    "median_over_time": lambda: Operation("median_over_time", median_over_time),
    "mean_over_image": lambda: Operation("mean_over_image", mean_over_image),
    "median_over_image": lambda: Operation("percentile_over_image", mean_over_image, percentile=0.5),
    "sum_over_image": lambda: Operation("sum_over_image", sum_over_image),
    "percentile_over_image": lambda percentile: Operation("percentile_over_image", mean_over_image, percentile=percentile),
    "threshold": lambda threshold, invert=False: Operation("threshold", threshold_operation, threshold=threshold, invert=invert),
    "percentile": lambda percentile: Operation("percentile", percentile_operation, percentile=percentile),
    "cumulative_threshold": lambda threshold: Operation("cumulative_threshold", cumulative_threshold, threshold=threshold),
    "sum": lambda **kwargs: Operation("sum", np.nansum, **kwargs),
    "std": lambda **kwargs: Operation("std", np.nanstd, **kwargs),
    "coefficient_of_variation": lambda **kwargs: Operation("coefficient_of_variation", coefficient_of_variation, **kwargs),
    "global_inhomogeneity_index": lambda **kwargs: Operation("global_inhomogeneity_index", global_inhomogeneity_index, **kwargs),
    "normalize_sum": lambda: Operation("normalize_sum", normalize_sum),
    "normalize_max": lambda: Operation("normalize_max", normalize_max),
    "centroid": lambda: Operation("centroid", centroid),
    "multiply": lambda: Operation("multiply", np.multiply),
    "divide": lambda: Operation("divide", np.divide),
    "subtract": lambda: Operation("subtract", np.subtract),
    "add": lambda: Operation("add", np.add),
    "slice_last": lambda index: Operation("slice_last", slice_last_operation, index=index),
    "slice_first": lambda index: Operation("slice_first", slice_first_operation, index=index),
    "stack": lambda: Operation("stack", custom_stack_operation),
    "size": lambda: Operation("size", lambda arr: len(arr)),
    "min": lambda **kwargs: Operation("min", min_operation, **kwargs), # return minimum of array
    "max": lambda **kwargs: Operation("max", max_operation, **kwargs),
    "minimum": lambda x2: Operation("minimum", minimum_operation, x2=x2), # return element-wise minimum
    "maximum": lambda x2: Operation("maximum", maximum_operation, x2=x2),
    "apply_mask": lambda: Operation("apply_mask", apply_mask),
    "costa_approach": lambda: Operation("costa_approach", costa_approach_operation)
}

class Preprocessor():
    def __init__(self, name: str, function, **kwargs):
        self.name = name
        self.function = function
        self.kwargs = kwargs

    def apply(self, interval, *input_data):
        return self.function(interval, *input_data, **self.kwargs)

    @staticmethod
    def initialize_preprocessor(preprocessor_name, **kwargs):
        if preprocessor_name in AVAILABLE_PREPROCESSORS:
            preprocessor_class = AVAILABLE_PREPROCESSORS[preprocessor_name]
            return preprocessor_class(**kwargs)
        else:
            raise ValueError(f"Invalid preprocessor name: {preprocessor_name}")

def resample_discrete(interval, input_data, n, m):
    """
    Resamples a 3D numpy array (t, x, y) to new dimensions (n, m) for the last two axes using weighted averaging.
    
    :param input_data: tuple of a 3D numpy array with shape (t, y, x) and a list of times length t.
    :param n, m: A int n, m representing new dimensions for the y (v/d) and x (r/l) axes.
    :return: A resampled 3D numpy array with shape (t, n, m).
    """
    np_array, time_list = input_data
    np_array = np.nan_to_num(np_array, nan=0)
    t, x, y = np_array.shape

    x_cell_size = 1.0 / x
    y_cell_size = 1.0 / y
    n_cell_size = 1.0 / n
    m_cell_size = 1.0 / m

    resampled_arr = np.zeros((t, n, m))

    for i in range(n):
        for j in range(m):
            x_start = i * n_cell_size
            x_end = (i + 1) * n_cell_size
            y_start = j * m_cell_size
            y_end = (j + 1) * m_cell_size

            weighted_sum = 0
            total_weight = 0

            for xi in range(x):
                for yj in range(y):
                    x_overlap = max(0, min(x_end, (xi + 1) * x_cell_size) - max(x_start, xi * x_cell_size))
                    y_overlap = max(0, min(y_end, (yj + 1) * y_cell_size) - max(y_start, yj * y_cell_size))
                    overlap_area = x_overlap * y_overlap

                    if overlap_area > 0:
                        weighted_sum += np_array[:, xi, yj] * overlap_area
                        total_weight += overlap_area

            resampled_arr[:, i, j] = weighted_sum / total_weight
    return resampled_arr, time_list

def filter_butterworth(interval, input_data, cutoff, filter_order):
    np_array, time_list = input_data

    fs = 1 / np.mean(np.diff(time_list))

    sos = butter(filter_order, cutoff, fs=fs, output='sos')
    filtered_data = sosfiltfilt(sos, np_array, axis=0)

    return filtered_data, time_list


def strip_all_nan(interval, input_data):
    np_array, time_list = input_data

    row_valid = ~np.all(np.isnan(np_array), axis=(0, 2))
    col_valid = ~np.all(np.isnan(np_array), axis=(0, 1))

    if not row_valid.any() or not col_valid.any():
        return np.full((len(time_list), 1, 1), np.nan), time_list

    r0, r1 = np.where(row_valid)[0][[0, -1]]
    c0, c1 = np.where(col_valid)[0][[0, -1]]

    return np_array[:, r0:r1 + 1, c0:c1 + 1], time_list

def breath_averaging(interval, input_data, center_inspiration, max_lag):
    def mean_shifted_copies(np_array, indices):
        shifted_arrays = []
        
        for idx in indices:
            shifted_array = np.roll(np_array, shift=idx, axis=0)
            shifted_arrays.append(shifted_array)
        
        mean_shifted_array = np.mean(shifted_arrays, axis=0)
        
        return mean_shifted_array

    np_array, time_list = input_data

    time_diff = np.diff(time_list)
    mean_time_diff = np.mean(time_diff)
    max_lag_indices = int(max_lag / mean_time_diff)

    flags = BaseItem.retrieve_and_sort_flags(interval)
    flags = BaseItem.filter_flags_by_type_alternatingly(
        flags, 
        [GlobalSettings.EOE_FLAG_TYPES, GlobalSettings.EOI_FLAG_TYPES] if center_inspiration else [GlobalSettings.EOI_FLAG_TYPES, GlobalSettings.EOE_FLAG_TYPES]
    )

    if len(flags) == 0:
        raise RuntimeError("Empty flags list")

    snippets = []
    indices_mid = []
    for i in range(len(flags) - 2):
        index_start, index_mid, index_end = BaseItem.get_indices_for_flags(time_list, flags[i], flags[i+1], flags[i+2])
        snippet = np_array[index_start:index_end]
        snippets.append(snippet)
        indices_mid.append(index_mid)

    ref_center = indices_mid[0]
    shifts = [index_mid - ref_center for index_mid in indices_mid]

    reference_snippet = mean_shifted_copies(np_array, shifts)

    lags = []
    for snippet, index_mid in zip(snippets, indices_mid):
        correlation = correlate(snippet, reference_snippet, mode='full')
        lag = np.argmax(correlation) - len(snippet) + 1
        lag = np.clip(lag, -max_lag_indices, max_lag_indices)
        adjusted_lag = index_mid - ref_center + lag
        lags.append(adjusted_lag)

    mean_signal = mean_shifted_copies(np_array, lags)
    return mean_signal

AVAILABLE_PREPROCESSORS = {
        "resample_over_image": lambda vd, rl: Preprocessor("resample_over_image", resample_discrete, n=vd, m=rl),
        "low_pass_filter": lambda cutoff, filter_order: Preprocessor("low_pass_filter", filter_butterworth, cutoff=cutoff, filter_order=filter_order),
        "strip_all_nan": lambda: Preprocessor("strip_all_nan", strip_all_nan)
}

class AnalysisItem():
    """ Object that defines and performs an analysis step. Inputs can be other analysis items or base items.
    Base items are created if required.
    Calculation is done by working along a chain of operations, which are chained to earlier operations or inputs.
    Results are saved within the object to prevent repeated calculations.
    """
    def __init__(self, name=None, base_item=None, parameters={}, interval=None, preprocessors=[], prerequisites=[], operations=[], title=None, unit=None, identifier=None, comment=None, export=False):
        # interval is name of interval or None
        # prerequisites is list of dicts with name, parameters as kwargs
        # operations is list of dicts with name, parameters as kwargs
        self.name = name # internal name, also used to create file names
        self.title = title # optional name string
        self.unit = unit # optional unit string
        self.identifier = identifier # optional unique identifier (LOINC, SNOMED CT) string
        self.comment = comment # optional comment string

        # only set if base item:
        self.base_item = base_item
        self.parameters = parameters
        self.interval = interval
        self.preprocessors = preprocessors

        self.hash = self.get_hash()

        # only set if composed item - these are dicts with infos
        self.prerequisites = prerequisites
        self.operations = operations
        
        # these are actual objects
        self.prerequisites_data = [None]*len(self.prerequisites)
        self.operations_data = [None]*len(self.operations)

        self.export = export

        self.manager = None # will be set later
        self.result = None # will be calculated on demand

    def set_manager(self, manager):
        self.manager = manager

    def get_hash(self):
        json_dict = {
                "name": self.name,
                "base_item": self.base_item,
                "parameters": self.parameters,
                "interval": self.interval,
                "preprocessors": self.preprocessors
                }
        return hashlib.sha256(json.dumps(json_dict, sort_keys=True).encode()).hexdigest()

    def to_dict(self):
        data = {
            "name": self.name,
            "title": self.title,
            "unit": self.unit,
            "identifier": self.identifier,
            "comment": self.comment,
            "base_item": self.base_item,
            "parameters": self.parameters,
            "interval": self.interval,
            "preprocessors": self.preprocessors,
            "prerequisites": self.prerequisites,
            "operations": self.operations,
            "export": self.export
        }
        return {key: value for key, value in data.items() if value}

    def get_results(self):
        # return stored results or calculate on demand
        # results cannot be None
        if self.result is None:
            self.calculate()
        # TODO: replace assert with exception
        assert self.result is not None
        return self.result

    def wipe(self):
        self.result = None
        self.prerequisites_data = [None for _ in self.prerequisites_data]
        self.operations_data = [None for _ in self.operations_data]

    def calculate(self):
        # analysis items can be pre-defined (by user e.g. via json) or dynamic (created as naked base items)
        for index, prerequisite_info in enumerate(self.prerequisites):
            # * has name? Get item with name. Must exist.
            # * has no name? create item. Get hash. Replace with existing if exists. Add to manager anyway.
            # generate new item including hash
            if prerequisite_info.get("name", None) is not None:
                item = self.manager.get_item_by_name(prerequisite_info['name'])
                if item is None:
                    raise RuntimeError(f"Item {prerequisite_info['name']} does not exist")
            else:
                item = AnalysisItem(**prerequisite_info)
                item = self.manager.analysis_items.get(item.hash, item)
                self.manager.add_analysis_item(item)
            self.prerequisites_data[index] = item.get_results()
            # all prerequisites are met.

        if self.name is None:
            # this is a base item
            interval = self.manager.interval_lookup_dict.get(self.interval, None)
            if interval is None:
                raise RuntimeError(f"Interval {self.interval} does not exist.")
            if interval.start_time is None or interval.end_time is None:
                self.result = np.full((1, 1, 1), np.nan)
                #raise RuntimeError(f"Interval {interval.name} is not set.") #TODO: raise exception but continue
            else:
                # self.result is the main processing stream, intermediate results are stored there
                self.result = self.manager.data_handler.get_impedance_data(interval.start_time, interval.end_time)

                for preprocessor_info in self.preprocessors:
                    self.result = Preprocessor.initialize_preprocessor(preprocessor_info['name'], **preprocessor_info.get('parameters', {})).apply(interval, self.result)

                base_item = BaseItem.initialize_base_item(self.base_item, **self.parameters)
                self.result = base_item.generate(self.result, interval)
        else:
            # this is a named/composed item
            self.result = self.prerequisites_data[0]

            # create operations and apply in order
            for index, operation_info in enumerate(self.operations):
                # operations have targets (refering to the index of input_data) or None (refering to main result stream)
                self.operations_data[index] = Operation.initialize_operation(operation_info['name'], **operation_info.get('parameters', {}))

                targets = operation_info.get('targets', [None])
                # replace with None object when read as str from json
                targets = [None if item == "None" else item for item in targets]
                # recreate target list with actual data
                targets = [self.result if target is None else self.prerequisites_data[target] for target in targets]
                
                self.result = self.operations_data[index].apply(*targets)

        # after calculation result must not be None or chaining results does not work
        assert self.result is not None

class AnalysisItemManager:
    def __init__(self, data_handler, settings_handler, error_handler):
        self.data_handler = data_handler
        self.settings_handler = settings_handler
        self.error_handler = error_handler
        self.interval_lookup_dict = settings_handler.get_interval_lookup_dict()
        self.analysis_items = {}
        self.exporter = Exporter(settings_handler)
        self.create_analysis_items_from_json_data()

    def add_analysis_item(self, analysis_item):
        # add item using hash as key
        self.analysis_items[analysis_item.hash] = analysis_item
        analysis_item.set_manager(self)

    def get_item_by_name(self, name):
        for item in self.analysis_items.values():
            if item.name == name:
                return item
        else:
            return None

    def create_analysis_items_from_json_data(self):
        data = self.settings_handler.json_data

        for analysis_item_data in data.get('analysis_items', []):
            # assert all names are unique
            analysis_item = AnalysisItem(**analysis_item_data)
            self.add_analysis_item(analysis_item)

    def save_analysis_items_to_json_data(self):
        self.settings_handler.json_data['analysis_items'] = [analysis_item.to_dict() for analysis_item in self.analysis_items.values()]

    def calculate_for_export(self, status_signal):
        # wipe results etc. if already calculated
        for item in self.analysis_items.values():
            item.wipe()
        export_items = {key: value for key, value in self.analysis_items.items() if value.export}
        self.exporter.wipe_exports()
        for export_item in export_items.values():
            status_signal.emit(f"calculating {export_item.name}...")
            _ = export_item.get_results() # results are saved in item
            self.exporter.add_to_export(export_item)
        self.exporter.save_files()

class Exporter:

    def __init__(self, settings_handler):
        self.settings_handler = settings_handler
        self.wipe_exports()

    def wipe_exports(self):
        self.export_queue = []
        self.export_queue.append(("peas_version", "", "", "", "", 0.22))

    def add_to_export(self, item):
        self.export_queue.append((item.name, item.title, item.unit, item.identifier, item.comment, np.squeeze(item.result)))

    def save_files(self):
        def round_to_sig_figs(number, sig_figs=4):
            if np.isnan(number) or number == 0:
                return number
            
            from math import log10, floor
            
            magnitude = floor(log10(abs(number)))
            scale = 10 ** (sig_figs - magnitude - 1)
            
            return round(number * scale) / scale

        output_path = Path(self.settings_handler.get_value("output_path"))
        output_path.mkdir(parents=True, exist_ok=True)

        scalar_items = []
        non_scalar_items = []

        for item_name, title, unit, identifier, comment, data in self.export_queue:
            if np.isscalar(data) or isinstance(data, np.ndarray):
                data = np.vectorize(round_to_sig_figs)(data)

            dim = np.array(data).ndim
            if dim == 0:  # Scalar value
                value = data.item() if isinstance(data, np.ndarray) else data
                scalar_items.append({
                    'item_name': item_name,
                    'title': title,
                    'value': value,
                    'unit': unit,
                    'identifier': identifier,
                    'comment': comment
                })
            elif dim in [1, 2]:  # 1D or 2D array
                filename = f'{item_name}.csv'
                np.savetxt(output_path / filename, data, delimiter=',', fmt='%.4g')
                non_scalar_items.append({
                    'item_name': item_name,
                    'title': title,
                    'value': filename,
                    'unit': unit,
                    'identifier': identifier,
                    'comment': comment
                })
            else:  # 3D or higher dimension
                filename = f'{item_name}.npy'
                np.save(output_path / filename, data)
                non_scalar_items.append({
                    'item_name': item_name,
                    'title': title,
                    'value': filename,
                    'unit': unit,
                    'identifier': identifier,
                    'comment': comment
                })

        all_items = scalar_items + non_scalar_items

        with open(output_path / 'scalars.csv', 'w', newline='') as f:
            fieldnames = ['item_name', 'title', 'value', 'unit', 'identifier', 'comment']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for item in all_items:
                writer.writerow(item)
