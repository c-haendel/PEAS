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
from scipy.interpolate import CubicSpline
from scipy.signal import butter, sosfiltfilt, find_peaks
from typing import List, Tuple

from utils.analysisutils import Flag, Function, FunctionType

class Detector():
    def __init__(self, settings_handler, name, detector_type):
        self.name = name
        self.settings_handler = settings_handler
        self.detector_type = detector_type
        self.settings = {}

    def get_setting(self, key):
        return self.settings.get(key, {}).get('value', None)
    def set_setting(self, key, new_value):
        if key not in self.settings:
            return

        expected_type = self.settings[key]['dtype']
        
        if isinstance(new_value, expected_type):
            self.settings[key]['value'] = new_value
        else:
            self.settings[key]['value'] = expected_type(new_value)
    def get_setting_metadata(self, key):
        if key in self.settings:
            return {k: self.settings[key][k] for k in ('dtype', 'unit', 'readonly')}
        return None
    def settings_to_state_dict(self):
        return {key: value['value'] for key, value in self.settings.items()}
    def update_from_state_dict(self, state_dict):
        for setting_key, _ in self.settings.items():
            for state_setting_key, state_setting_value in state_dict.items():
                if setting_key == state_setting_key:
                    self.set_setting(setting_key, state_setting_value)
    def detect(self, data, flags):
        # dummy function overwritten by subclasses
        return []

    @staticmethod
    def initialize_detector(detector_name, settings_handler, **kwargs):
        if detector_name in AVAILABLE_DETECTORS:
            detector_class = AVAILABLE_DETECTORS[detector_name]
            return detector_class(settings_handler, detector_name, detector_class.DETECTOR_TYPE, **kwargs)
        else:
            raise ValueError(f"Invalid detector name: {detector_name}")

class MinMaxMinDetector(Detector):
    """ Simplest detector, taking the global maximum as well as the global minima before and after
    """
    DETECTOR_TYPE = "maneuver"

    def __init__(self, settings_handler, name, detector_type, **kwargs):
        super().__init__(settings_handler, name, detector_type, **kwargs)

    def detect(self, data: tuple, flags: List[Flag]) -> List[Flag]:

        if data[0] is None or len(data[0]) == 0: return []

        flags = []
        original_data, original_time_array = data
        summed_data = np.nansum(original_data, axis=(1, 2))

        max_idx = np.argmax(summed_data)
        max_val = summed_data[max_idx]
        max_time = original_time_array[max_idx]

        flags.append(Flag("eoi_time", Function(FunctionType.LINE, x=max_time)))  # Vertical line
        flags.append(Flag("eoi_z", Function(FunctionType.LINE, slope=0, intercept=max_val)))  # Horizontal line

        min_idx_before_max = np.argmin(summed_data[:max_idx])
        min_time_before_max = original_time_array[min_idx_before_max]

        flags.append(Flag("eoe_time", Function(FunctionType.LINE, x=min_time_before_max)))  # Vertical line

        min_idx_after_max = max_idx + np.argmin(summed_data[max_idx:])
        min_val = summed_data[min_idx_after_max]
        min_time = original_time_array[min_idx_after_max]

        # Create Flags for minimum after maximum
        flags.append(Flag("eoe_time", Function(FunctionType.LINE, x=min_time)))  # Vertical line
        flags.append(Flag("eoe_z", Function(FunctionType.LINE, slope=0, intercept=min_val)))  # Horizontal line

        return flags

class SimplePeakFlowDetector(Detector):
    """ Detector for forced expiration.
Create horizontal flag "eofi_z" at global maximum. Create a vertical flag "eoe_time" at the global minimum before the maximum. Create a vertical flag "eofe_time" and a horizontal flag "eofe_z" at the minimum after the maximum.
Perform cubic spline interpolation for the signal and resample according to sliding window. Using a centered sliding window of a specified size, detect the time where the slope is lowest (meaning the steepest drop). Create a point flag "peak_expiratory_flow" at that point. Also create a flag "peak_expiratory_flow_tangent" that goes through this point and has the determined slope. Finally, find the intersection of this tangent and the horizontal line at the maximum. At this point, place a vertical flag "eofi_time".
    """
    DETECTOR_TYPE = "maneuver"

    def __init__(self, settings_handler, name, detector_type, **kwargs):
        super().__init__(settings_handler, name, detector_type, **kwargs)
        self.settings = {
            'sliding_window_size': {
                'value': kwargs.get('sliding_window_size', 0.01),
                'dtype': float,
                'unit': 's',
            },
            'back_extrapolated_fraction': {
                'value': 0,
                'dtype': float,
                'unit': '',
                'readonly': True,
            }
        }

    def detect(self, data: Tuple[np.ndarray, List[float]], flags: List[Flag]) -> List[Flag]:

        if data[0] is None or len(data[0]) == 0: return []

        flags = []
        #original_sampling_rate = self.settings_handler.get_value("source_frequency")
        original_data, original_time_array = data
        window_size = self.get_setting('sliding_window_size')

        summed_data = np.nansum(original_data, axis=(1, 2))

        if len(summed_data) <= 2:
            return []
        cs = CubicSpline(original_time_array, summed_data)

        # Define the resampled frequency and time array based on the sliding window size
        # sliding window has width of 2 sampling intervals
        resampled_sampling_rate = 1 / window_size * 2
        resampled_time_array = np.linspace(min(original_time_array), max(original_time_array), int(len(original_time_array) * resampled_sampling_rate / len(original_time_array)))

        resampled_data = cs(resampled_time_array)

        max_idx = np.argmax(resampled_data)
        max_val = resampled_data[max_idx]

        flags.append(Flag("eofi_z", Function(FunctionType.LINE, slope=0, intercept=max_val)))  # Horizontal line

        min_idx_before_max = np.argmin(resampled_data[:max_idx])
        min_time_before_max = (min_idx_before_max / resampled_sampling_rate)
        min_time_before_max = resampled_time_array[min_idx_before_max]

        flags.append(Flag("eoe_time", Function(FunctionType.LINE, x=min_time_before_max)))  # Vertical line

        min_idx_after_max = max_idx + np.argmin(resampled_data[max_idx:])
        min_val = resampled_data[min_idx_after_max]
        min_time = (min_idx_after_max / resampled_sampling_rate)
        min_time = resampled_time_array[min_idx_after_max]

        flags.append(Flag("eofe_time", Function(FunctionType.LINE, x=min_time)))  # Vertical line
        flags.append(Flag("eofe_z", Function(FunctionType.LINE, slope=0, intercept=min_val)))  # Horizontal line

        slopes = np.gradient(resampled_data, resampled_time_array)

        # Find the point of lowest slope (steepest drop)
        lowest_slope_idx = np.argmin(slopes)
        lowest_slope_time = (lowest_slope_idx / resampled_sampling_rate)
        lowest_slope_time = resampled_time_array[lowest_slope_idx]
        lowest_slope_value = resampled_data[lowest_slope_idx]

        # Create Flags for the point of lowest slope and its tangent
        flags.append(Flag("peak_expiratory_flow", Function(FunctionType.POINT, x=lowest_slope_time, y=lowest_slope_value)))
        flags.append(Flag("peak_expiratory_flow_tangent", Function(FunctionType.LINE, slope=slopes[lowest_slope_idx], intercept=-lowest_slope_time * slopes[lowest_slope_idx] + lowest_slope_value)))

        # Find intersection of tangent and horizontal line at maximum
        intersection_time = (max_val - lowest_slope_value) / slopes[lowest_slope_idx] + lowest_slope_time
        flags.append(Flag("eofi_time", Function(FunctionType.LINE, x=intersection_time)))  # Vertical line

        # Calculate quality measure 'pre-expired fraction'
        # Find the nearest data point to intersection_time
        nearest_idx_to_intersection = np.argmin(np.abs(resampled_time_array - intersection_time))
        value_at_intersection = resampled_data[nearest_idx_to_intersection]

        # Calculate the fraction dropped at intersection_time
        self.set_setting('back_extrapolated_fraction', (value_at_intersection - max_val) / (min_val - max_val))

        return flags

class ReverseCooldownDetector(Detector):
    """The detect function takes the summed signal. It determines all minima and maxima in the signal. It then eliminates all maxima where another maximum exists within the time specified by "cooldown". The same is done for minima. Finally, minima and maxima must be alternating, meaning if there are multiple maxima between two minima, all but the last maximum must be eleminated. Minima analogously.
With this list of minima and maxima, create Flag objects with vertical lines named "eoe_time" for minima and "eoi_time" for maxima. Append these to the flags list and return it.
    """
    DETECTOR_TYPE = "breath"

    def __init__(self, settings_handler, name, detector_type, **kwargs):
        super().__init__(settings_handler, name, detector_type, **kwargs)
        self.settings = {
            'cooldown': {
                'value': kwargs.get('cooldown', 4),
                'dtype': float,
                'unit': 's',
            },
        }

    def get_setting(self, key):
        return self.settings[key]['value']

    def detect(self, data: Tuple[np.ndarray, List[float]], flags: List[Flag]) -> List[Flag]:

        if data[0] is None or len(data[0]) == 0: return []

        flags = []
        original_data, original_time_array = data

        # Sum the data over the x and y dimensions to get a 1D array
        summed_data = np.nansum(original_data, axis=(1, 2))

        # Find all minima and maxima
        minima_indices = (np.diff(np.sign(np.diff(summed_data))) > 0).nonzero()[0] + 1
        maxima_indices = (np.diff(np.sign(np.diff(summed_data))) < 0).nonzero()[0] + 1

        # Get the cooldown setting
        cooldown = self.get_setting('cooldown')
        cooldown_idx = int(cooldown * len(original_time_array) / (original_time_array[-1] - original_time_array[0]))

        # Eliminate peaks within cooldown
        def filter_within_cooldown(indices):
            filtered_indices = []
            last_valid_idx = -cooldown_idx - 1
            for idx in indices:
                if idx - last_valid_idx > cooldown_idx:
                    last_valid_idx = idx
                    filtered_indices.append(idx)
            return np.array(filtered_indices)

        minima_indices = filter_within_cooldown(minima_indices)
        maxima_indices = filter_within_cooldown(maxima_indices)

        # Ensure minima and maxima are alternating
        all_indices = np.concatenate((minima_indices, maxima_indices))
        all_labels = np.concatenate((['min']*len(minima_indices), ['max']*len(maxima_indices)))

        # Sort indices and labels by time
        sort_order = np.argsort(all_indices)
        all_indices = all_indices[sort_order]
        all_labels = all_labels[sort_order]

        # Filter to ensure alternation
        filtered_indices = []
        last_label = None
        for idx, label in zip(all_indices, all_labels):
            if last_label is None or last_label != label:
                filtered_indices.append(idx)
                last_label = label

        # Extract the final lists of minima and maxima indices
        minima_indices = np.array([idx for idx in filtered_indices if idx in minima_indices])
        maxima_indices = np.array([idx for idx in filtered_indices if idx in maxima_indices])

        # Ignore the first maximum if it comes before the first minimum
        # Ignore the last minimum if it comes after the last maximum
        if len(minima_indices) > 0 and len(maxima_indices) > 0:
            if maxima_indices[0] < minima_indices[0]:
                maxima_indices = maxima_indices[1:]
            if minima_indices[-1] > maxima_indices[-1]:
                minima_indices = minima_indices[:-1]

        # Create flags for minima and maxima
        for idx in minima_indices:
            min_time = original_time_array[idx]
            flags.append(Flag("eoe_time", Function(FunctionType.LINE, x=min_time)))  # Vertical line

        for idx in maxima_indices:
            max_time = original_time_array[idx]
            flags.append(Flag("eoi_time", Function(FunctionType.LINE, x=max_time)))  # Vertical line

        return flags

class ManualBreathDetector(Detector):
    """ Pseudo-detector for regular breathing.
    """
    DETECTOR_TYPE = "breath"

    def __init__(self, settings_handler, name, detector_type, **kwargs):
        super().__init__(settings_handler, name, detector_type, **kwargs)
        self.settings = {
            'breaths': {
                'value': kwargs.get('breaths', 3),
                'dtype': int,
                'unit': '',
            },
        }
    def detect(self, data, flags):

        if data[0] is None or len(data[0]) == 0: return []

        time_step = len(data[1]) / self.settings_handler.get_value("source_frequency") / (self.get_setting('breaths')*2 + 1)

        flag_list = []
        for i in range(1, self.get_setting('breaths')*2 + 1):
            flag_list.append(Flag("eoe_time" if i%2==0 else "eoi_time", Function(FunctionType.LINE, x=i*time_step + data[1][0])))

        return flag_list

class LowPassDetector(Detector):
    DETECTOR_TYPE = "breath"

    def __init__(self, settings_handler, name, detector_type, **kwargs):
        super().__init__(settings_handler, name, detector_type, **kwargs)
        self.settings = {
            'filter_order': {
                'value': kwargs.get('filter_order', 4),
                'dtype': int,
                'unit': '',
            },
            'cutoff': {
                'value': kwargs.get('cutoff', 1.0),
                'dtype': float,
                'unit': 'Hz',
            },
            'timeout': {
                'value': kwargs.get('timeout', 2.0),
                'dtype': float,
                'unit': 's',
            }
        }

    def detect(self, data: Tuple[np.ndarray, List[float]], flags: List[Flag]) -> List[Flag]:

        if data[0] is None or len(data[0]) == 0: return []

        original_data, original_time_array = data
        flags = []

        if len(original_data) <= 1 or len(original_time_array) <= 1:
            return flags

        # Sum the data over the x and y dimensions to get a 1D array
        summed_data = np.nansum(original_data, axis=(1, 2))

        # Low-pass filter settings
        filter_order = self.get_setting('filter_order')
        cutoff = self.get_setting('cutoff')
        fs = 1 / (original_time_array[1] - original_time_array[0])  # Sampling frequency
        timeout = self.get_setting('timeout')

        # Apply Butterworth filter
        sos = butter(filter_order, cutoff, fs=fs, output='sos')
        try:
            filtered_data = sosfiltfilt(sos, summed_data)
        except ValueError:
            return flags

        # Find minima and maxima
        minima_indices, _ = find_peaks(-filtered_data)
        maxima_indices, _ = find_peaks(filtered_data)

        # Filter out minima and maxima that are too close (less than 'timeout' apart)
        for i in reversed(range(len(minima_indices) - 1)):
            if original_time_array[minima_indices[i + 1]] - original_time_array[minima_indices[i]] < timeout:
                max_between = [x for x in maxima_indices if minima_indices[i] < x < minima_indices[i + 1]]
                if max_between:
                    maxima_indices = np.delete(maxima_indices, np.where(maxima_indices == max_between[0]))
                minima_indices = np.delete(minima_indices, i)

        for i in reversed(range(len(maxima_indices) - 1)):
            if original_time_array[maxima_indices[i + 1]] - original_time_array[maxima_indices[i]] < timeout:
                min_between = [x for x in minima_indices if maxima_indices[i] < x < maxima_indices[i + 1]]
                if min_between:
                    minima_indices = np.delete(minima_indices, np.where(minima_indices == min_between[0]))
                maxima_indices = np.delete(maxima_indices, i)

        # Ensure the first is a minimum and the last is a maximum
        if minima_indices.size == 0 or maxima_indices.size == 0 or minima_indices[0] > maxima_indices[-1]:
            minima_indices = np.array([])
            maxima_indices = np.array([])
        while minima_indices.size > 0 and maxima_indices.size > 0 and maxima_indices[0] < minima_indices[0]:
            maxima_indices = maxima_indices[1:]
        while minima_indices.size > 0 and maxima_indices.size > 0 and minima_indices[-1] > maxima_indices[-1]:
            minima_indices = minima_indices[:-1]

        # Create Flag objects
        for idx in minima_indices:
            flags.append(Flag("eoe_time", Function(FunctionType.LINE, x=original_time_array[idx])))
        for idx in maxima_indices:
            flags.append(Flag("eoi_time", Function(FunctionType.LINE, x=original_time_array[idx])))

        return flags

AVAILABLE_DETECTORS = {
        #"reverse_cooldown": ReverseCooldownDetector,
        "low_pass": LowPassDetector,
        "min_max_min": MinMaxMinDetector,
        #"manual_breath": ManualBreathDetector,
        "simple_peak_flow": SimplePeakFlowDetector
        }
