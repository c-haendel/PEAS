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

# for data manipulation
import math
import numpy as np
# for file handling
from draeger import load_bin
from sentec import load_zri, load_eit
#for reconstruction
import pyeit.mesh as mesh
import pyeit.eit.protocol as protocol
from pyeit.mesh.shape import thorax
import pyeit.eit.greit as greit
from pyeit.io.daeger_eit import DAEGER_EIT as draeger_eit

AVAILABLE_RECONSTRUCTION_ALGORITHMS = {
        "GREIT": greit
        }

class EITDataHandler():
    def __init__(self, error_handler):
        self.error_handler = error_handler
        self.wipe_data()

    def timestamp2index(self, timestamp):
        """ return first index after timestamp
        """
        after = np.squeeze(np.argwhere(self.timestamps >= timestamp))
        if after.size<=1: return len(self.timestamps)-1 # somehow 0 and 1 can be empty
        else: return after[0]

    def wipe_data(self):
        self.data = np.array([])
        self.data_raw = np.array([])
        self.timestamps = []

    def load_reconstructed_input_file(self, filename, append=False):
        if filename.suffix == ".bin":
            try:
                data, _, _, sampling_frequency = load_bin(filename, max_channels_medibus=0)
                if append:
                    self.data = np.concatenate((self.data, data))
                else:
                    self.data = data
                self.timestamps = np.arange(0, self.data.shape[0])*1./sampling_frequency # these are not real measured timestamps but evenly spaced points in time
                return sampling_frequency
            except Exception as e:
                self.error_handler.handle_exception(e)
                return 0
        elif filename.suffix == ".zri":
            try:
                data, timestamps, _, _ = load_zri(filename)
                sampling_frequency = len(data) / (timestamps[-1] - timestamps[0]).total_seconds()
                if append:
                    self.data = np.concatenate((self.data, data))
                else:
                    self.data = data
                self.timestamps = np.arange(0, self.data.shape[0])*1./sampling_frequency # these are not real measured timestamps but evenly spaced points in time
                return sampling_frequency
            except Exception as e:
                self.error_handler.handle_exception(e)
                return 0
        elif filename.suffix == ".npz":
            try:
                with np.load(filename) as f:
                    self.data = f["data"]
                    sampling_frequency = f["framerate"].item()
                self.timestamps = np.arange(0, self.data.shape[0])*1./sampling_frequency # these are not real measured timestamps but evenly spaced points in time
                return sampling_frequency
            except Exception as e:
                self.error_handler.handle_exception(e)
                return 0
        else:
            try:
                raise RuntimeError(f"Unknown filetype {filename.suffix}")
            except Exception as e:
                self.error_handler.handle_exception(e)
                return 0

    def load_raw_input_file(self, filename, append=False):
        if filename.suffix == ".eit":
            try:
                model = draeger_eit(filename)
                data_raw = model.load()
                if append:
                    self.data_raw = np.concatenate((self.data_raw, data_raw))
                else:
                    self.data_raw = data_raw
                self.timestamps = np.arange(0, self.data_raw.shape[0])*1./model.info['framerate'] # these are not real measured timestamps but evenly spaced points in time
                return model.info['framerate']
            except:
                try:
                    data_raw, timestamps, _, _ = load_eit(filename)
                    sampling_frequency = len(data_raw) / (timestamps[-1] - timestamps[0]).total_seconds()
                    if append:
                        self.data_raw = np.concatenate((self.data_raw, data_raw))
                    else:
                        self.data_raw = data_raw
                    self.timestamps = np.arange(0, self.data.shape[0])*1./sampling_frequency # these are not real measured timestamps but evenly spaced points in time
                    return sampling_frequency
                except Exception as e:
                    self.error_handler.handle_exception(e)
                    return 0
        elif filename.suffix == ".npz":
            try:
                with np.load(filename) as f:
                    self.data_raw = f["data"]
                    self.timestamps = np.arange(0, self.data_raw.shape[0])*1./f['framerate'] # these are not real measured timestamps but evenly spaced points in time
                    return f["framerate"].item()
            except Exception as e:
                self.error_handler.handle_exception(e)
                return 0
        else:
            try:
                raise RuntimeError(f"Unknown filetype {filename.suffix}")
            except Exception as e:
                self.error_handler.handle_exception(e)
                return 0

    def reconstruct(self, **kwargs):
        """ Wrapper for reconstruct_multi_frame. Reconstruct all images from raw data within the datahandler.
        """

        reconstruction_algorithm = kwargs.get('reconstruction_algorithm', None)
        reconstruction_settings = {}

        if reconstruction_algorithm == "GREIT":
            reconstruction_settings['greit_p'] = kwargs.get('greit_p', 0.5)
            reconstruction_settings['greit_lambda'] = kwargs.get('greit_lambda', 0.01)
            reference_timestamp = kwargs.get('reconstruction_reference', 0)
            reconstruction_settings['reference_index'] = self.timestamp2index(reference_timestamp)
            if reference_timestamp == 0: # if reference is set to 0, take index of frame with median voltage sum
                sorted_indices = np.argsort(np.sum(self.data_raw))
                reconstruction_settings['reference_index'] = sorted_indices[len(sorted_indices)//2]
                # take minimum instead?
                #refIndex = np.argmin(np.sum(self.data_raw, axis=1))
        else:
            raise RuntimeError(f"Unknown reconstruction algorithm {reconstruction_algorithm}")

        self.data = self.reconstruct_multi_frame(self.data_raw, reconstruction_algorithm, **reconstruction_settings)

    def reconstruct_multi_frame(self, vx, reconstruction_algorithm, **kwargs):
        """ Reconstruct all images from passed raw data.
        """
        base_conductivity = 1 # hardcoded
        # setup mesh (with thorax shape), protocol
        n_el = 1 + math.sqrt(1 + vx[1])
        mesh_obj = mesh.create(n_el, h0=0.1, fd=thorax)
        protocol_obj = protocol.create(n_el, dist_exc=1, step_meas=1, parser_meas="std") # TODO: verify settings for sentec *.eit
        # setup solver
        if reconstruction_algorithm == "GREIT":
            v0 = vx[kwargs.get('reference_index', 0)]
            eit = greit.GREIT(mesh_obj, protocol_obj) # default grid size 32*32
            eit.setup(p=kwargs['greit_p'], lamb=kwargs['greit_lambda'], perm=1, jac_normalized=True)
        else:
            self.error_handler.handle_exception(RuntimeError(f"Solver {reconstruction_algorithm} not found."))
            return None

        conductivityList = []
        for i in range(len(vx)):
            ds = eit.solve(vx[i], v0, normalize=True)
            _, _, ds = eit.mask_value(ds, mask_value=np.nan)
            conductivityList.append(np.real(ds))
        impedance = 1/(np.asarray(conductivityList)+base_conductivity)
        # Correction for inverted v/d-axis
        impedance = np.flip(impedance, axis=1)
        return impedance

    def get_impedance_data(self, start_timestamp, end_timestamp=None):
        if self.data.size == 0:
            return None, None
        start_idx = self.timestamp2index(start_timestamp)
        if end_timestamp is None:
            return self.data[start_idx,:,:], self.timestamps[start_idx]
        else:
            end_idx = self.timestamp2index(end_timestamp)
            return self.data[start_idx:end_idx,:,:], self.timestamps[start_idx:end_idx]

    def get_plot_data(self, voltage_mode=False):
        if voltage_mode:
            # return raw voltage data
            if self.data_raw.size == 0:
                return None, None
            plotData = np.nansum(self.data_raw, axis=(1))
        else:
            # return reconstructed data
            if self.data.size == 0:
                return None, None
            plotData = np.nansum(self.data, axis=(1,2))
        return plotData, self.timestamps
