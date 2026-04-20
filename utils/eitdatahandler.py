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
import numpy as np
from pathlib import Path
# for file handling
from draeger import load_bin
from sentec import load_zri, load_eit
#for reconstruction
import pyeit.mesh as mesh
import pyeit.eit.protocol as protocol
from pyeit.mesh.shape import thorax
import pyeit.eit.greit as greit
from pyeit.io.daeger_eit import DAEGER_EIT as draeger_eit

class EITDataHandler():
    RECONSTRUCTION_ALGORITHMS = {
            "GREIT": greit
            }
    RECONSTRUCTED_FORMATS = {
        "Draeger reconstructed": ".bin",
        "Sentec reconstructed": ".zri",
        "numpy": ".npz",
    }
    RAW_FORMATS = {
        "Draeger/Sentec raw": '.eit',
        "numpy": '.npz',
    }

    def __init__(self, error_handler):
        self.error_handler = error_handler
        self.wipe_data()

    def _detect_raw_from_contents(self, path: Path) -> bool:
        import numpy as np
        # TODO: replace magic numbers with a more robust detection method
        # 208/928 are the known raw frame sizes for Draeger/Sentec respectively
        with np.load(path) as f:
            return f['data'].shape[1] in (208, 928)

    def is_raw_file(self, path: Path) -> bool:
        suffix = path.suffix.lower()
        in_raw = suffix in self.RAW_FORMATS.values()
        in_reconstructed = suffix in self.RECONSTRUCTED_FORMATS.values()

        if in_raw and not in_reconstructed:
            return True
        if in_reconstructed and not in_raw:
            return False
        return self._detect_raw_from_contents(path)

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
        # TODO: use RECONSTRUCTED_FORMATS
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
        # TODO: use RAW_FORMATS
        def detect_eit_format(path):
            with open(path, "rb") as f:
                b80 = f.read(80)
            if b"---Draeger EIT-Software" in b80:
                return "draeger-eit"

            d4 = b80[:4]
            if len(d4) < 4:
                return "unknown"

            if d4[0:3] == b"\x00\x00\x00" and d4[3] in (2, 3, 4):
                return "lq4pre" if d4[3] == 4 else f"lq{d4[3]}"
            if d4 == b"\x04\x00\x00\x00":
                return "lq4"
            if d4 == b"\x05\x00\x00\x00":
                return "lq5"
            return "unknown"
        if filename.suffix == ".eit":
            try:
                format_str = detect_eit_format(filename)
                if format_str == "draeger-eit":
                    model = draeger_eit(filename)
                    data_raw = model.load()
                    if append:
                        self.data_raw = np.concatenate((self.data_raw, data_raw))
                    else:
                        self.data_raw = data_raw

                    self.timestamps = np.arange(0, self.data_raw.shape[0])*1./model.info['framerate'] # these are not real measured timestamps but evenly spaced points in time
                    return model.info['framerate']
                elif format_str[:2] == "lq":
                    data_raw, trel, _ = load_eit(filename)
                    sampling_frequency = 1.0 / np.median(np.diff(np.array(trel)))
                    if append:
                        self.data_raw = np.concatenate((self.data_raw, data_raw))
                    else:
                        self.data_raw = data_raw
                    self.timestamps = np.arange(0, self.data_raw.shape[0])*1./sampling_frequency # these are not real measured timestamps but evenly spaced points in time
                    return sampling_frequency
                else:
                    raise Exception("Unknown file format")
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

        source_frequency = kwargs.get('source_frequency', None)

        if reconstruction_algorithm == "GREIT":
            reconstruction_settings['greit_p'] = kwargs.get('greit_p', 0.5)
            reconstruction_settings['greit_lambda'] = kwargs.get('greit_lambda', 0.01)
            reconstruction_settings['greit_normalize'] = kwargs.get('greit_normalize', False)
            reference_timestamp = kwargs.get('reconstruction_reference', 0)
            reconstruction_settings['reference_index'] = self.timestamp2index(reference_timestamp)
            if reference_timestamp == 0: # if reference is set to 0, take index of frame with median voltage sum
                sorted_indices = np.argsort(np.sum(self.data_raw))
                reconstruction_settings['reference_index'] = sorted_indices[len(sorted_indices)//2]
                # take minimum instead?
                #refIndex = np.argmin(np.sum(self.data_raw, axis=1))

        lowpass_enabled = kwargs.get('lowpass_enabled', False)
        lowpass_cutoff = kwargs.get('lowpass_cutoff', 8.0)
        lowpass_filter_order = kwargs.get('lowpass_filter_order', 4)

        self.data = self.reconstruct_multi_frame(
            self.data_raw, reconstruction_algorithm,
            source_frequency=source_frequency,
            lowpass_enabled=lowpass_enabled,
            lowpass_cutoff=lowpass_cutoff,
            lowpass_filter_order=lowpass_filter_order,
            **reconstruction_settings
        )

    def reconstruct_multi_frame(self, vx, reconstruction_algorithm, **kwargs):
        """ Reconstruct all images from passed raw data.
        """
        from scipy.signal import butter, sosfiltfilt

        base_conductivity = 1 # hardcoded
        source_frequency = kwargs.get('source_frequency', None)
        lowpass_enabled = kwargs.get('lowpass_enabled', False)
        lowpass_cutoff = kwargs.get('lowpass_cutoff', 8.0)
        lowpass_filter_order = kwargs.get('lowpass_filter_order', 4)

        if lowpass_enabled and source_frequency is not None and source_frequency > 0:
            vx = self._apply_lowpass_filter(vx, source_frequency, lowpass_cutoff, lowpass_filter_order)

        # setup mesh (with thorax shape), protocol
        if vx.shape[1] == 208:
            # Dräger Pulmovista
            n_el = 16
            pair_distance = 1
        elif vx.shape[1] == 928:
            # Sentec BB2
            n_el = 32
            pair_distance = 5
        else:
            n_el = None
            pair_distance = None
            self.error_handler.handle_exception(RuntimeError(f"Unknown raw data structure with length {vx.shape[1]}."))
        mesh_obj = mesh.create(n_el, h0=0.1, fd=thorax)
        protocol_obj = protocol.create(n_el, dist_exc=pair_distance, step_meas=pair_distance, parser_meas="rotate_meas")
        # setup solver
        if reconstruction_algorithm == "GREIT":
            v0 = vx[kwargs.get('reference_index', 0)]
            eit = greit.GREIT(mesh_obj, protocol_obj)
            eit.setup(p=kwargs['greit_p'], lamb=kwargs['greit_lambda'], perm=1, jac_normalized=True)
        else:
            self.error_handler.handle_exception(RuntimeError(f"Solver {reconstruction_algorithm} not found."))
            return None
        conductivityList = []
        for i in range(len(vx)):
            ds = eit.solve(vx[i], v0, normalize=kwargs.get('greit_normalize', False))
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

    def _apply_lowpass_filter(self, vx, source_frequency, cutoff, filter_order):
        from scipy.signal import butter, sosfiltfilt

        try:
            sos = butter(filter_order, cutoff, btype='low', fs=source_frequency, output='sos')
        except (TypeError, ValueError):
            try:
                sos = butter(filter_order, cutoff/(source_frequency/2), btype='low')
            except Exception:
                return vx
        if vx.ndim < 2:
            try:
                return sosfiltfilt(sos, vx)
            except ValueError:
                return vx
        result = np.zeros_like(vx)
        for i in range(vx.shape[1]):
            try:
                result[:, i] = sosfiltfilt(sos, vx[:, i])
            except ValueError:
                result[:, i] = vx[:, i]
        return result
