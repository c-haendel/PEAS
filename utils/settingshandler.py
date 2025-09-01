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
import shutil
from pathlib import Path
from pyqtgraph.parametertree import Parameter

from utils.interval import Interval
from utils.errorhandler import CriticalError
from utils.globalsettings import GlobalSettings

class SettingsHandler():
    def __init__(self, error_handler):

        self.error_handler = error_handler

        # build base parameter
        self.parameter_base_dict = [
            {'name': 'input', 'title': 'input', 'type': 'group', 'children': [
                {'name': 'raw_filename', 'title': 'filename raw', 'type': 'file', 'nameFilter': 'Draeger EIT (*.eit);;numpy (*.npz)'},
                {'name': 'reconstruction', 'title': 'reconstruction', 'type': 'group', 'children': [
                    {'name': 'reconstruction_algorithm', 'title': 'reconstruction algorithm', 'type': 'list', 'limits': ['GREIT'], 'expanded': False, 'children': [
                        {'name': 'reconstruction_reference', 'title': 'reconstruction reference', 'type': 'float', 'value': 0, 'suffix': 's'},
                        {'name': 'greit_p', 'title': 'GREIT p', 'type': 'float', 'value': 0.5},
                        {'name': 'greit_lambda', 'title': 'GREIT lambda', 'type': 'float', 'value': 0.01},
                    ]},
                    {'name': 'reconstruct', 'title': 'reconstruct', 'type': 'action'}
                ]},
                {'name': 'reconstructed_filename', 'title': 'filename reconstructed', 'type': 'file', 'nameFilter': 'Draeger binary (*.bin);;Sentec zri (*.zri);;numpy (*.npz)'},
                {'name': 'source_frequency', 'title': 'sampling frequency', 'type': 'float', 'suffix': 'Hz'},
                {'name': 'analysis', 'title': 'analysis', 'type': 'group', 'children': [
                    {'name': 'analysis_template', 'title': 'analysis template', 'type': 'file', 'nameFilter': 'analysis template (*.json)', 'value': GlobalSettings.ANALYSIS_ITEMS_FILE},
                    {'name': 'edit_analyses', 'title': 'edit analyses', 'type': 'action'},
                ]},
            ]},
            {'name': 'intervals', 'title': 'intervals', 'type': 'group', 'children': []},
            {'name': 'output', 'type': 'group', 'children': [
                {'name': 'output_path', 'title': 'output path', 'type': 'file', 'fileMode': 'Directory', 'options': 'ShowDirsOnly'},
                {'name': 'export_results', 'title': 'export results', 'type': 'action'},
            ]},
        ]
        self.json_data = {}
        self.parameter = Parameter.create(name='params', type='group', children=self.parameter_base_dict)

        self.read_analysis_template()
        self.create_intervals_from_json_data()
        self.rebuild_intervals_parameter()

    def get_interval_lookup_dict(self):
        return {i.name: i for i in self.interval_list}

    def export_analysis_json(self):
        """
        Read analysis template file and write it to output, assuming file is always up-to-date.
        """
        shutil.copy2(self.get_value('analysis_template'), Path(self.get_value('output_path')) / "analysis_template.json")

    def read_analysis_template(self):
        def replace_codes(d):
            """
            Recursively replaces 'inf' and '-inf' strings with np.inf and -np.inf in a nested dictionary.
            """
            if isinstance(d, dict):
                return {key: replace_codes(value) for key, value in d.items()}
            elif isinstance(d, list):
                return [replace_codes(item) for item in d]
            elif d == "inf":
                return np.inf
            elif d == "-inf":
                return -np.inf
            elif d == "None":
                return None
            else:
                return d
        try:
            with open(self.get_value('analysis_template'), 'r') as fp:
                self.json_data = json.load(fp)
        except Exception as e:
            raise CriticalError(f"Failed to read the JSON file: {e}")
        self.json_data = replace_codes(self.json_data)

    def create_intervals_from_json_data(self):
        # dynamic elements: breath detector selection; breath detector settings
        # create list of interval-objects from interval_list
        self.interval_list = []
        for interval_dict in self.json_data.get("intervals", []):
            self.interval_list.append(Interval.from_dict(self, interval_dict))

    def save_intervals_to_json_data(self):
        self.json_data['intervals'] = [interval.to_dict() for interval in self.interval_list]

    def write_analysis_template(self):
        def replace_codes(d):
            """
            Recursively replaces np.inf, -np.inf, and None with strings for JSON serialization.
            """
            if isinstance(d, dict):
                return {key: replace_codes(value) for key, value in d.items()}
            elif isinstance(d, list):
                return [replace_codes(item) for item in d]
            elif d == np.inf:
                return "inf"
            elif d == -np.inf:
                return "-inf"
            elif d is None:
                return "None"
            else:
                return d

        try:
            with open(self.get_value('analysis_template'), 'w') as fp:
                json.dump(replace_codes(self.json_data), fp, indent=4)
        except Exception as e:
            raise CriticalError(f"Failed to save JSON file: {e}")

    def write_state_file(self):
        # TODO: state file must include reconstruction settings
        def convert_float32_to_float(obj):
            if isinstance(obj, dict):
                return {k: convert_float32_to_float(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_float32_to_float(i) for i in obj]
            elif isinstance(obj, np.float32):
                return float(obj)
            return obj
        # Create a Path object for the output directory
        output_path = Path(self.get_value("output_path"))
        # Create the output directory if it does not exist
        output_path.mkdir(parents=True, exist_ok=True)

        # save all modifyable interval settings:
        # * interval times
        # * detector + detector settings
        state_list = []
        for interval in self.interval_list:
            state_list.append(interval.to_state_dict())
        with open(output_path / GlobalSettings.STATE_FILE, 'w') as file:
            json.dump(convert_float32_to_float(state_list), file)

    def read_state_file(self):
        # Create a Path object for the output directory
        path = Path(self.get_value("output_path")) / GlobalSettings.STATE_FILE
        for interval in self.interval_list:
            interval.flags = []
            self.start_time = None
            self.end_time = None
        if path.is_file():
            try:
                # Open and read the JSON file
                with open(path, 'r') as file:
                    state_data = json.load(file)
                for interval in self.interval_list:
                    for state_interval in state_data:
                        if interval.name == state_interval.get('name', None):
                            interval.update_from_state_dict(state_interval)
                self.sync_tree_to_interval() # read settings to be displayed
            except Exception as e:
                self.error_handler.handle_exception(e)
                return None
        else:
            return None
    
    def rebuild_intervals_parameter(self):
        # create intervals parameter
        intervals_parsed = []
        for interval in self.interval_list:
            intervals_parsed.append(interval.to_parameter_dict())
        intervals_parsed.append({'name': 'edit_intervals', 'title': 'edit intervals', 'type': 'action'})
        new_intervals_parameter = Parameter.create(name='intervals', type='group', children=intervals_parsed)
        # replace parameter in parameter tree
        old_intervals_parameter = self.parameter.param('intervals')
        old_intervals_parameter.clearChildren()
        for child in new_intervals_parameter.children():
            child_config = child.saveState()
            old_intervals_parameter.addChild(child_config)

    def param_recursive(self, target_name, root_param=None):
        if root_param is None:
            root_param = self.parameter
        if root_param.name() == target_name:
            return root_param

        for child in root_param.children():
            result = self.param_recursive(target_name, child)
            if result is not None:
                return result
        return None

    def sync_interval_to_tree(self):
        # settings of intervals and detectors are stored in the respective objects, not the parameter tree. Manual sync after change in tree.
        for interval in self.interval_list:
            for key, _ in interval.detector.settings.items():
                interval.detector.settings[key]['value'] = self.get_value(f"{interval.name}_detector_{key}")

    def sync_tree_to_interval(self):
        # settings of intervals and detectors are stored in the respective objects, not the parameter tree. Manual sync after change in interval settings.
        for interval in self.interval_list:
            for key, setting in interval.detector.settings.items():
                self.set_value(f"{interval.name}_detector_{key}", setting['value'])

    def set_value(self, name, value):
        param = self.param_recursive(name)
        if param is not None:
            param.setValue(value)

    def set_append_state(self, state):
        """ Show/hide append button if file is loaded
            0: hide all, 1: show reconstructed append, 2: show raw append
        """
        # hide all
        parent = self.parameter.param('input')
        for param_name in ['reconstructed_filename_append', 'raw_filename_append']:
            try:
                parent.removeChild(parent.child(param_name))
            except KeyError:
                pass
        if state == 1:
            # show reconstructed append
            parent.insertChild(3, {'name': 'reconstructed_filename_append', 'title': 'append filename reconstructed', 'type': 'file', 'nameFilter': 'Draeger binary (*.bin);;Sentec zri (*.zri)'})
        elif state == 2:
            # show raw append
            parent.insertChild(1, {'name': 'raw_filename_append', 'title': 'append filename raw', 'type': 'file', 'nameFilter': 'Draeger EIT (*.eit)'})

    def get_value(self, name):
        param = self.param_recursive(name)
        if param is not None:
            return param.value()
        else:
            return None

    def get_param_dict(self, param):
        if param is None:
            return {}

        result = {}
        result[param.name()] = param.value()
        
        if param.hasChildren():
            for child in param.children():
                result.update(self.get_param_dict(child))
        
        return result
