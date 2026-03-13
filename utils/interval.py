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

from utils.detector import Detector, DETECTORS
from utils.analysisutils import Flag

class Interval():
    def __init__(self, settings_handler, name, title, start_time, end_time, detector_name, screenshot):
        self.name = name
        self.settings_handler = settings_handler
        self.title = title
        self.start_time = start_time
        self.end_time = end_time
        self.detector = Detector.initialize_detector(detector_name, self.settings_handler)
        self.flags = []
        self.screenshot = screenshot

    @classmethod
    def from_dict(cls, settings_handler, input_dict):
        return cls(settings_handler, input_dict['name'], input_dict['title'], input_dict.get('start_time', None), input_dict.get('end_time', None), input_dict['detector'], input_dict.get('screenshot', False))

    def to_dict(self):
        return {
                "name": self.name,
                "title": self.title,
                "detector": self.detector.name,
                "screenshot": self.screenshot
        }

    def run_detector(self, data):
        # data is complete signal in the interval
        if self.detector is not None and data is not None and len(data)>1:
            self.flags = self.detector.detect(data, self.flags)

    def to_parameter_dict(self):
        # Initialize the 'detector' dictionary with its own 'children' list
        detector_dict = {
            'name': f"{self.name}_detector",
            'title': 'detector',
            'type': 'list',
            'limits': list(DETECTORS.keys()),
            'value': self.detector.name,
            'expanded': False,
            'children': []
        }

        # Add additional settings from the detector as children of the 'detector' dictionary
        for key, meta in self.detector.settings.items():
            value = meta['value']
            dtype = meta.get('dtype', float)  # Default type is float if not specified
            unit = meta.get('unit', '')  # Default unit is empty if not specified
            readonly = meta.get('readonly', False)  # Default is modifyable if not specified

            setting_dict = {
                'name': f"{self.name}_detector_{key}",
                'title': key,
                'type': dtype.__name__,  # Convert the Python data type to its string name
                'value': value,
                'suffix': unit,
                'readonly': readonly,
            }
            detector_dict['children'].append(setting_dict)

        # Create the main parameter dictionary with the 'detector' dictionary as one of its children
        interval_parsed = {
            'name': self.name,
            'title': self.title,
            'type': 'group',
            'children': [detector_dict]
        }

        return interval_parsed

    def to_state_dict(self):
        state_dict = {
            'name': self.name,
            'start_time': self.start_time,
            'end_time': self.end_time,
            'detector': self.detector.name,
            'detector_settings': self.detector.settings_to_state_dict(),
            'flags': [flag.to_dict() for flag in self.flags]
        }
        return state_dict

    def update_from_state_dict(self, state_dict):
        # name is already identical
        self.start_time = state_dict.get('start_time', self.start_time)
        self.end_time = state_dict.get('end_time', self.end_time)
        # detector is replaced no matter what
        if state_dict.get('detector', None) is not None:
            self.detector = Detector.initialize_detector(state_dict['detector'], self.settings_handler)
            self.detector.update_from_state_dict(state_dict.get('detector_settings', {}))
        self.flags = [Flag.from_dict(flag_dict) for flag_dict in state_dict.get('flags', [])]
