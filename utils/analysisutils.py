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

from typing import Union

class FunctionType:
    POINT = "point"
    LINE = "line"

class Function:
    def __init__(self, function_type: str, **kwargs):
        self.function_type = function_type
        self.params = kwargs

    def evaluate(self, x: float) -> Union[float, None]:
        if self.function_type == FunctionType.POINT:
            if x == self.params['x']:
                return self.params['y']
            else:
                return None

        elif self.function_type == FunctionType.LINE:
            if 'slope' in self.params:
                return self.params['slope'] * x + self.params['intercept']
            else:
                return None if x != self.params['x'] else float('inf')

    def to_dict(self):
        return {'function_type': self.function_type, 'params': self.params}

    @classmethod
    def from_dict(cls, input_dict):
        return cls(input_dict['function_type'], **input_dict['params'])

class Flag:
    def __init__(self, flag_type, function: Function):
        # eoi: end of inspiration
        # eofi: end of forced inspiration
        self.flag_type = flag_type # eoi_time, eoe_time, eoi_z, eoe_z, peak_flow, other
        self.function = function
        self.plot_item = None

    def __del__(self):
        if self.plot_item is not None and self.plot_item.scene() is not None:  # Check if line is still part of a scene
            self.plot_item.scene().removeItem(self.plot_item)

    @classmethod
    def from_dict(cls, input_dict):
        function = Function.from_dict(input_dict['function'])
        return cls(input_dict['flag_type'], function)

    def to_dict(self):
        return {'flag_type': self.flag_type, 'function': self.function.to_dict()}

    def export(self):
        # Logic to export analysis_items
        pass
