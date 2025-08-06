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

class GlobalSettings:
    PLOT_COLOR_DATA_Z = (255, 255, 255) # white
    PLOT_COLOR_DATA_V = (180, 180, 180) # gray
    PLOT_COLOR_REGION_PEN = (200, 50, 50)
    PLOT_COLOR_REGION_BRUSH = (200, 50, 50, 50)
    PLOT_COLOR_REGION_HOVERBRUSH = (200, 50, 50, 80)
    PLOT_COLOR_EOI = (50, 50, 200, 255) # blue
    PLOT_COLOR_EOE = (125, 50, 125, 255) # purple
    PLOT_COLOR_FLAG_DEFAULT = (255, 255, 255, 100)
    PLOT_Z_DATA = 5
    PLOT_Z_REGION = 1
    PLOT_Z_FLAG = 3
    PLOT_SCREENSHOT_WIDTH = 1000

    PLOT_PG_COLORMAP = "CET-C3"

    EOE_FLAG_TYPES = ["eoe_time", "eofe_time"]
    EOI_FLAG_TYPES = ["eoi_time", "eofi_time"]

    ANALYSIS_ITEMS_FILE = "analysis_template.json"
    STATE_FILE = "statefile.json"
    ERROR_LOG = "error.log"
