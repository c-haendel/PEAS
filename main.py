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

import sys
from PyQt5 import QtWidgets, QtCore

from utils.eitdatahandler import EITDataHandler
from utils.settingshandler import SettingsHandler
from utils.analysisitem import AnalysisItemManager
from utils.errorhandler import ErrorHandler
from utils.gui import GUI

if __name__ == '__main__':
    QtWidgets.QApplication.setHighDpiScaleFactorRoundingPolicy(QtCore.Qt.HighDpiScaleFactorRoundingPolicy.PassThrough)
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    app = QtWidgets.QApplication(sys.argv)

    error_handler = ErrorHandler()
    data_handler = EITDataHandler(error_handler)
    settings_handler = SettingsHandler(error_handler)
    # analysis item manager needs data handler to directly access data for calculations
    # analysis item manager needs settings handler to access output path and interval list
    analysis_item_manager = AnalysisItemManager(data_handler, settings_handler, error_handler)
    file_path_arg = sys.argv[1] if len(sys.argv) > 1 else None
    win = GUI(data_handler, settings_handler, analysis_item_manager, error_handler, file_path_arg)
    win.show()
    sys.exit(app.exec_())
