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
import argparse
from pathlib import Path

from PyQt5 import QtWidgets, QtCore

from utils.eitdatahandler import EITDataHandler
from utils.settingshandler import SettingsHandler
from utils.analysisitem import AnalysisItemManager
from utils.errorhandler import ErrorHandler
from utils.gui import GUI
from utils.calculationthread import CalculationThread

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', nargs='?', help='Input file path')
    parser.add_argument('--run', action='store_true', help='Run computation headlessly and exit')
    parser.add_argument('--template', type=str, help='Path to analysis template JSON file (default: analysis_template.json)')
    args = parser.parse_args()

    QtWidgets.QApplication.setHighDpiScaleFactorRoundingPolicy(QtCore.Qt.HighDpiScaleFactorRoundingPolicy.PassThrough)
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    app = QtWidgets.QApplication(sys.argv)

    error_handler = ErrorHandler()
    data_handler = EITDataHandler(error_handler)
    settings_handler = SettingsHandler(error_handler, template_path=args.template)
    analysis_item_manager = AnalysisItemManager(data_handler, settings_handler, error_handler)

    if args.run:
        if not args.file:
            print("Error: --run requires a file argument")
            sys.exit(1)

        calculation_thread = CalculationThread(data_handler, analysis_item_manager, error_handler)
        calculation_thread.error_occurred.connect(lambda e: (print(f"Error: {e}"), app.quit()))
        calculation_thread.start()

        def run_headless():
            file_path = Path(args.file)
            is_raw = data_handler.is_raw_file(file_path)

            settings_handler.set_value("output_path", file_path.with_suffix(""))
            settings_handler.set_append_state(2 if is_raw else 1)
            settings_handler.read_state_file()

            def on_load_done(_):
                calculation_thread.calculation_complete.disconnect(on_load_done)
                if is_raw:
                    calculation_thread.calculation_complete.connect(on_reconstruct_done)
                    calculation_thread.enqueue_task("reconstruct", settings_handler.get_param_dict(
                        settings_handler.param_recursive("reconstruction_algorithm")))
                else:
                    run_export()

            def on_reconstruct_done(_):
                calculation_thread.calculation_complete.disconnect(on_reconstruct_done)
                run_export()

            def run_export():
                settings_handler.write_state_file()
                settings_handler.export_analysis_json()
                calculation_thread.calculation_complete.connect(on_analyze_done)
                calculation_thread.enqueue_task("analyze", None)

            def on_analyze_done(_):
                calculation_thread.calculation_complete.disconnect(on_analyze_done)
                calculation_thread.stop()
                app.quit()

            calculation_thread.calculation_complete.connect(on_load_done)
            calculation_thread.enqueue_task("load_raw" if is_raw else "load_reconstructed", file_path)

        QtCore.QTimer.singleShot(0, run_headless)
        sys.exit(app.exec_())

    win = GUI(data_handler, settings_handler, analysis_item_manager, error_handler, args.file)
    win.show()
    sys.exit(app.exec_())
