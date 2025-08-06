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

from PyQt5 import QtCore
import queue

class CalculationThread(QtCore.QThread):
    # tasks return dict of menu items that must be updated
    calculation_complete = QtCore.pyqtSignal(dict)
    error_occurred = QtCore.pyqtSignal(object)
    status = QtCore.pyqtSignal(str)

    def __init__(self, data_handler, analysis_item_manager, error_handler):
        super().__init__()
        self.task_queue = queue.Queue()
        self.running = True
        self.data_handler = data_handler
        self.analysis_item_manager = analysis_item_manager
        self.error_handler = error_handler

    def run(self):
        while self.running:
            try:
                task_id, task_params = self.task_queue.get(timeout=1)
                if task_id == 'load_reconstructed':
                    self.status.emit("loading file...")
                    result = self.task_load_reconstructed(task_params, append=False)
                elif task_id == 'append_reconstructed':
                    self.status.emit("loading file...")
                    result = self.task_load_reconstructed(task_params, append=True)
                elif task_id == 'load_raw':
                    self.status.emit("loading file...")
                    result = self.task_load_raw(task_params, append=False)
                elif task_id == 'append_raw':
                    self.status.emit("loading file...")
                    result = self.task_load_raw(task_params, append=True)
                elif task_id == 'reconstruct':
                    self.status.emit("reconstructing file...")
                    result = self.task_reconstruct(task_params)
                elif task_id == 'analyze':
                    self.status.emit("calculating results...")
                    result = self.task_analyze(task_params)
                else:
                    raise ValueError(f'No task {task_id} defined')

                self.status.emit("done.")
                self.calculation_complete.emit(result)
            except queue.Empty:
                pass
            except Exception as e:
                print("exception occured in calculation thread")
                # Handle exceptions and signal the main thread
                self.error_occurred.emit(e)
                self.stop()

    def task_load_reconstructed(self, params, append):
        if not append: self.data_handler.wipe_data()
        sf = self.data_handler.load_reconstructed_input_file(params, append)
        return {'source_frequency': sf}

    def task_load_raw(self, params, append):
        if not append: self.data_handler.wipe_data()
        sf = self.data_handler.load_raw_input_file(params, append)
        return {'source_frequency': sf}

    def task_reconstruct(self, params):
        self.data_handler.reconstruct(**params)
        return {}

    def task_analyze(self, params):
        self.analysis_item_manager.calculate_for_export(self.status)
        return {}

    def enqueue_task(self, task_id, task_params):
        self.task_queue.put((task_id, task_params))

    def stop(self):
        self.running = False
