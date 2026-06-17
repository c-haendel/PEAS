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

import typing
from pathlib import Path
import numpy as np
from PyQt5 import QtWidgets, QtCore
import pyqtgraph as pg
import pyqtgraph.exporters as exporters
from pyqtgraph.parametertree import ParameterTree

from utils.calculationthread import CalculationThread
from utils.detector import Detector
from utils.analysisutils import FunctionType
from utils.globalsettings import GlobalSettings


class CustomPlotWidget(pg.PlotWidget):
    regionSelected = QtCore.pyqtSignal(tuple)  # custom signal

    def __init__(self, *args, **kwargs):
        super(CustomPlotWidget, self).__init__(*args, **kwargs)
        self.linear_region = None
        self.start_pos = None
        self.end_pos = None

    def setLinearRegion(self, linear_region):
        self.linear_region = linear_region

    def mousePressEvent(self, event):
        if self.linear_region is None:
            return
        mousePoint = self.plotItem.vb.mapSceneToView(event.pos())
        self.start_pos = mousePoint.x()
        self.end_pos = self.start_pos
        self.linear_region.setRegion([self.start_pos, self.end_pos])
        self.linear_region.show()

    def mouseMoveEvent(self, event):
        self.scene().sigMouseMoved.emit(event.pos())  # Manually emit the signal
        if self.linear_region is None or self.start_pos is None:
            return
        mousePoint = self.plotItem.vb.mapSceneToView(event.pos())
        self.end_pos = mousePoint.x()
        self.linear_region.setRegion([self.start_pos, self.end_pos])

    def mouseReleaseEvent(self, event):
        if self.linear_region is None or self.start_pos is None:
            return
        mousePoint = self.plotItem.vb.mapSceneToView(event.pos())
        self.end_pos = mousePoint.x()
        if self.end_pos < self.start_pos:
            self.start_pos, self.end_pos = self.end_pos, self.start_pos
        self.linear_region.setRegion([self.start_pos, self.end_pos])
        self.regionSelected.emit((self.start_pos, self.end_pos))  # emit custom signal
        self.start_pos = None
        self.end_pos = None


class GUI(QtWidgets.QMainWindow):
    def __init__(
        self,
        data_handler,
        settings_handler,
        analysis_item_manager,
        error_handler,
        file_path_arg,
    ):
        super(GUI, self).__init__()
        self.data_handler = data_handler
        self.settings_handler = settings_handler
        self.analysis_item_manager = analysis_item_manager
        self.error_handler = error_handler
        self.calculation_thread = CalculationThread(
            self.data_handler, self.analysis_item_manager, self.error_handler
        )
        self.calculation_thread.calculation_complete.connect(
            self.calculation_complete_handler
        )
        self.calculation_thread.error_occurred.connect(self.handle_error)
        self.calculation_thread.status.connect(self.set_status_bar)
        self.calculation_thread.start()

        self.showMaximized()
        self.initUI()

        # handle file argument if present
        if file_path_arg is not None:
            # TODO: set as raw filename if raw
            self.settings_handler.set_value("reconstructed_filename", file_path_arg)

    @staticmethod
    def find_param_item(tree_widget, param_name, parent_item=None):
        # this is necessary to set focus to element by name (instead of label/title)
        if parent_item is None:
            count = tree_widget.topLevelItemCount()
            for i in range(count):
                item = tree_widget.topLevelItem(i)
                found_item = GUI.find_param_item(tree_widget, param_name, item)
                if found_item:
                    return found_item
        else:
            if hasattr(parent_item, "param") and parent_item.param.name() == param_name:
                return parent_item
            for i in range(parent_item.childCount()):
                child_item = parent_item.child(i)
                found_item = GUI.find_param_item(tree_widget, param_name, child_item)
                if found_item:
                    return found_item
        return None

    def initUI(self):
        self.setWindowTitle("PEAS EIT Analysis Software")

        self.current_interval = None
        self.voltage_mode = False

        # menu
        self.measurementMenu = ParameterTree()
        self.measurementMenu.setParameters(
            self.settings_handler.parameter, showTop=False
        )
        self.measurementMenu.setMinimumWidth(300)  # TODO: check this number

        # overview plot
        self.overviewPlotWidget = pg.PlotWidget(title="")
        self.overviewPlotWidget.setMaximumHeight(150)
        self.overviewPlotWidget.enableAutoRange()
        self.overviewPlotWidget.setMouseEnabled(x=False, y=False)
        self.overviewPlotWidget.showAxis("bottom", False)
        self.overviewPlotWidget.showAxis("left", False)
        self.overviewPlot = self.overviewPlotWidget.plot([])
        self.overviewPlot.setZValue(GlobalSettings.PLOT_Z_DATA)

        # detail plot
        self.detailPlotWidget = CustomPlotWidget(title="")
        self.detailPlotWidget.setAutoVisible(y=True)
        self.detailPlotWidget.enableAutoRange(axis="y")  # maybe remove this
        self.detailPlotWidget.setMouseEnabled(x=True, y=False)  # maybe remove this
        self.detailPlotWidget.setLabel("bottom", "t", units="s")
        self.detailPlotWidget.setLabel("left", "Z", units="AU")
        self.detailPlot = self.detailPlotWidget.plot([])
        self.detailPlot.setZValue(GlobalSettings.PLOT_Z_DATA)

        # region connecting above plots
        self.region = pg.LinearRegionItem(clipItem=self.overviewPlot)
        self.region.setZValue(GlobalSettings.PLOT_Z_REGION)
        self.overviewPlotWidget.addItem(self.region, ignoreBounds=True)

        # cursor
        self.overviewCursor = pg.InfiniteLine(movable=False)
        self.overviewPlotWidget.addItem(self.overviewCursor, ignoreBounds=True)
        self.detailCursor = pg.InfiniteLine(movable=False)
        self.detailPlotWidget.addItem(self.detailCursor, ignoreBounds=True)

        # intervalregions: display the active region
        self.detailIntervalRegion = pg.LinearRegionItem(
            pen=pg.mkPen(GlobalSettings.PLOT_COLOR_REGION_PEN),
            brush=pg.mkBrush(GlobalSettings.PLOT_COLOR_REGION_BRUSH),
            hoverBrush=pg.mkBrush(GlobalSettings.PLOT_COLOR_REGION_HOVERBRUSH),
            clipItem=self.detailPlot,
        )
        self.detailIntervalRegion.setZValue(GlobalSettings.PLOT_Z_REGION)
        self.detailPlotWidget.addItem(self.detailIntervalRegion)
        self.detailPlotWidget.setLinearRegion(self.detailIntervalRegion)
        self.detailIntervalRegion.hide()
        self.overviewIntervalRegion = pg.LinearRegionItem(
            movable=False,
            pen=pg.mkPen(GlobalSettings.PLOT_COLOR_REGION_PEN),
            brush=pg.mkBrush(GlobalSettings.PLOT_COLOR_REGION_BRUSH),
            hoverBrush=pg.mkBrush(GlobalSettings.PLOT_COLOR_REGION_HOVERBRUSH),
            clipItem=self.overviewPlot,
        )
        self.overviewIntervalRegion.setZValue(GlobalSettings.PLOT_Z_REGION)
        self.overviewPlotWidget.addItem(self.overviewIntervalRegion)
        self.overviewIntervalRegion.hide()

        # build nested box layout
        menuLayout = QtWidgets.QVBoxLayout()
        menuLayout.addWidget(self.measurementMenu)
        plotLayout = QtWidgets.QVBoxLayout()
        plotLayout.addWidget(self.overviewPlotWidget, 1)
        plotLayout.addWidget(self.detailPlotWidget, 1)
        outerLayout = QtWidgets.QHBoxLayout()
        outerLayout.addLayout(menuLayout)
        outerLayout.addLayout(plotLayout, 1)

        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)
        self.centralWidget().setLayout(outerLayout)

        self.set_status_bar("")

        ### wire up interactive elements
        ## plot interconnection
        # when moving region in overviewPlot, change detailPLot
        self.region.sigRegionChanged.connect(self.update_detail_plot)
        # when detailPlot changes, change region in overviewPlot
        self.detailPlotWidget.sigRangeChanged.connect(self.update_region)
        ## menu changes
        self.settings_handler.parameter.sigTreeStateChanged.connect(
            self.react_to_changed_parameter_value
        )
        ## menu focus
        self.measurementMenu.itemSelectionChanged.connect(
            self.react_to_changed_parameter_selection
        )
        ## cursor
        self.overviewPlot.scene().sigMouseMoved.connect(self.overview_cursor_moved)
        self.proxy = pg.SignalProxy(
            self.detailPlotWidget.scene().sigMouseMoved,
            rateLimit=60,
            slot=self.detail_cursor_moved,
        )
        ## intervalregions
        self.detailPlotWidget.regionSelected.connect(self.interval_region_set)
        # self.detailIntervalRegion.sigRegionChanged.connect(self.interval_region_set)

    def lock_gui(self):
        self.setDisabled(True)
        self.detailPlot.hide()
        self.overviewPlot.hide()

    def calculation_complete_handler(self, result):
        self.set_bounds_on_load()
        self.detailPlot.show()
        self.overviewPlot.show()
        self.setEnabled(True)
        for key, value in result.items():
            self.settings_handler.set_value(key, value)

    def handle_error(self, exception):
        self.error_handler.handle_exception(exception)

    def set_status_bar(self, text):
        self.statusBar().showMessage(text)

    # for cursor, plots and regions
    def update_detail_plot(self):
        minX, maxX = self.region.getRegion()
        self.detailPlotWidget.setXRange(minX, maxX, padding=0)

    def update_region(self, _, viewRange):
        rgn = viewRange[0]
        self.region.setRegion(rgn)

    def overview_cursor_moved(self, point):
        self.cursor_moved(point, self.overviewPlotWidget)

    def detail_cursor_moved(self, point):
        self.cursor_moved(
            point[0], self.detailPlotWidget
        )  # manually emitted signal is tuple - first element is Point

    def cursor_moved(self, point, plotWidget):
        if plotWidget.sceneBoundingRect().contains(point):
            position = plotWidget.getPlotItem().vb.mapSceneToView(point)
            self.overviewCursor.setPos(position.x())
            self.detailCursor.setPos(position.x())
            self.set_status_bar("{:.3f} - {:.3f}".format(position.x(), position.y()))

    def interval_region_set(self, region):
        if self.current_interval is None:
            self.detailIntervalRegion.hide()
            self.overviewIntervalRegion.hide()
        else:
            if region[0] != region[1]:
                self.current_interval.start_time = region[0]
                self.current_interval.end_time = region[1]
                try:
                    self.current_interval.run_detector(
                        self.data_handler.get_impedance_data(
                            self.current_interval.start_time, self.current_interval.end_time
                        )
                    )
                    # if selected interval has non-zero length but no flags detected: set interval None
                    if len(self.current_interval.flags)==0:
                        self.current_interval.start_time = None
                        self.current_interval.end_time = None
                except Exception as e:
                    self.error_handler.handle_exception(e)
            else:
                self.current_interval.start_time = None
                self.current_interval.end_time = None
                self.current_interval.flags = []
            self.settings_handler.sync_tree_to_interval()  # for detector outputs to be reflected in the GUI
        # manual trigger of react_to_changed_parameter_selection() to refresh other plot
        self.react_to_changed_parameter_selection()
        # TODO: trigger to save state file

    def update_plots(self):
        # get preprocessing info
        z_sum, t = self.data_handler.get_plot_data(voltage_mode=self.voltage_mode)
        self.overviewPlot.setData(t, z_sum)
        self.detailPlot.setData(t, z_sum)

    def set_voltage_mode(self, activate):
        self.voltage_mode = activate
        if activate:
            # voltage mode
            self.detailPlotWidget.setLabel("left", "U", units="AU")
            self.detailPlot.setPen(pg.mkPen(GlobalSettings.PLOT_COLOR_DATA_V))
            self.overviewPlot.setPen(pg.mkPen(GlobalSettings.PLOT_COLOR_DATA_V))
        else:
            # impedance mode
            self.detailPlotWidget.setLabel("left", "Z", units="AU")
            self.detailPlot.setPen(pg.mkPen(GlobalSettings.PLOT_COLOR_DATA_Z))
            self.overviewPlot.setPen(pg.mkPen(GlobalSettings.PLOT_COLOR_DATA_Z))

    def set_bounds_on_load(self):
        plotData = self.data_handler.timestamps
        self.overviewIntervalRegion.setBounds((0, plotData[-1]))
        # TODO error on next line
        # self.detailIntervalRegion.setBounds((0,plotData[-1]))
        self.update_plots()

    # for menu changes
    def setup_with_analysis_template(self):
        self.settings_handler.read_analysis_template()
        self.settings_handler.create_intervals_from_json_data()
        self.settings_handler.rebuild_intervals_parameter()
        self.analysis_item_manager.interval_lookup_dict = self.settings_handler.get_interval_lookup_dict()
        self.analysis_item_manager.analysis_items = {}
        self.analysis_item_manager.create_analysis_items_from_json_data()

    def react_to_changed_parameter_value(self, param, changes):
        for param, change, value in changes:
            # do nothing on item delete
            if change == "parent":
                return False
            path = self.settings_handler.parameter.childPath(param)
            if "reconstructed_filename" in path:
                if value is None: break  # catch null changes
                self.lock_gui()
                self.set_voltage_mode(False)
                self.calculation_thread.enqueue_task("load_reconstructed", Path(value))
                self.settings_handler.set_value(
                    "output_path", Path(value).with_suffix("")
                )
                self.settings_handler.set_value("raw_filename", None)
                self.settings_handler.set_append_state(1)
                self.settings_handler.read_state_file()
            elif "reconstructed_filename_append" in path:
                if value is None: break  # catch null changes
                self.lock_gui()
                self.set_voltage_mode(False)
                self.calculation_thread.enqueue_task("append_reconstructed", Path(value))
                self.settings_handler.read_state_file()
            elif "raw_filename" in path:
                if value is None: break  # catch null changes
                self.lock_gui()
                self.set_voltage_mode(True)
                self.calculation_thread.enqueue_task("load_raw", Path(value))
                self.settings_handler.set_value(
                    "output_path", Path(value).with_suffix("")
                )
                self.settings_handler.set_value("reconstructed_filename", None)
                self.settings_handler.set_append_state(2)
                self.settings_handler.read_state_file()
            elif "raw_filename_append" in path:
                if value is None: break  # catch null changes
                self.lock_gui()
                self.set_voltage_mode(True)
                self.calculation_thread.enqueue_task("append_raw", Path(value))
                self.settings_handler.read_state_file()
            elif "reconstruct" in path:
                self.lock_gui()
                self.set_voltage_mode(False)
                recon_params = self.settings_handler.get_param_dict(
                    self.settings_handler.param_recursive("reconstruction_algorithm")
                )
                recon_params["source_frequency"] = self.settings_handler.get_value("source_frequency")
                self.calculation_thread.enqueue_task("reconstruct", recon_params)
            elif "save_analysis_state" in path:
                self._save_analysis_state()
            elif "export_results" in path:
                # TODO: grey out export button in voltage mode
                if not self.voltage_mode:
                    self._save_analysis_state()
                    self.export_interval_plots()
                    self.lock_gui()
                    self.calculation_thread.enqueue_task("analyze", None)
            elif "analysis_template" in path:
                self.setup_with_analysis_template()
            elif "edit_analyses" in path:
                data_source = self.analysis_item_manager.analysis_items
                if len(data_source)>0:
                    self.setup_window = SetupWindow(data_source)
                    result = self.setup_window.exec()
                    if result == QtWidgets.QDialog.Accepted:
                        self.analysis_item_manager.save_analysis_items_to_json_data()
                        self.settings_handler.write_analysis_template()
                    self.setup_with_analysis_template() # (discard changes and) refresh GUI
            elif "edit_intervals" in path:
                data_source = self.settings_handler.interval_list
                if len(data_source)>0:
                    self.setup_window = SetupWindow(data_source)
                    result = self.setup_window.exec()
                    if result == QtWidgets.QDialog.Accepted:
                        self.settings_handler.save_intervals_to_json_data()
                        self.settings_handler.write_analysis_template()
                    self.setup_with_analysis_template() # (discard changes and) refresh GUI

            # current_interval is always set correctly according to focus
            # all settings concerning intervals must be set explicitly because data is in objects, not tree
            elif (
                self.current_interval is not None and self.current_interval.name in path
            ):
                # detector or detector settings changed
                if any(item.endswith("_detector") for item in path):
                    # detector has changed
                    if path[-1].endswith("_detector"):
                        self.current_interval.detector = Detector.initialize_detector(
                            value, self.settings_handler
                        )
                        self.settings_handler.rebuild_intervals_parameter()
                        # set focus using custom find function
                        self.measurementMenu.setCurrentItem(
                            GUI.find_param_item(
                                self.measurementMenu,
                                self.settings_handler.param_recursive(path[-1]).name(),
                            )
                        )
                    # detector setting has changed
                    else:
                        self.settings_handler.sync_interval_to_tree()
                    # anything has changed: run detector
                    self.current_interval.run_detector(
                        self.data_handler.get_impedance_data(
                            self.current_interval.start_time,
                            self.current_interval.end_time,
                        )
                    )
                    self.settings_handler.sync_tree_to_interval()  # for detector outputs to be reflected in the GUI
                    self.hide_all_flags()
                    self.plot_flags_for_interval(self.current_interval)
                # interval settings changed
                else:
                    # currently no settings available
                    pass

    def react_to_changed_parameter_selection(self):
        for selectedItem in self.measurementMenu.selectedItems():
            path = self.settings_handler.parameter.childPath(selectedItem.param)
            self.current_interval = None
            self.hide_all_flags()
            self.detailIntervalRegion.hide()
            self.overviewIntervalRegion.hide()
            # cycle through all intervals and highlight accordingly
            for interval in self.settings_handler.get_interval_lookup_dict().values():
                if interval.name in path:
                    self.current_interval = interval
                    if (
                        interval.start_time is not None
                        and interval.end_time is not None
                    ):
                        # set all elements according to current interval
                        self.detailIntervalRegion.setRegion(
                            (interval.start_time, interval.end_time)
                        )
                        self.overviewIntervalRegion.setRegion(
                            (interval.start_time, interval.end_time)
                        )
                        self.plot_flags_for_interval(self.current_interval)
                        # show all element
                        self.detailIntervalRegion.show()
                        self.overviewIntervalRegion.show()

    def plot_flag_function(self, flag):
        function = flag.function
        function_type = function.function_type
        plot_item = None

        if function_type == FunctionType.POINT:
            x = function.params["x"]
            y = function.params["y"]
            plot_item = self.detailPlotWidget.plot([x], [y], pen=None, symbol="o")

        elif function_type == FunctionType.LINE:
            if flag.flag_type == "eoi_time" or flag.flag_type == "eofi_time":
                pen = pg.mkPen(GlobalSettings.PLOT_COLOR_EOI, width=2)  # blue
            elif flag.flag_type == "eoe_time" or flag.flag_type == "eofe_time":
                pen = pg.mkPen(GlobalSettings.PLOT_COLOR_EOE, width=2)  # purple
            else:
                pen = pg.mkPen(
                    GlobalSettings.PLOT_COLOR_FLAG_DEFAULT
                )  # transparent gray

            if "slope" in function.params:
                slope = function.params["slope"]
                intercept = function.params["intercept"]
                # Calculate the angle in degrees
                angle = np.arctan(slope) * (180 / np.pi)
                # Create an InfiniteLine with the calculated angle and intercept
                plot_item = pg.InfiniteLine(angle=angle, pos=(0, intercept), pen=pen)
                self.detailPlotWidget.addItem(plot_item)
            else:
                x_val = function.params["x"]
                # Create a vertical InfiniteLine
                plot_item = pg.InfiniteLine(angle=90, pos=x_val, pen=pen)
                self.detailPlotWidget.addItem(plot_item)

        plot_item.setZValue(GlobalSettings.PLOT_Z_FLAG)
        flag.plot_item = plot_item

    def plot_flags_for_interval(self, selected_interval):
        for flag in selected_interval.flags:
            if flag.plot_item is None:
                self.plot_flag_function(flag)
            else:
                flag.plot_item.show()

    def hide_all_flags(self):
        for interval in self.settings_handler.interval_list:
            for flag in interval.flags:
                if flag.plot_item is not None:
                    flag.plot_item.hide()

    def _save_analysis_state(self):
        self.settings_handler.write_state_file()
        self.settings_handler.export_analysis_json()

    def export_interval_plots(self):
        def export_interval_plot(interval):
            output_path = Path(self.settings_handler.get_value('output_path')) / (interval.name + ".png")
            exporter = exporters.ImageExporter(self.detailPlotWidget.plotItem)
            exporter.parameters()['width'] = GlobalSettings.PLOT_SCREENSHOT_WIDTH 
            exporter.export(str(output_path))
        # save current view
        original_x_limits = self.detailPlotWidget.plotItem.vb.viewRange()[0]
        for interval in self.settings_handler.interval_list:
            if interval.screenshot and interval.start_time:
                self.detailPlotWidget.plotItem.vb.setXRange(interval.start_time, interval.end_time)
                self.detailIntervalRegion.setRegion(
                    (interval.start_time, interval.end_time)
                )
                self.detailIntervalRegion.show()
                self.plot_flags_for_interval(interval)

                export_interval_plot(interval)

                self.hide_all_flags()
                self.detailIntervalRegion.hide()
        # restore current view
        self.current_interval = None
        self.detailPlotWidget.plotItem.vb.setXRange(original_x_limits[0], original_x_limits[1])


class ContentType:
    """Enum-like class to define cell editing types"""
    TEXT = "text"
    CHECKBOX = "checkbox"
    DROPDOWN = "dropdown"
    DICT = "dict"
    DICTLIST = "dictlist"
    CUSTOM = "custom"


class SetupWindow(QtWidgets.QDialog):
    # Column definitions for intervals
    intervals_columns = {
        "name": {"attr": "name", "content_type": ContentType.TEXT, "table_view": True},
        "title": {"attr": "title", "content_type": ContentType.TEXT, "table_view": True},
        "save screenshot": {"attr": "screenshot", "content_type": ContentType.CHECKBOX, "table_view": True}
    }

    # Column definitions for analysis items
    analyses_columns = {
        "name": {"attr": "name", "content_type": ContentType.TEXT, "table_view": True},
        "title": {"attr": "title", "content_type": ContentType.TEXT, "table_view": True},
        "unit": {"attr": "unit", "content_type": ContentType.TEXT, "table_view": False},
        "identifier": {"attr": "identifier", "content_type": ContentType.TEXT, "table_view": False},
        "comment": {"attr": "comment", "content_type": ContentType.TEXT, "table_view": False},
        "base item": {"attr": "base_item", "content_type": ContentType.TEXT, "table_view": False},
        "interval": {"attr": "interval", "content_type": ContentType.TEXT, "table_view": False},
        "export": {"attr": "export", "content_type": ContentType.CHECKBOX, "table_view": True},
        'parameters': {"attr": "parameters", "content_type": ContentType.DICT, "table_view": False},
        'preprocessors': {"attr": "preprocessors", "content_type": ContentType.DICTLIST, "table_view": False},
        'prerequisites': {"attr": "prerequisites", "content_type": ContentType.DICTLIST, "table_view": False},
        'operations': {"attr": "operations", "content_type": ContentType.DICTLIST, "table_view": False}
    }

    def __init__(self, data_source):
        """
        Initialize SetupWindow
        
        :param data_source: Source of data (settings_handler or analysis_item_manager)
        """
        super().__init__()
        self.data_source = data_source

        # Determine columns
        if isinstance(self.data_source, dict):
            self.columns = self.analyses_columns
            self.title = "Analysis Items Setup"
        else:
            self.columns = self.intervals_columns
            self.title = "Intervals Setup"

        self.resize(800, 600)

        self.setWindowTitle(self.title)

        self.table = QtWidgets.QTableWidget()
        vbox = QtWidgets.QVBoxLayout()
        vbox.addWidget(self.table)
        button_hbox = QtWidgets.QHBoxLayout()
        self.cancel_button = QtWidgets.QPushButton("Cancel", self) # reload intervals/analyses from json
        self.save_button = QtWidgets.QPushButton("Save", self) # save intervals/analyses to json

        self.cancel_button.clicked.connect(self.reject)
        self.save_button.clicked.connect(self.accept)

        button_hbox.addWidget(self.cancel_button)
        button_hbox.addWidget(self.save_button)
        vbox.addLayout(button_hbox)

        self.setLayout(vbox)

        # Setup table with data
        self.setup_dynamic_table(
            self.table, 
            self.columns, 
            self.data_source
        )

    def setup_dynamic_table(self, table: QtWidgets.QTableWidget, 
                             columns: typing.Dict[str, typing.Dict[str, typing.Any]], 
                             data_source: typing.Union[list, dict]):
        """
        Dynamically setup table with given columns and data, supporting various edit types.
        Supports both list and dict as data_source.
        """

        table.verticalHeader().setVisible(False)

        # Normalize keys based on data_source type
        if isinstance(data_source, dict):
            row_keys = list(data_source.keys())
        else:
            row_keys = list(range(len(data_source)))

        # filter columns to be displayed
        columns = {
            key: value
            for key, value in columns.items()
            if value.get("table_view") is True
        }

        table.setColumnCount(len(columns))
        table.setHorizontalHeaderLabels(list(columns.keys()))
        table.setRowCount(len(row_keys) + 1)  # One extra for the Add button

        for row_index, key in enumerate(row_keys):
            item = data_source[key]
            for col_num, col_def in enumerate(columns.values()):
                if col_def.get('content_type', None) == ContentType.CHECKBOX:
                    checkbox = QtWidgets.QCheckBox()
                    checkbox.setChecked(bool(getattr(item, col_def['attr'], False)))
                    table.setCellWidget(row_index, col_num, checkbox)
                    checkbox.setDisabled(True)
                else:
                    table_item = QtWidgets.QTableWidgetItem(str(getattr(item, col_def['attr'], None) or ""))
                    if col_def.get('edit_mode') == 'multiline':
                        table_item.setFlags(table_item.flags() | QtCore.Qt.ItemIsEditable)
                    table.setItem(row_index, col_num, table_item)
                    table_item.setFlags(table_item.flags() & ~QtCore.Qt.ItemIsEditable)

        # Handle double-click to open dialog
        def handle_double_click(row, _):
            if row < len(row_keys):
                actual_key = row_keys[row]
                self.open_edit_dialog(actual_key, data_source)

        table.cellDoubleClicked.connect(handle_double_click)

        # Add button at bottom
        add_button = QtWidgets.QPushButton("Add")
        add_button.clicked.connect(lambda: self.open_edit_dialog(None, data_source))
        add_row_index = table.rowCount() - 1
        table.setSpan(add_row_index, 0, 1, table.columnCount())
        table.setCellWidget(add_row_index, 0, add_button)

        table.resizeColumnsToContents()

    def open_edit_dialog(self, key, data_source):
        """
        Opens the edit dialog for a specific row
        
        :param row: Row index of the item to edit
        :param data_source: List of objects to populate the edit dialog
        """
        if key is None:
            if isinstance(data_source, dict):
                item = type(data_source.values()[0])()
            else:
                item = type(data_source[0])()
        else:
            item = data_source[key]
        dialog = EditDialog(item, self.columns)
        dialog.exec_()

        self.table.clear()
        self.setup_dynamic_table(
            self.table,
            self.columns,
            self.data_source
        )


class EditDialog(QtWidgets.QDialog):
    def __init__(self, item, columns):
        super().__init__()
        self.item = item
        self.columns = columns
        self.setWindowTitle("Edit Item")
        self.inputs = {}

        layout = QtWidgets.QVBoxLayout()
        form_layout = QtWidgets.QGridLayout()

        row = 0
        for label, config in self.columns.items():
            attr = config["attr"]
            value = getattr(self.item, attr, None)

            form_layout.addWidget(QtWidgets.QLabel(label), row, 0)

            if config["content_type"] == ContentType.TEXT:
                widget = QtWidgets.QLineEdit(str(value) if value is not None else "")
            elif config["content_type"] == ContentType.CHECKBOX:
                widget = QtWidgets.QCheckBox()
                widget.setChecked(bool(value))
            elif config["content_type"] in (ContentType.DICT, ContentType.DICTLIST):
                widget = self.build_tree_widget(value)
                widget.setMinimumHeight(150)
            else:
                continue

            self.inputs[attr] = widget
            form_layout.addWidget(widget, row, 1)
            row += 1

        layout.addLayout(form_layout)

        self.save_button = QtWidgets.QPushButton("Save")
        self.save_button.clicked.connect(self.save_changes)
        layout.addWidget(self.save_button)

        self.setLayout(layout)

    def build_tree_widget(self, data):
        tree = QtWidgets.QTreeWidget()
        tree.setColumnCount(2)
        tree.setHeaderLabels(["Key", "Value"])
        tree.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        tree.customContextMenuRequested.connect(lambda pos: self.open_context_menu(tree, pos))

        def populate_tree(parent, data):
            if isinstance(data, dict):
                for k, v in data.items():
                    item = QtWidgets.QTreeWidgetItem([str(k), str(v) if not isinstance(v, (dict, list)) else ""])
                    item.setFlags(item.flags() | QtCore.Qt.ItemIsEditable)
                    parent.addChild(item)
                    if isinstance(v, (dict, list)):
                        populate_tree(item, v)
            elif isinstance(data, list):
                for i, entry in enumerate(data):
                    item = QtWidgets.QTreeWidgetItem([f"[{i}]", ""])
                    item.setFlags(item.flags() | QtCore.Qt.ItemIsEditable)
                    parent.addChild(item)
                    populate_tree(item, entry)

        if data:
            populate_tree(tree.invisibleRootItem(), data)

        tree.expandAll()
        return tree

    def open_context_menu(self, tree, pos):
        item = tree.itemAt(pos)
        menu = QtWidgets.QMenu()

        add_action = menu.addAction("Add Child")
        delete_action = menu.addAction("Delete")

        action = menu.exec_(tree.viewport().mapToGlobal(pos))

        if action == add_action:
            new_item = QtWidgets.QTreeWidgetItem(["new_key", ""])
            new_item.setFlags(new_item.flags() | QtCore.Qt.ItemIsEditable)
            if item:
                item.addChild(new_item)
                item.setExpanded(True)
            else:
                tree.invisibleRootItem().addChild(new_item)

        elif action == delete_action and item:
            parent = item.parent()
            if parent:
                parent.removeChild(item)
            else:
                tree.invisibleRootItem().removeChild(item)

    def save_changes(self):
        def coerce_value(value):
            if value == "":
                return None
            for type_ in (int, float, str):
                try:
                    return type_(value)
                except ValueError:
                    continue
            return value

        for label, config in self.columns.items():
            attr_name = config["attr"]
            widget = self.inputs.get(attr_name)

            if isinstance(widget, QtWidgets.QLineEdit):
                value = widget.text()
                setattr(self.item, attr_name, coerce_value(value))

            elif isinstance(widget, QtWidgets.QCheckBox):
                setattr(self.item, attr_name, widget.isChecked())

            elif isinstance(widget, QtWidgets.QTreeWidget):
                def parse_item(item):
                    # Check if this has children
                    if item.childCount() == 0:
                        return coerce_value(item.text(1))
                    elif all(item.child(i).text(0).startswith("[") for i in range(item.childCount())):
                        # Looks like a list
                        return [parse_item(item.child(i)) for i in range(item.childCount())]
                    else:
                        return {
                            item.child(i).text(0): parse_item(item.child(i))
                            for i in range(item.childCount())
                        }

                root = widget.invisibleRootItem()
                if config["content_type"] == ContentType.DICT:
                    parsed = {}
                    for i in range(root.childCount()):
                        child = root.child(i)
                        parsed[child.text(0)] = parse_item(child)
                    setattr(self.item, attr_name, parsed)

                elif config["content_type"] == ContentType.DICTLIST:
                    parsed = []
                    for i in range(root.childCount()):
                        child = root.child(i)
                        parsed.append(parse_item(child))
                    setattr(self.item, attr_name, parsed)

        self.accept()
