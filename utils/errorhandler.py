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

import traceback
from PyQt5.QtWidgets import QMessageBox, QDialog, QVBoxLayout, QLabel, QTextEdit, QPushButton, QScrollArea
from utils.globalsettings import GlobalSettings
import datetime
import sys


class CriticalError(Exception):
    """Exception raised for critical errors that require program termination."""
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class ErrorHandler:
    def __init__(self):
        self.log_file = GlobalSettings.ERROR_LOG
        sys.excepthook = self.excepthook

    def log_error(self, exception):
        # Log the error to the console
        print(f"ERROR: {exception}")
        # Write the error to a log file
        with open(self.log_file, "a") as file:
            file.write(f"{datetime.datetime.now()}: {exception}\n")

    def show_error(self, exception_type, exception_message, exception_traceback, verbose=True):
        # Create a dialog to show the error
        dialog = QDialog()
        dialog.setWindowTitle("Error")
        
        layout = QVBoxLayout()

        # Add the error type and message
        error_label = QLabel(f"<b>Exception Type:</b> {exception_type}<br><b>Exception Message:</b> {exception_message}")
        layout.addWidget(error_label)
        
        # Add a scrollable text area for the traceback
        if verbose:
            traceback_label = QLabel("<b>Traceback:</b>")
            layout.addWidget(traceback_label)
            
            traceback_text = QTextEdit()
            traceback_text.setReadOnly(True)
            traceback_text.setPlainText(exception_traceback)
            
            scroll_area = QScrollArea()
            scroll_area.setWidget(traceback_text)
            scroll_area.setWidgetResizable(True)
            scroll_area.setMinimumHeight(200)
            scroll_area.setMinimumWidth(400)
            layout.addWidget(scroll_area)
        
        # Add an OK button to close the dialog
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(dialog.accept)
        layout.addWidget(ok_button)
        
        dialog.setLayout(layout)
        dialog.exec_()

    def handle_exception(self, exception):
        # Get the exception type, message, and traceback
        exception_type = type(exception).__name__
        exception_message = str(exception)
        # Format the traceback to include file name and line number
        formatted_traceback = ''.join(traceback.format_exception(type(exception), exception, exception.__traceback__))

        # Log the exception
        self.log_error(", ".join((exception_type, exception_message, formatted_traceback)))

        # Display the detailed error message to the user
        self.show_error(exception_type, exception_message, formatted_traceback)

        # If it's a critical error, terminate the program
        if isinstance(exception, CriticalError):
            QMessageBox.critical(None, "Critical Error", "A critical error occurred. The application will now terminate.")
            sys.exit(1)

    def excepthook(self, exc_type, exc_value, exc_traceback):
        # Handle uncaught exceptions
        exception = exc_value
        self.handle_exception(exception)
