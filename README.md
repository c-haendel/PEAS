# PEAS - Parametric EIT Analysis Software

<p align="center">
  <img src="icon.png" alt="PEAS Icon" width="304" height="304"/>
</p>

**P**arametric **E**IT **A**nalysis **S**oftware

PEAS is a PyQt5-based GUI application for analyzing electrical impedance tomography (EIT) recordings. It provides a modular, configurable framework for research protocols with automatic breath/maneuver detection, EIT reconstruction, and customizable analysis pipelines.

## Key Features

* Analyze electrical impedance tomography (EIT) recordings
* Modular structure, freely configurable for every research protocol
* Simple and responsive graphical user interface
* Define measurement intervals and edit settings in GUI
* Platform independent (Windows, Linux, macOS)
* Load raw voltage files (Dräger *.eit, Sentec *.eit)
* Reconstruction with PyEIT (GREIT algorithm)
* Load reconstructed EIT images (Dräger *.bin, Sentec *.zri)
* Concatenate multiple input files
* Automatic detection of breaths and breathing maneuvers
* User-defined outputs composed of numerous `base items` and `operations`
* Export results as CSV or NPY

## Installation

### Prerequisites

* Python 3.8 or higher
* pip package manager

### Dependencies

Install all dependencies using pip:

```bash
pip install .
```

Or install from the repository:

```bash
pip install -e .
```

The following packages are required:
* **pyqtgraph** - Plotting and parameter tree UI
* **PyQt5** - GUI framework
* **numpy** - Numerical computation
* **scipy** - Signal processing (butter, sosfiltfilt, find_peaks, CubicSpline)
* **pyeit** - EIT image reconstruction (GREIT algorithm)

### Optional: Virtual Environment

```bash
python -m venv peas_env
source peas_env/bin/activate  # On Windows: peas_env\Scripts\activate
pip install -e .
```

## Quick Start

### GUI Mode

Run the application with a graphical interface:

```bash
python main.py
```

Or open a file directly:

```bash
python main.py path/to/file.eit
```

### Headless Mode

For batch processing without GUI:

```bash
python main.py path/to/file.eit --run
```

This will automatically load, reconstruct (if needed), analyze, and export results, then exit.

### Custom Template

Use a custom analysis template:

```bash
python main.py path/to/file.eit --template path/to/custom_template.json
```

## Configuration

PEAS uses JSON template files to define analysis intervals and items. The following templates are provided:

* **analysis_template.json** - Standard template with regular/deep breathing and vital capacity maneuvers
* **analysis_template_peep.json** - PEEP trial template (Costa et al. approach)

When defining new templates, it is recommended to start from an existing one. You can also activate/edit items from within the software (experimental feature).

## User Interface

The graphical interface is divided into three areas:

* **Left: Menu Bar** - Allows input of parameters and selection of intervals. Divided into sections typically handled sequentially during analysis.
* **Top: Overview Plot Window** - Displays the entire signal and the position of the plot window. Can be moved via drag-and-drop.
* **Bottom: Plot Window** - Shows a zoomed-in view of the signal. Zoom in time dimension with mouse wheel; Y-dimension is always auto-scaled. Interval input via click-and-drag.

The plots display a summation signal: voltage sum in gray (raw data) or impedance sum in white (reconstructed images).

## Operation

### Performing Analyses

1. **Loading Data**
   * **Voltage data**: Set reconstruction options and reference time point. If no reference is defined (value = 0), the time of the median voltage sum is chosen. Images are reconstructed automatically.
   * **Reconstructed images**: Load pre-reconstructed EIT images directly.
   
   The sampling rate is automatically determined.

2. **Determining Intervals**
   * Select an interval by clicking on the corresponding item in the menu bar or associated parameter.
   * Select the detector and adjust automatically displayed parameters.
   * Define interval start and end via click-and-drag on the plot window while an interval is selected.
   * After selecting an interval or changing the detector, detection is automatically executed and results are displayed.

3. **Saving, Calculation, and Export**
   * Output path is automatically chosen based on the filename but can be set manually.
   * Click "Export Results" to save intervals, settings, and calculate all output values.

   **Export format:**
   * 0-dimensional quantities → `scalars.csv`
   * 1- or 2-dimensional quantities → `{item_name}.csv`
   * 3+ dimensional quantities → `{item_name}.npy`

### Creating/Modifying Analysis Items

Analysis items are defined in the JSON template file. Each item can be:

* **Base item**: Direct computation from raw data with optional preprocessing
* **Custom item**: Composed of prerequisites (other items) and operations

Refer to [TECHNICAL.md](TECHNICAL.md) for complete documentation of all available base items, operations, and preprocessors.

### Example: Cumulative Functional ROI

This example creates a cumulative functional region of interest containing 90% of tidal variation:

```json
{
  "name": "rb_fROI_cumulative",
  "title": "cumulative functional ROI (90% tidal variation)",
  "prerequisites": [
    {
      "base_item": "tidal_image",
      "interval": "rb"
    }
  ],
  "operations": [
    {"name": "mean_over_time"},
    {"name": "maximum", "parameters": {"x2": 0}},
    {"name": "normalize_sum"},
    {"name": "cumulative_threshold", "parameters": {"threshold": 0.9}}
  ],
  "export": true
}
```

## File Formats

### Input Formats

* **.eit** - Raw voltage data (Dräger or Sentec format)
* **.bin** - Reconstructed EIT images (Dräger format)
* **.zri** - Reconstructed EIT images (Sentec format)
* **.npz** - Numpy arrays (custom format)

### Output Formats

* **.csv** - Comma-separated values for 0D, 1D, and 2D data
* **.npy** - Numpy binary format for 3D+ data

## Technical Reference

For detailed technical documentation including all detectors, base items, operations, and preprocessors, see [TECHNICAL.md](TECHNICAL.md).

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

## Author

Copyright © 2025 Claas Händel
