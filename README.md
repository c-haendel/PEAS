# PEAS Documentation

# General/core features

**P**arametric **E**IT **A**nalysis **S**oftware

* analyse electrical impedance tomography (EIT) recordings
* modular structure, freely configurable for every research protocol
* simple and responsive graphical user inerface
* define measurement intervals and edit settings in GUI
* platform independent
* load raw voltage files (Dräger \*.eit)
* reconstruction with PyEIT (GREIT-algorithm)
* load reconstructed EIT images (Dräger \*.bin, Sentec \*.zri)
* concatenate multiple input files
* different options for automatic detection of breaths/breathing maneovres
* user defined outputs composed of numerous `base items` and `operations`
* export results as csv or npy


# Usage

Run in a directory with a file called analysis_template.json provided here.
An alternative analysis_template_peep.json is also provided as an example for a PEEP trial like Costa et al.
When defining new templates, it is suggested to start with an already existing one. It is also possible to activate/edit items from within the software, however this is experimental.


# User Interface

The graphical interface is divided into three areas:
* Left: Menu Bar allows for input of parameters and selection of intervals. The menu bar is divided into different sections, which are typically handled sequentially during an analysis.
* Top: Overview Plot Window always displays the entire signal and the position of the plot window. It can be moved via drag-and-drop.
* Bottom: Plot Window shows a zoomed-in view of the signal. Zooming (in the time dimension) is possible with the mouse wheel; in the y-dimension, the plot is always auto-scaled. Interval input is done via click-and-drag.

The plots each display a summation signal. If these are voltage values, the sum of the voltages is shown in gray. If the images are reconstructed, the impedance sum is shown in white.


# Operation

## Performing Analyses

1. Loading Data.

Data can be loaded in two ways:
* Loading voltage data. Setting reconstruction options and a reference time point. (If no reference time point is defined (value is 0), the time of the median voltage sum is chosen.) Reconstruction of images.
* Loading already reconstructed images.
The sampling rate is automatically determined.

2. Determining Intervals.

* The selection of an interval is made by clicking on the corresponding item in the menu bar or any associated parameter.
* Selection of the detector and the (automatically displayed) parameters of the detector.
* The definition of the start and end of an interval is done via click-and-drag on the plot window while an interval is selected.
* After selecting an interval or changing the detector, the detector is automatically executed, and the corresponding results are displayed in the plot window.

3. Saving, Calculation, and Output.

The output path is automatically chosen based on the opened file name, but can also be set manually.
A mouse click on export results saves the defined intervals and settings and calculates all defined output values in the output directory.

* All 0-dimensional quantities are saved in scalars.csv.
* Any 1- or 2-dimensional quantities are saved in a csv file named after the analysis items.
* Any 3- or more-dimensional quantities are saved in an npy file named after the analysis items.

## Creating/Modifying analysis items

This is initially done via the analysis_items.json. Please refer to the existing entries for structure.


# Technical Details

## detectors

Detectors can be of various types:

* maneuver: Detection of forced inspiration and expiration.
* breath: Detection of cyclic breathing patterns.

Detectors may also have parameters and output quantities, which are displayed in the graphical interface.
The result of a detection is stored internally as flags, which are displayed as points or lines in the plot window:

* eoe_time, eofe_time: End of an expiration or forced expiration is shown in purple.
* eoi_time, eofi_time: End of an inspiration or forced inspiration is shown in blue.
* Other lines are shown in gray.
* Points are marked with a circle.

### low_pass

Detector for cyclic breathing patterns. Suitable for resting and deep breathing.
A Butterworth low-pass filter is applied forward and backward to the summation signal. Then, local maxima and minima are determined. Starting from the end of the interval, minima/maxima within timeout after a minimum/maximum are removed.

Parameters:
* filter_order: Filter order (Default: 4)
* cutoff: Cutoff frequency of the filter (Default: 1.0 Hz)
* timeout: (Default: 2 s)

Placed flags:
* eoe_time: End of expiration. The first flag is always an eoe_time flag.
* eoi_time: End of inspiration. The last flag is always an eoi_time flag.

### min_max_min

Detector for forced inspiration and expiration. Particularly suitable for slow forced maneuvers.
It detects the global maximum in the interval as well as the global minima before and after it.

Placed flags:
* eoe_time: Time of end of expiration (before maximum)
* eofi_time: Time of end of forced inspiration
* eofi_z: Impedance value at the end of forced inspiration
* eofe_time: Time of end of forced expiration (after maximum)
* eofe_z: Impedance value at the end of forced expiration

### simple_peak_flow

Detector for forced inspiration and expiration. Particularly suitable for fast forced maneuvers. The implementation is based on the definition of the beginning of expiration in spirometry.

Parameters:
* sliding_window_size: Time interval in seconds between which the flow (slope of the impedance curve) is calculated. (Default: 0.01 s)

Output:
* pre_expired_fraction: The fraction of the volume already expired at the time of the start of forced expiration, measured as a fraction of the forced expired volume.

The sampling rate of the summation signal is first converted via cubic spline interpolation (new sampling rate is 2/sliding_window_size). The global maximum in the interval is detected, as well as the global minima before and after it. A line is drawn through the minimum slope after this maximum (peak expiratory flow). The intersection of this line with a horizontal line through the maximum gives the start time of forced expiration.

Placed flags:
* eoe_time: Time of end of expiration (before maximum)
* eofi_time: Time of the beginning of forced expiration as described above
* eofi_z: Impedance value at the end of forced inspiration
* eofe_time: Time of the end of forced expiration (after maximum)
* eofe_z: Impedance value at the end of forced expiration
* peak_expiratory_flow: Time and impedance value of the peak expiratory flow
* peak_expiratory_flow_tangent: Line through the peak expiratory flow

## Structure of analysis items

There are two types of analysis items:

* Automatically generated base items
* (User-defined) custom analysis items.

### base items

Parameters:
* base_item: One of the names of available base items listed below.
* parameters: Parameters passed to the base item.
* preprocessors: Optional list of preprocessors.
* interval: Name of the interval to which the base item should be applied.

Output:
* Numpy array with dimensions (samples/breaths, vd, rl[, i])

Where the indices mean:
* samples/breaths: This index describes the result in the time dimension, usually as the number of samples or breaths.
* vd: This index describes the result in the spatial dimension in the ventral-dorsal direction.
* rl: This index describes the result in the spatial dimension in the right-left direction.
* i: This optional index is used to access multiple results of a base item (e.g., `expiratory_time_constant`).

### Custom Analysis Items
* name: Arbitrary unique name. Used as the filename during export, so it must not contain special characters except "_".
* prerequisites: List of analysis items whose results are needed for the calculation.
* operations: Optional list of operations.

## preprocessing

### resample_over_image vd, rl

Discrete change in spatial resolution, particularly reduction in resolution.
The old and new (scaled) data grids are overlaid. Each new data point is the weighted mean of the corresponding old data points, weighted by the area overlap.

Inputs:
* Numpy array with dimensions (samples, vd0, rl0)

Parameters:
* `vd1`, `rl1`: New dimensions in the ventro-dorsal and right-left directions.

Outputs:
* Numpy array with dimensions (samples, vd1, rl1)

### low_pass_filter

Apply butterworth bidirectional low pass filter.

Parameters:
* `cutoff`: filter cutoff
* `filter_order`: filter order

## base items

### tidal_image (analog: expiratory_tidal_image)

Calculates tidal images for all detected breaths by subtracting the end-expiration from the previous end-inspiration.

Parameters:
* None

Output:
* numpy array with dimensions (breaths, vd, rl)

### expired_volume (analog: inspired_volume)

Difference in impedance between the time t after end-inspiration and end-expiration.

Parameters:

* `t`: Time in seconds after end-inspiration (Default: 1).

Output:

* Numpy array with dimensions (breaths, vd, rl)

### time_to_expire (analog: time to inspire)

Time in seconds between end-inspiration and expiration of fraction (part of the expiratory tidal difference).
Beforehand, a resampling conversion (cubic spline interpolation) to the sampling rate sf takes place.

Parameters:

* `fraction`: Value between 0 and 1.
* `sf`: Sampling rate in Hertz (Default: 100).

Output:

* Numpy array with dimensions (breaths, vd, rl)

### flow_when_expired (analog: flow_when_inspired)

Flow in AU/s at the time of expiration of fraction (part of the expiratory tidal difference).
Beforehand, a resampling conversion (cubic spline interpolation) to the sampling rate sf takes place.
Expiration has negative flow.

Parameters:
* `fraction`: Value between 0 and 1.
* `sf`: Sampling rate in Hertz (Default: 100).

Output:
* Numpy array with dimensions (breaths, vd, rl)

### peak_expiratory_flow

Maximum flow in AU/s between end-inspiration and end-expiration.
Beforehand, a resampling conversion (cubic spline interpolation) to the sampling rate sf takes place.
Expiration has negative flow.

Parameters:
* `sf`: Sampling rate in Hertz (Default: 100).

Output:
* Numpy array with dimensions (breaths, vd, rl)

### peak_expiratory_flow_time

Time of peak flow between end-inspiration and end-expiration in seconds since the end-inspiration.
Beforehand, a resampling conversion (cubic spline interpolation) to the sampling rate sf takes place.

Parameters:
* `sf`: Sampling rate in Hertz (Default: 100).

Output:
* Numpy array with dimensions (breaths, vd, rl)

### mean_expiratory_flow

Flow in AU/s between two fractions of expiration supplied with list.
Beforehand, a resampling conversion (cubic spline interpolation) to the sampling rate sf takes place.

Parameters:
* `fraction`: List of two values between 0 and 1.
* `sf`: Sampling rate in Hertz (Default: 100).

Output:
* Numpy array with dimensions (breaths, vd, rl)

### end_expiratory_lung_impedance (analog: end_inspiratory_lung_impedance)

Impedance at the time of end-expiration.

Parameters:
* None

Output:
* Numpy array with dimensions (breaths, vd, rl)

### expiratory_concavity

Measure of regional concavity of expiratory (positive) part of flow-volume-loop: Area enclosed by all data points in during expiration divided by convex hull of the same points. Results in value between 0 (most concave) and 1.
Beforehand, a resampling conversion (cubic spline interpolation) to the sampling rate sf takes place.

Parameters:
* `sf`: Sampling rate in Hertz (Default: 100).

Output:
* Numpy array with dimensions (breaths, vd, rl)

### expiratory_time_constant

Pixel-wise expiratory time constant. Data outside valid_range are discarded. An exponential regression is performed with the function 'z0 * np.exp(-t/tau) + c'. All variables are limited to positive values.

Parameters:
* `valid_range`: Data range used for regression (part of the expiratory tidal difference data range) (Default: [0.25, 0.75]).

Output:
* Numpy array with dimensions (breaths, vd, rl, i)

i: Index of the fit variable, i=0: z0; i=1: tau; i=2: c; i=3: R².

### custom_fits

Pixel-wise regression of any function. Data outside valid_range is discarded.

Parameters:
* `valid_range`: Data range used for regression (part of the expiratory tidal difference data range) (Default: [0.25, 0.75]).
* `func_str`: String representation of the regression function, e.g., 'z0 * np.exp(-t/tau) + c'
* `variables`: Names of the variables in func_str, starting with the dependent variable, e.g., ['t', 'z0', 'tau', 'c']
* `bounds` etc: See scipy.optimize.curve_fit, e.g., [0, "inf"]

Output:
* Numpy array with dimensions (breaths, vd, rl, i)

i: Index of the fit variable (according to the order of variables without the dependent variable, starting from 0), followed by the R² value of the fit.

Note: Symmetrical quantities (like in bi-exponential regression) are not sorted!

### breath_times

Length of the times of the globally detected breaths in seconds.
Since the detected breaths always end with an inspiration, there is no expiration time or breath cycle time for the last breath.

Parameters:
* None

Output:
* Numpy array with dimensions (breaths, 1, 1, i)

i: i=0: Inspiration time, i=1: Expiration time, i=2: Total (breath cycle length)

### passthrough

Data (after preprocessing) is output 1:1.

Parameters:
* None

Output:
* Numpy array with dimensions (samples, vd, rl)

### interval_data

Access interval attributes.

### detector_data

Access detector attributes.

Parameters:
* attr_name: requested attribute the detector exposes

## operations

Alle Operationen ignorieren nicht-numerische Daten.

### mean_over_time (analog: median_over_time)

Calculate mean over time dimension.

Input:
* numpy array with dimensions (t, vd, rl[, i])

Output:
* numpy array with dimensions (1, vd, rl[, i])

### mean_over_image (analog: median_over_image, sum_over_image, percentile_over_image)

Calculate mean over spacial dimensions.

Input:
* numpy array with dimensions (t, vd, rl[, i])
* (percentile)

Output:
* numpy array with Dimensions (t, 1, 1[, i])

### threshold

Threshold function. If a value is above the threshold it will be turned to 1, else to 0.

Input:
* numpy array with dimensions (t, vd, rl)

Parameters:
* `threshold`: absolute threshold.

Output:
* numpy array with dimensions (t, vd, rl)

### percentile

Relative threshold function. If a value is above the percentile of data, it will be turned to 1, else to 0.

Input:
* numpy array with dimensions (t, vd, rl)

Parameters:
* `percentile`: percentile, 0-100.

Output:
* numpy array with dimensions (t, vd, rl)

### cumulative_threshold

Cumulative threshold function. The minimal number of pixels summing to threshold will be set to 1, others to 0.

Input:
* numpy array with dimensions (1, vd, rl)

Parameters:
* `threshold`: absolute threshold.

Output:
* numpy array with dimensions (1, vd, rl)

### sum (analog: std)

Sum.

Input:
* numpy array

Parameters:
* `axis`, ...: see numpy.nansum.

Output:
* numpy array with fewer dimensions

### coefficient_of_variation (analog: global_inhomogeneity_index)

Standard deviation divided by mean.

Input:
* numpy array

Output:
* numpy array with dimensions (0)

### normalize_sum

Normalize to sum 1.

Input:
* numpy array with dimensions (1, vd, rl)

Output:
* numpy array with dimensions (1, vd, rl)

### normalize_max

Normalize to maximum 1.

Input:
* numpy array with dimensions (1, vd, rl)

Output:
* numpy array with dimensions (1, vd, rl)

### centroid

Calculate centroid of distribution along spacial dimensions.

Input:
* numpy array with dimensions (1, vd, rl)

Output:
* numpy array with dimensions (1, 1, 1, i=2)

i=0: centroid in ventro-dorsal direction; i=1: centroid in right-left direction.

### multiply, divide, subtract, add

Multiplication of `analysis items`.

Input:
* 2 numpy arrays with dimensions (t, vd, rl)

Parameters:
* `targets`: List of targets of the operation, e.g. [0, 1] for first and second element of the defined `prerequisites`. (Default 0)
* see numpy.multiply, numpy.divide, numpy.subtract, numpy.add

Output:
* numpy array with dimensions (t, vd, rl)

### slice_last

Selection of a result/slicing the last dimension for when `base item` or `operation` have multiple results (e.g. `expiratory_time_constant`; `centroid`)

Input:
* numpy array with dimensions (t, vd, rl, i)

Parameter:
* `index`: index of the result, starting with 0.

Output:
* numpy array with dimensions (t, vd, rl)

### slice_first

Slice into first (=time) dimension of data.

Input:
* numpy array with dimensions (t, vd, rl)

Parameter:
* `index`: index of desired slice

Output:
* numpy array with dimensions (1, vd, rl) or similar

### min, max
Minimum/maximum of the data.

Input:
* numpy array

Parameter:
* `axis`, ...: see numpy.nanmin, numpy.nanmax.

Output:
* numpy array with dimensions (0), ...

### minimum (analog: maximum)

Clip the data, replace all values over `x2` mit `x2`.

Input:
* numpy array

Parameter:
* `x2`: reference value.

Output:
* numpy array

### apply_mask

Apply a boolean mask, turning masked values to numpy.nan.

Input:
* numpy array
* boolean numpy array (mask)

Output:
* numpy array

### costa_approach

Operation to create a table of loss of compliance compared to maximum compliance in both directions of PEEP like described by Costa et al.

Input:
* numpy array (stack of tidal images) with dimensions (n, vd, rl)

Output:
* numpy array of dimensions (n, 2) (collapse, overdistension)


