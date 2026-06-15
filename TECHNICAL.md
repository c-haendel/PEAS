# PEAS Technical Reference

This document provides detailed technical documentation for all detectors, preprocessors, base items, and operations available in PEAS.

## Table of Contents

1. [Data Dimension Convention](#data-dimension-convention)
2. [Detectors](#detectors)
3. [Preprocessors](#preprocessors)
4. [Base Items](#base-items)
5. [Operations](#operations)
6. [Analysis Item Structure](#analysis-item-structure)

---

## Data Dimension Convention

All data arrays follow a consistent dimension convention:

**(samples/breaths, vd, rl [, i])**

* **samples/breaths**: Time dimension (number of samples or breaths)
* **vd**: Spatial dimension in ventral-dorsal direction
* **rl**: Spatial dimension in right-left direction
* **i**: Optional index for multiple results (e.g., fit parameters)

---

## Detectors

Detectors identify breathing patterns and maneuvers in EIT recordings. They are classified into two types:

* **breath**: Detection of cyclic breathing patterns
* **maneuver**: Detection of forced inspiration and expiration

Detector results are stored as flags, displayed in the plot window:
* **eoe_time, eofe_time** (End of Expiration): Purple vertical lines
* **eoi_time, eofi_time** (End of Inspiration): Blue vertical lines
* **Other lines**: Gray
* **Points**: Marked with circles

### Available Detectors

#### low_pass (breath detector)

Detector for cyclic breathing patterns. Suitable for resting and deep breathing.

**Algorithm:**
1. Apply Butterworth low-pass filter (forward and backward) to the summation signal
2. Detect local maxima and minima
3. Remove minima/maxima within `timeout` after another minimum/maximum (starting from end of interval)

**Parameters:**
| Parameter | Type | Default | Unit | Description |
|-----------|------|---------|------|-------------|
| `filter_order` | int | 4 | - | Butterworth filter order |
| `cutoff` | float | 1.0 | Hz | Filter cutoff frequency |
| `timeout` | float | 2.0 | s | Minimum time between peaks |

**Flags placed:**
* `eoe_time`: End of expiration (first flag is always eoe_time)
* `eoi_time`: End of inspiration (last flag is always eoi_time)

---

#### min_max_min (maneuver detector)

Detector for forced inspiration and expiration. Particularly suitable for slow forced maneuvers.

**Algorithm:**
1. Find global maximum in the interval
2. Find global minimum before the maximum
3. Find global minimum after the maximum

**Parameters:** None

**Flags placed:**
* `eoe_time`: Time of end of expiration (before maximum)
* `eofi_time`: Time of end of forced inspiration (at maximum)
* `eofi_z`: Impedance value at end of forced inspiration (horizontal line)
* `eofe_time`: Time of end of forced expiration (after maximum)
* `eofe_z`: Impedance value at end of forced expiration (horizontal line)

---

#### simple_peak_flow (maneuver detector)

Detector for forced inspiration and expiration. Particularly suitable for fast forced maneuvers. Implementation is based on spirometry definition of expiration onset.

**Algorithm:**
1. Resample summation signal via cubic spline interpolation (new rate = 2/sliding_window_size)
2. Detect global maximum and minima before/after
3. Find point of steepest negative slope (peak expiratory flow)
4. Draw tangent through peak flow point
5. Intersection of tangent with horizontal line at maximum gives forced expiration start

**Parameters:**
| Parameter | Type | Default | Unit | Description |
|-----------|------|---------|------|-------------|
| `sliding_window_size` | float | 0.01 | s | Time interval for flow calculation |

**Output:**
* `back_extrapolated_fraction`: Fraction of volume already expired at forced expiration onset (readonly)

**Flags placed:**
* `eoe_time`: Time of end of expiration (before maximum)
* `eofi_time`: Time of forced expiration onset (tangent intersection)
* `eofi_z`: Impedance at end of forced inspiration (horizontal line)
* `eofe_time`: Time of end of forced expiration (after maximum)
* `eofe_z`: Impedance at end of forced expiration (horizontal line)
* `peak_expiratory_flow`: Point at peak expiratory flow
* `peak_expiratory_flow_tangent`: Tangent line through peak flow

---

### Disabled Detectors

The following detectors exist in the codebase but are currently disabled:

* **reverse_cooldown**: Breath detector with cooldown-based peak filtering
* **manual_breath**: Pseudo-detector for manual breath placement

---

## Preprocessors

Preprocessors are applied to data before base item calculation.

### resample_over_image

Discrete change in spatial resolution, particularly for resolution reduction.

**Algorithm:**
1. Overlay old and new (scaled) data grids
2. Each new data point is the weighted mean of corresponding old data points
3. Weights are based on area overlap

**Input:** Numpy array with dimensions (samples, vd₀, rl₀)

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `vd` | int | New ventral-dorsal dimension |
| `rl` | int | New right-left dimension |

**Output:** Numpy array with dimensions (samples, vd, rl)

---

### low_pass_filter

Apply Butterworth bidirectional low-pass filter over time dimension.

**Input:** Numpy array with dimensions (samples, vd, rl)

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `cutoff` | float | Filter cutoff frequency (Hz) |
| `filter_order` | int | Butterworth filter order |

**Output:** Numpy array with dimensions (samples, vd, rl)

---

### strip_all_nan

Remove rows and columns that contain only NaN values.

**Input:** Numpy array with dimensions (samples, vd, rl)

**Parameters:** None

**Output:** Numpy array with dimensions (samples, vd', rl') where vd' ≤ vd and rl' ≤ rl

---

## Base Items

Base items compute fundamental measurements from EIT data. They can be used directly or as prerequisites for custom analysis items.

### Volume and Impedance

#### tidal_image / expiratory_tidal_image

Calculates tidal images for all detected breaths.

**Algorithm:**
* `tidal_image`: End-inspiration impedance minus preceding end-expiration impedance
* `expiratory_tidal_image`: End-inspiration impedance minus following end-expiration impedance

**Parameters:** None

**Input:** Interval with breath flags (eoe_time, eoi_time)

**Output:** Numpy array with dimensions (breaths, vd, rl)

---

#### expired_volume / inspired_volume

Difference in impedance between time t after end-inspiration and end-expiration.

**Algorithm:**
* `expired_volume`: Z(end-expiration) - Z(end-inspiration + t)
* `inspired_volume`: Z(end-expiration + t) - Z(end-inspiration)

**Parameters:**
| Parameter | Type | Default | Unit | Description |
|-----------|------|---------|------|-------------|
| `t` | float | 1.0 | s | Time after end-inspiration |

**Input:** Interval with breath flags

**Output:** Numpy array with dimensions (breaths, vd, rl)

---

### Time-Based Measurements

#### time_to_expire / time_to_inspire

Time required to expire/inspire a fraction of the tidal volume.

**Algorithm:**
1. Resample breath data via cubic spline interpolation to sampling rate `sf`
2. Find time when fraction of tidal change is reached

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fraction` | float | - | Value between 0 and 1 |
| `sf` | float | 100 | Sampling rate in Hz |

**Input:** Interval with breath flags

**Output:** Numpy array with dimensions (breaths, vd, rl)

**Units:** seconds

---

#### breath_times

Length of inspiration, expiration, and total breath cycle times.

**Algorithm:** Calculate time differences between consecutive flags.

**Parameters:** None

**Input:** Interval with breath flags (eoe_time, eoi_time)

**Output:** Numpy array with dimensions (breaths, 1, 1, 3)

**Index i:**
* i=0: Inspiration time
* i=1: Expiration time
* i=2: Total breath cycle time

**Note:** Since detected breaths always end with inspiration, there is no expiration time for the last breath.

---

### Flow-Based Measurements

#### flow_when_expired / flow_when_inspired

Flow at the time when a fraction of tidal volume has been expired/inspired.

**Algorithm:**
1. Resample breath data via cubic spline interpolation to sampling rate `sf`
2. Calculate flow (time derivative) at the time when fraction is reached

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fraction` | float | - | Value between 0 and 1 |
| `sf` | float | 100 | Sampling rate in Hz |

**Input:** Interval with breath flags

**Output:** Numpy array with dimensions (breaths, vd, rl)

**Units:** AU/s (Arbitrary Units per second; negative for expiration)

---

#### peak_expiratory_flow

Maximum flow magnitude during expiration.

**Algorithm:**
1. Resample breath data via cubic spline interpolation to sampling rate `sf`
2. Calculate flow (time derivative)
3. Find maximum (absolute) flow

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sf` | float | 100 | Sampling rate in Hz |

**Input:** Interval with breath flags

**Output:** Numpy array with dimensions (breaths, vd, rl)

**Units:** AU/s (negative for expiration)

---

#### peak_expiratory_flow_time

Time of peak flow during expiration, measured from start of expiration.

**Algorithm:**
1. Resample breath data via cubic spline interpolation to sampling rate `sf`
2. Calculate flow (time derivative)
3. Find time of maximum (absolute) flow

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sf` | float | 100 | Sampling rate in Hz |

**Input:** Interval with breath flags

**Output:** Numpy array with dimensions (breaths, vd, rl)

**Units:** seconds (from end-inspiration)

---

#### mean_expiratory_flow

Mean flow between two fractions of expiration.

**Algorithm:**
1. Resample breath data via cubic spline interpolation to sampling rate `sf`
2. Find times for both fractions
3. Calculate mean flow between these times

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fraction` | list | - | Two values between 0 and 1, e.g., [0.2, 0.8] |
| `sf` | float | 100 | Sampling rate in Hz |

**Input:** Interval with breath flags

**Output:** Numpy array with dimensions (breaths, vd, rl)

**Units:** AU/s (negative for expiration)

---

### Absolute Impedance

#### end_expiratory_lung_impedance / end_inspiratory_lung_impedance

Lung impedance at end-expiration or end-inspiration.

**Parameters:** None

**Input:** Interval with breath flags

**Output:** Numpy array with dimensions (breaths, vd, rl)

---

### Advanced Analysis

#### expiratory_concavity

Measure of regional concavity of the expiratory flow-volume loop.

**Algorithm:**
1. Resample breath data via cubic spline interpolation to sampling rate `sf`
2. For each pixel, construct flow-volume curve
3. Calculate area enclosed by the curve
4. Calculate area of convex hull of the same points
5. Concavity = 1 - (area / convex_hull_area)

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sf` | float | 100 | Sampling rate in Hz |

**Input:** Interval with breath flags

**Output:** Numpy array with dimensions (breaths, vd, rl)

**Range:** 0 to 1 (0 = most concave, 1 = convex/linear)

---

#### expiratory_time_constant

Pixel-wise expiratory time constant using exponential regression.

**Algorithm:**
1. Identify data within valid_range (fraction of tidal change)
2. Fit exponential decay: z(t) = z₀ · exp(-t/τ) + c
3. All parameters constrained to positive values
4. Calculate R² for goodness of fit

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `valid_range` | list | [0.25, 0.75] | Data range as fraction of tidal change |

**Input:** Interval with breath flags

**Output:** Numpy array with dimensions (breaths, vd, rl, 4)

**Index i:**
* i=0: z₀ (initial impedance)
* i=1: τ (time constant in seconds)
* i=2: c (offset)
* i=3: R² (goodness of fit)

---

#### custom_fit

Pixel-wise regression with arbitrary user-defined function.

**Algorithm:**
1. Identify data within valid_range
2. Parse function string and perform curve fitting
3. Apply bounds constraints
4. Calculate R² for goodness of fit

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `valid_range` | list | [0.25, 0.75] | Data range as fraction of tidal change |
| `func_str` | str | - | Function string, e.g., `'z0 * np.exp(-t/tau) + c'` |
| `variables` | list | - | Variable names, e.g., `['t', 'z0', 'tau', 'c']` |
| `bounds` | tuple | (-∞, ∞) | Parameter bounds for scipy.optimize.curve_fit |
| `p0` | array | None | Initial parameter guesses |

**Input:** Interval with breath flags

**Output:** Numpy array with dimensions (breaths, vd, rl, n+1)

**Index i:**
* i=0 to n-1: Fitted parameters (in order of variables, excluding dependent variable)
* i=n: R² (goodness of fit)

**Note:** Symmetrical quantities (e.g., in bi-exponential regression) are not sorted!

---

### Utility Items

#### thorax_roi

Returns a predefined thorax region of interest mask.

**Parameters:** None

**Input:** None (standalone)

**Output:** Boolean numpy array with dimensions (1, 32, 32)

---

#### interval_data

Access interval attributes.

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `attr_name` | str | Name of interval attribute to access |

**Input:** Interval

**Output:** Depends on attribute

---

#### detector_data

Access detector attributes or settings.

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `attr_name` | str | Name of detector attribute or setting |

**Input:** Interval with detector

**Output:** Detector attribute value or setting value

---

#### passthrough

Output data as-is (after preprocessing).

**Parameters:** None

**Input:** Interval data

**Output:** Numpy array with dimensions (samples, vd, rl)

---

## Operations

Operations transform data from prerequisites. All operations ignore non-numeric (NaN) data.

### Aggregation Operations

#### mean_over_time / median_over_time

Calculate mean or median over the time dimension.

**Input:** Numpy array with dimensions (t, vd, rl [, i])

**Parameters:** None (standard numpy axis=0)

**Output:** Numpy array with dimensions (1, vd, rl [, i])

---

#### mean_over_image / median_over_image / sum_over_image / percentile_over_image

Calculate statistics over spatial dimensions (vd, rl).

**Input:** Numpy array with dimensions (t, vd, rl [, i])

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `percentile` | float | - | For percentile_over_image only |

**Output:** Numpy array with dimensions (t, 1, 1 [, i])

---

### Thresholding Operations

#### threshold

Binary threshold function.

**Algorithm:** Value > threshold → 1, else → 0

**Input:** Numpy array with dimensions (t, vd, rl)

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `threshold` | float | Absolute threshold value |
| `invert` | bool | If True, reverse comparison |

**Output:** Boolean array with dimensions (t, vd, rl)

---

#### percentile

Relative threshold based on data percentile.

**Algorithm:** Calculate percentile value from data, then threshold above it.

**Input:** Numpy array with dimensions (t, vd, rl)

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `percentile` | float | Percentile value (0-100) |

**Output:** Boolean array with dimensions (t, vd, rl)

---

#### cumulative_threshold

Select minimal number of pixels that sum to threshold.

**Algorithm:**
1. Sort pixels by value (descending)
2. Accumulate until threshold is reached
3. Set selected pixels to 1, others to 0

**Input:** Numpy array with dimensions (1, vd, rl)

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `threshold` | float | Target cumulative sum |

**Output:** Boolean array with dimensions (1, vd, rl)

---

### Statistical Operations

#### sum / std

Sum or standard deviation over specified axes.

**Input:** Numpy array (any dimensions)

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `axis` | int/tuple | Axis or axes to reduce (see numpy.nansum/nanstd) |

**Output:** Numpy array with reduced dimensions

---

#### coefficient_of_variation

Standard deviation divided by absolute mean.

**Input:** Numpy array with dimensions (t, vd, rl)

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `axis` | tuple | (1, 2) | Axes to apply operation over |

**Output:** Scalar or array with dimensions (t,) depending on axis

---

#### global_inhomogeneity_index

Global Inhomogeneity Index (GI) - measure of impedance distribution heterogeneity.

**Algorithm:** GI = Σ|zᵢ - median(z)| / Σzᵢ

**Input:** Numpy array with dimensions (t, vd, rl)

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `axis` | tuple | (1, 2) | Axes to apply operation over |

**Output:** Array with reduced dimensions

---

#### min / max

Minimum or maximum value.

**Input:** Numpy array (any dimensions)

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `axis` | int/tuple | Axis or axes to reduce (see numpy.nanmin/nanmax) |

**Output:** Scalar or reduced array

---

### Normalization Operations

#### normalize_sum

Normalize data to sum to 1.

**Input:** Numpy array with dimensions (1, vd, rl)

**Parameters:** None

**Output:** Numpy array with dimensions (1, vd, rl), sum = 1

---

#### normalize_max

Normalize data to maximum value of 1.

**Input:** Numpy array with dimensions (1, vd, rl)

**Parameters:** None

**Output:** Numpy array with dimensions (1, vd, rl), max = 1

---

### Spatial Operations

#### centroid

Calculate centroid of distribution along spatial dimensions.

**Algorithm:**
1. Calculate weighted mean position
2. Normalize to fraction of dimension size (0-1)

**Input:** Numpy array with dimensions (t, vd, rl)

**Parameters:** None

**Output:** Numpy array with dimensions (t, 1, 1, 2)

**Index i:**
* i=0: Centroid in ventral-dorsal direction (fraction, 0-1)
* i=1: Centroid in right-left direction (fraction, 0-1)

---

### Mathematical Operations

#### multiply / divide / subtract / add

Element-wise arithmetic operations.

**Input:** 2 numpy arrays with compatible dimensions (broadcasting applies)

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `targets` | list | [None] | Prerequisite indices to operate on |

**Output:** Numpy array with broadcasted dimensions

---

#### minimum / maximum

Element-wise minimum or maximum with a reference value.

**Input:** Numpy array (any dimensions)

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `x2` | float | Reference value for clipping |

**Output:** Numpy array with same dimensions (values clipped to x2)

---

### Data Manipulation Operations

#### slice_last

Select a single element from the last dimension.

**Use case:** Extract specific result from multi-output base items (e.g., τ from expiratory_time_constant).

**Input:** Numpy array with dimensions (..., i)

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `index` | int | Index of result (0-based, supports negative indexing) |

**Output:** Numpy array with dimensions (...)

---

#### slice_first

Select a single element from the first (time) dimension.

**Input:** Numpy array with dimensions (t, ...)

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `index` | int | Index of time slice (0-based, supports negative indexing) |

**Output:** Numpy array with dimensions (1, ...)

---

#### stack

Stack multiple arrays along the first (time) dimension.

**Input:** Multiple numpy arrays with compatible spatial dimensions

**Parameters:** None

**Output:** Numpy array with dimensions (n, ...) where n = number of inputs

---

#### size

Return the number of elements along the first dimension.

**Input:** Numpy array

**Parameters:** None

**Output:** Integer (length of first dimension)

---

#### apply_mask

Apply boolean mask, turning masked values to NaN.

**Input:** 
* Data array with dimensions (t, vd, rl)
* Boolean mask array with matching spatial dimensions

**Parameters:** None

**Output:** Numpy array with same dimensions (masked values → NaN)

---

#### costa_approach

Costa et al. approach for PEEP titration analysis. Calculates loss of compliance compared to maximum compliance in both directions.

**Algorithm:**
1. Find maximum compliance pixel-wise
2. Calculate normalized difference from maximum
3. Separate into collapse (before max) and overdistension (after max)

**Input:** Numpy array with dimensions (n, vd, rl) - stack of tidal images

**Parameters:** None

**Output:** Numpy array with dimensions (n, vd, rl, 2)

**Index i:**
* i=0: Collapse (loss of compliance before maximum)
* i=1: Overdistension (loss of compliance after maximum)

---

## Analysis Item Structure

### JSON Template Format

Analysis items are defined in JSON format:

```json
{
  "name": "unique_item_name",
  "title": "Display title",
  "unit": "Unit string (optional)",
  "identifier": "LOINC/SNOMED CT code (optional)",
  "comment": "Description (optional)",
  
  // For base items:
  "base_item": "base_item_name",
  "parameters": { "param1": value1, ... },
  "interval": "interval_name",
  "preprocessors": [
    { "name": "preprocessor_name", "parameters": { ... } }
  ],
  
  // For custom items:
  "prerequisites": [
    { "name": "prerequisite_name" },
    { "base_item": "base_item_name", "interval": "interval_name", ... }
  ],
  "operations": [
    { "name": "operation_name", "parameters": { ... }, "targets": [0, 1] }
  ],
  
  "export": true
}
```

### Field Descriptions

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `name` | str | Yes | Unique identifier (no special chars except `_`) |
| `title` | str | No | Human-readable display name |
| `unit` | str | No | Unit of measurement |
| `identifier` | str | No | Standard identifier (LOINC, SNOMED CT) |
| `comment` | str | No | Additional notes |
| `base_item` | str | Base items only | Name of base item class |
| `parameters` | dict | No | Parameters for base item/operation |
| `interval` | str | Base items only | Name of interval to use |
| `preprocessors` | list | No | List of preprocessors to apply |
| `prerequisites` | list | Custom items only | Input items for calculation |
| `operations` | list | No | List of operations to apply |
| `export` | bool | No | Whether to export result |

### Processing Flow

1. **Base items**: Load interval data → Apply preprocessors → Calculate base item
2. **Custom items**: Resolve prerequisites → Apply operations in sequence → Result
3. **Export**: Items with `export: true` are saved to output directory

---

## References

* **GREIT Algorithm**: Schweiger et al., "GREIT: a unified approach for EIT reconstruction", 2007
* **Costa Approach**: Costa et al., "Bedside estimation of PEEP-induced lung recruitment", 2019
