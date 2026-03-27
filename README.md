# TMD R Translation

This is a working R translation of the MATLAB TMD3.0 code. Being mainly interested in the predicted tides, the translation focuses on the core functionality for that, and does not include all of the analysis and plotting utilities in the original MATLAB code.

To test it, a netcdf file is needed that is compatible with the TMD3.0 format. See also https://nl.mathworks.com/matlabcentral/fileexchange/133417-tide-model-driver-tmd-version-3-0.

If Rcpp is available, to attemps speedup, but the code falls back to pure R if Rcpp is not available or if compilation fails. 

To test it, clone this repository, and run:

```r
source("load_tmd_translation_current.R")
```

## Included functionality

- `tmd_data()`
- `tmd_interp()`
- `tmd_predict()`
- `tmd_conlist()`
- `tmd_ellipse()`
- `tidal_range()`
- `tmd_astrol()`
- `tmd_constit()`
- `tmd_nodal()`
- `tmd_harp()`
- `tmd_infer_minor()`
- `tmd_ll2ps()`
- `tmd_ps2ll()`

## Using `tmd_predict()`

`tmd_predict()` is the main function for computing predicted tides from a TMD-compatible NetCDF file.

Basic form:

```r
z <- tmd_predict(filename, lat, lon, time)
```

Main arguments:

- `filename`: path to the TMD-compatible `.nc` file
- `lat`, `lon`: latitude and longitude locations
- `time`: a date/time input accepted by the translation
- `ptype`: optional solution type, one of `"h"`, `"z"`, `"u"`, `"U"`, `"v"`, `"V"`
- `input_mode`: either `"paired"` or `"grid"`
- `time_mode`: either `"all"` or `"track"`

### Single-point example

```r
z <- tmd_predict(
  "Subset_NS_EC_forTMD30_21_2_2025.nc",
  52.5,
  2.95,
  as.Date("2025-01-01")
)
```

### Paired points example

Use `input_mode = "paired"` when latitude and longitude vectors should be matched element by element:

```r
z <- tmd_predict(
  "Subset_NS_EC_forTMD30_21_2_2025.nc",
  c(52.0, 52.5, 53.0),
  c(2.8, 2.9, 3.0),
  as.Date("2025-01-01"),
  input_mode = "paired"
)
```

### Grid example

Use `input_mode = "grid"` when latitude and longitude vectors should be expanded into a full grid:

```r
z <- tmd_predict(
  "Subset_NS_EC_forTMD30_21_2_2025.nc",
  seq(52, 55, 0.1),
  seq(2, 4, 0.1),
  as.Date("2025-01-01"),
  input_mode = "grid"
)
```

### Multiple times for multiple locations

By default, `time_mode = "all"`, which evaluates all locations for all times.

Example:

```r
z <- tmd_predict(
  "Subset_NS_EC_forTMD30_21_2_2025.nc",
  c(52.5, 52.6),
  c(2.95, 3.00),
  c(
    as.POSIXct("2025-01-01 00:00:00", tz = "UTC"),
    as.POSIXct("2025-01-01 06:00:00", tz = "UTC")
  ),
  input_mode = "paired",
  time_mode = "all"
)
```

This returns all location/time combinations, i.e. a `2 x 2` result in this example.

### Drift-track example

Use `time_mode = "track"` when each location should be paired with the corresponding time:

```r
z <- tmd_predict(
  "Subset_NS_EC_forTMD30_21_2_2025.nc",
  c(52.5, 52.6),
  c(2.95, 3.00),
  c(
    as.POSIXct("2025-01-01 00:00:00", tz = "UTC"),
    as.POSIXct("2025-01-01 06:00:00", tz = "UTC")
  ),
  input_mode = "paired",
  time_mode = "track"
)
```

This returns one value per `(lat[i], lon[i], time[i])` combination.

### Current velocity example

For zonal or meridional velocity predictions, specify `ptype`:

```r
u <- tmd_predict(
  "Subset_NS_EC_forTMD30_21_2_2025.nc",
  52.5,
  2.95,
  as.Date("2025-01-01"),
  ptype = "u"
)
```

Notes:

- `ptype = "h"` is the default.
- `ptype = "z"` is treated the same as `"h"`.
- `input_mode = "paired"` is the default.
- `time_mode = "all"` is the default.

## Input mode

The public functions support:

- `input_mode = "paired"`
- `input_mode = "grid"`

This applies to:

- `tmd_predict()`
- `tmd_interp()`
- `tmd_ellipse()`

### `input_mode = "paired"`

Latitude and longitude vectors are interpreted elementwise.

### `input_mode = "grid"`

Latitude and longitude vectors are expanded into a full grid.

The default is `input_mode = "paired"`.

## Time mode

`tmd_predict()` supports:

- `time_mode = "all"`
- `time_mode = "track"`

### `time_mode = "all"`

All locations are evaluated for all times.

For example, 2 locations and 3 times produce 6 predicted values arranged as a location-by-time result.

### `time_mode = "track"`

Each location is paired with the corresponding time.

For example, 2 locations and 2 times produce 2 predicted values:

- `(lat[1], lon[1], time[1])`
- `(lat[2], lon[2], time[2])`

The default is `time_mode = "all"`.

## Optional Rcpp acceleration

This version can use `Rcpp` to speed up the hottest numerical kernels:

- multi-slice bilinear interpolation
- harmonic synthesis in `tmd_harp()`

If `Rcpp` is installed and the C++ file can be compiled, the accelerated path is used automatically.

If `Rcpp` is not installed, or compilation fails, the code falls back to the working pure-R implementation automatically.

The C++ backend file is:

- `tmd_rcpp_backend.cpp`

The R wrapper is:

- `tmd_rcpp_overrides.R`

## Benchmarking

You can benchmark the current implementation with:

```r
source("M:/TMD3.0 translation/benchmark_tmd_translation_current.R")
benchmark_tmd_translation()
```

This compares:

- the pure-R loader
- the Rcpp-capable loader

and reports:

- min / median / mean runtime
- median speedup
- whether the outputs match numerically
- the maximum absolute difference between outputs

You can increase the number of repetitions, for example:

```r
benchmark_tmd_translation(
  filename = "Subset_NS_EC_forTMD30_21_2_2025.nc",
  date = as.Date("2025-01-01"),
  n_iter = 5
)
```

## Main implementation files

- `load_tmd_translation_current.R`
- `benchmark_tmd_translation_v2.R`
- `tmd_helpers.R`
- `tmd_data_interp.R`
- `tmd_harmonics.R`
- `tmd_analysis.R`
- `tmd_final_overrides.R`
- `tmd_input_mode_overrides.R`
- `tmd_minor_overrides.R`
- `tmd_predict_overrides.R`
- `tmd_time_mode_overrides.R`
- `tmd_rcpp_overrides.R`
- `tmd_rcpp_overrides_v2.R`
- `tmd_rcpp_backend.cpp`
