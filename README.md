# CLOUD18_calibration

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://wimolimo.github.io/CLOUD18_calibration.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wimolimo.github.io/CLOUD18_calibration.jl/dev/)
[![Build Status](https://github.com/wimolimo/CLOUD18_calibration.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/wimolimo/CLOUD18_calibration.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/wimolimo/CLOUD18_calibration.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/wimolimo/CLOUD18_calibration.jl)

A Julia package for calibrating Proton Transfer Reaction Time-of-Flight Mass Spectrometry (PTR-ToF-MS) data from the CLOUD18 campaign. This package provides tools for both dry and humidity-dependent calibration of organic compound traces, normalizing to primary ions, using hexanone as reference. Also, it is corrected for duty-cycle and inlet-loss.

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Calibration Workflow](#calibration-workflow)
- [Configuration](#configuration)
- [Usage Examples](#usage-examples)
- [Output Files](#output-files)
- [Calibration Methods](#calibration-methods)
- [Troubleshooting](#troubleshooting)
- [Dependencies](#dependencies)
- [Authors](#authors)

## Overview

CLOUD18_calibration processes PTR-ToF-MS measurement data by:
1. Normalizing signals to primary ion intensities (e.g. NH4⁺ or the full primary ion list of NH4+ and clusters)
2. Applying dry calibration factors from standard gas measurements
3. Optionally applying humidity-dependent calibration for compounds sensitive to water vapor
4. Converting raw detector signals (dcps) to concentrations (ppb)

The package uses hexanone (C6H12O) as the primary reference compound and supports calibration of hundreds of organic compounds simultaneously.

## Features

- **Flexible Calibration**: Choose between dry-only or humidity-dependent calibration
- **Primary Ion Normalization**: Automatically corrects for ionization source variations
- **Duty Cycle Correction**: Accounts for mass-dependent detection efficiency
- **Composition-Based Calibration**: Different strategies for compounds based on oxygen content:
  - 0 O: No calibration (zero sensitivity)
  - 1 O: Humidity-dependent calibration (equilibrium)
  - ≥2 O: Dry calibration (kinetic limit)
- **Interactive Selection of included Data Points**: Manually exclude outliers from calibration fits
- **Interactive Selection of Fit Functions**: Choose between Power or linear fit (soon exponential and double exponential)
- **Uncertainty Propagation**: Error estimates for calibration factors
- **Export to CLOUD Format**: Direct export of calibrated traces with metadata headers

## Installation

### Prerequisites
- Julia ≥ 1.11
- TOFTracer2 package (included in sources)
- Python with matplotlib (for plotting via PyPlot)

### Setup

1. Clone the repository:
```bash
git clone https://github.com/wimolimo/CLOUD18_calibration.jl.git
cd CLOUD18_calibration
```

2. Activate and instantiate the Julia environment:
```julia
julia> ]  # Enter package mode
pkg> activate .
pkg> instantiate
pkg> precompile
```

3. Verify TOFTracer2 dependency path in `Project.toml`:
```toml
[sources]
TOFTracer2 = {path = "..\\TOF-Tracer2-dev"}
```

## Calibration Workflow

### Complete 3-Step Process

#### Step 1: Humidity Dependence Calibration (Optional)
Generate humidity-dependent parameters from controlled standard gas measurements at varying absolute humidity:

```julia
# Run first to create humcalibfile
include("src/humidity_dependence_calibration.jl")
```

**Output**: `fitParameters_relative.txt` with double exponential fit parameters (p1-p5) for each compound

#### Step 2: Trace Calibration (Main Step)
Calibrate all measurement data using dry and humidity-dependent parameters:

```julia
include("src/run_CLOUD18_calibration.jl")
```

The script will guide you through:
- Selecting primary ions (full list or NH4+ only)
- Choosing humidity-dependent vs. dry calibration
- Selecting calibration data points interactively
- Choosing fit function type for reference vs PI signal (linear, power)

**Inputs**:
- Dry calibration data (HDF5 or CSV)
- Humidity calibration parameters (from Step 1)
- Measurement result files (HDF5)
- Licor humidity data (TXT files)

**Outputs**:
- Calibration factor traces (PNG/PDF plots)
- Directly calibrated concentration traces (PNG/PDF plots)
- Exported traces (CSV in CLOUD format)

#### Step 3: Inlet Loss Correction
Apply transmission corrections for sampling line losses.
##########################################################################################################

## Configuration

Edit `src/run_CLOUD18_calibration.jl` to configure your calibration:

### File Paths
```julia
# Directory structure
dir_CLOUD18 = joinpath(@__DIR__, "..", "..")
dir_calib_data = joinpath(dir_CLOUD18, "CLOUD18_data", "Calibration")
dir_licor_data = joinpath(dir_CLOUD18, "CLOUD18_data", "Licor")
#rename licor file to use with licor_restofthename.txt

# Calibration data files
drycalibsfile = joinpath(dir_calib_data, "dry_std", "results", "_result.hdf5")
humcalibfile = joinpath(dir_calib_data, "Humidity-dependent_std", "results", "fitParameters_relative.txt")

# Measurement data to calibrate
resultfp = joinpath(dir_CLOUD18, "CLOUD18_data", "Nonanal", "2025-11-25")
resultfiles = ["$(resultfp)/results/_result.hdf5"]
```

### Calibration Parameters
```julia
ionization = "NH4+"  # or "H+"
refMass = massLibrary.HEXANONE_nh4[1]  # Reference compound mass
refName = "C6H12O"  # Reference compound formula ############################################################################
```

### Export Settings
```julia
exportTraces = true

HeaderForExportDict = Dict(
    "title" => "CLOUD18 Nonanal Calibration",
    "level" => 2,
    "version" => "01",
    "authorname_mail" => "Your Name your.email@example.com",
    "units" => "ppt",
    "addcomment" => "Humidity-dependent calibration with hexanone reference...",
    "threshold" => 0,
    "nrrows_addcomment" => 4
)
```

## Usage Examples

### Example 1: Full Humidity-Dependent Calibration
```julia
# 1. Set up configuration
config = CalibrationConfig(
    dir_licor_data, 
    humcalibfile, 
    drycalibsfile, 
    resultfp, 
    resultfiles, 
    "NH4+",
    [],  # Primary ions selected interactively
    refMass, 
    refName, 
    true,  # Export traces
    HeaderForExportDict
)

# 2. Run calibration
calibrate_traces_main(config)

# 3. When prompted:
# - Choose 'f' for full primary ion list
# - Choose 'y' for humidity-dependent calibration
# - Select outliers to exclude from fit (click on plot and press 'c')
# - Choose 'power' or 'linear' fit function
```

### Example 2: Dry Calibration Only
```julia
# Same as above, but when prompted:
# - Choose 'n' for dry calibration only (no humidity dependence)
```

### Example 3: Multiple Result Files
```julia
# Calibrate multiple measurements with same mass list
resultfiles = [
    "$(resultfp)/part1/results/_result.hdf5",
    "$(resultfp)/part2/results/_result.hdf5",
    "$(resultfp)/part3/results/_result.hdf5"
]
```

## Output Files

### Generated in Dry Calibration Directory
- `dryCalibs.png` - Primary ions and reference mass time series
- `Hexanone_VS_PIs.png/pdf` - Dry calibration fit plot with uncertainty band
- `Hexanone_VS_PIs_params.csv` - Fit parameters for reuse

### Generated in Result Directory
- `CalibrationTraces.png/pdf` - Time series of calibration factors for each compound
- `DirectlyCalibratedTraces.png/pdf` - Calibrated concentration time series
- `CLOUD_PTR_NH4_ambient_vXX.txt` - Exported calibrated traces (if exportTraces=true)

### File Formats
- **HDF5**: Binary format for measurement results (masses, traces, timestamps)
- **CSV/TXT**: Tabular calibration parameters and exported data
- **PNG/PDF**: Publication-ready plots with legends and axis labels

## Calibration Methods

### Dry Calibration (Hexanone vs Primary Ions)
Normalizes hexanone signal to primary ion intensity using fit functions:
- **Power function**: `f(x) = a * x^b` (recommended for most cases)
- **Linear function**: `f(x) = a * x + b`
- **Exponential**: `f(x) = a * exp(-b*x) + c` (not implemented yet)
- **Double exponential**: `f(x) = a₁ * exp(-b₁*x) + a₂ * exp(-b₂*x) + c` (not implemented yet)

### Humidity-Dependent Calibration
For compounds with 1 oxygen atom and compounds from the gas standard, applies additionally humidity correction:
```
Sensitivity(AH) = f_dry * f_humid(AH)
```
where `f_humid` is a double exponential function of absolute humidity (AH), normalized to the zero humidity point of hexanone.

### Composition-Based Strategy
- **Undefined masses**: Dry calibration with hexanone
- **0 oxygen**: Not calibrated (zero sensitivity assumed)
- **1 oxygen**: Humidity-dependent calibration (equilibrium regime)
- **≥2 oxygen**: Dry calibration (kinetic limit regime)
- **Gas standard compounds**: Individual humidity-dependent calibration

### Uncertainty Estimation
Propagates errors from:
- Fit parameter uncertainties of the reference vs PI fit (standard errors)
- Humidity calibration scatter (not implemented yet)

## Troubleshooting

### Common Issues

**Error: "Reference compound not found in calibration data"**
- Check that `refName` matches a formula in `humcalibfile` (e.g., "C6H12O.NH4+" or "C6H12O")
- The package uses atomic composition matching, so "C6H12O.NH4+" = "C6H12H+O.NH4+"

**Plot window not responding**
- Ensure PyPlot is properly configured with matplotlib backend
- Try closing all figures and restarting Julia
- Check X11 forwarding if working remotely

**No calibration data exported**
- Verify `exportTraces = true` in configuration
- Check that filtered traces meet export criteria (C≥1, N=1)
- Adjust `findVaryingMasses` parameters if too few traces selected

**Licor data interpolation fails**
- Ensure Licor files match naming pattern `licor_*.txt`
- Check that Licor timestamps overlap with measurement period
- Verify column name: `H₂O_(mmol_mol⁻¹)` in header

### Interactive Selection Tips
- Zoom and pan before selecting points
- Press 'c' on keyboard while hovering over outliers
- Select exactly the number of points requested
- Program continues on automatcally after selecting the number of points requested

## Dependencies

### Julia Packages
- **TOFTracer2**: Core mass spectrometry data processing
- **CSV, DataFrames**: Tabular data handling
- **HDF5**: Binary result file format
- **PyPlot/PythonPlot**: Plotting via matplotlib
- **LsqFit**: Nonlinear curve fitting
- **Statistics, Dates**: Standard library utilities

### External
- Python ≥ 3.7 with matplotlib
- Licor data acquisition system (optional, for humidity-dependent calibration)

## Authors

- **Timo Wittler** - wittler.timo@uibk.ac.at
- **Clea Ruth** - clea.ruth@student.uibk.ac.at

Institute of Ion Physics and Applied Physics  
University of Innsbruck, Austria

## License

See [LICENSE](LICENSE) file for details.

## Citation

If you use this package in your research, please cite:
```
Wittler, T., Ruth, C. (2026). CLOUD18_calibration: Calibration tools for PTR-ToF-MS 
data from the CLOUD18 campaign. https://github.com/wimolimo/CLOUD18_calibration.jl
```

## Related Projects

- [TOFTracer2](https://github.com/wimolimo/TOF-Tracer2-dev): Core PTR-ToF-MS data processing library
- [CLOUD Experiment](https://cloud.web.cern.ch/): CERN's Cosmics Leaving OUtdoor Droplets facility

---

For questions, issues, or contributions, please open an issue on GitHub or contact the authors directly.
