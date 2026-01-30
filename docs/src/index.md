```@meta
CurrentModule = CLOUD18_calibration
```

# CLOUD18_calibration

Documentation for [CLOUD18_calibration](https://github.com/wimolimo/CLOUD18_calibration.jl).

## Main Calibration Function

```@docs
CLOUD18_calibration.CalibrateTraces.calibrate_traces_main
```

## Building Calibration

```@docs
CLOUD18_calibration.CalibrateTraces.build_calibration_traces
CLOUD18_calibration.CalibrateTraces.calc_fhex
CLOUD18_calibration.CalibrateTraces.compute_summed_primary_ions
```

## Data Loading & Processing

```@docs
CLOUD18_calibration.CalibrateTraces.load_hexVSpis_params
CLOUD18_calibration.CalibrateTraces.load_and_merge_results
CLOUD18_calibration.CalibrateTraces.load_licor_data
CLOUD18_calibration.CalibrateTraces.interpolate_licor_to_ptr_time
```

## Plotting

```@docs
CLOUD18_calibration.CalibrateTraces.plot_calibration_traces
CLOUD18_calibration.CalibrateTraces.plot_directly_calibrated_traces
CLOUD18_calibration.CalibrateTraces.scatterDryCalibs2
```

## Helper Functions

```@docs
CLOUD18_calibration.CalibrateTraces.parse_formula_to_composition
CLOUD18_calibration.CalibrateTraces.find_formula_index
CLOUD18_calibration.CalibrateTraces.plot_fit
CLOUD18_calibration.CalibrateTraces.dryCal_selectPIandRefDataInteractive
CLOUD18_calibration.CalibrateTraces.export_calibrated_traces
```

## Humidity Dependence Calibration

```@docs
CLOUD18_calibration.HumidityDependenceCalibration.getHumiditySensitivity
CLOUD18_calibration.HumidityDependenceCalibration.DoubleExponential_and_fit
CLOUD18_calibration.HumidityDependenceCalibration.get_ion_metadata
CLOUD18_calibration.HumidityDependenceCalibration.export_sensitivities
CLOUD18_calibration.HumidityDependenceCalibration.print_relative_error_summary
CLOUD18_calibration.HumidityDependenceCalibration.plot_relative_normalization
```

## Inlet Loss Correction

```@docs
CLOUD18_calibration.InletLossCorrection.run_inlet_loss_correction
```
