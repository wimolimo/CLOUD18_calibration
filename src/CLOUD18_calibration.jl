"""
    CLOUD18_calibration
A Julia package for calibrating CLOUD18 data.
"""
module CLOUD18_calibration

include("calibrate_traces.jl")
using .CalibrateTraces: calibrate_traces_main
export CalibrateTraces, calibrate_traces_main

include("InletLossCorrection.jl")
using .InletLossCorrection: run_inlet_loss_correction
export InletLossCorrection, run_inlet_loss_correction

end # module CLOUD18_calibration

