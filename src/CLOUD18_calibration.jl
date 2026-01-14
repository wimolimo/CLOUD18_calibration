"""
    CLOUD18_calibration
A Julia package for calibrating CLOUD18 data.
"""
module CLOUD18_calibration

include("calibrate_traces.jl")
using .CalibrateTraces: calibrate_traces_main, CalibrationConfig
export CalibrateTraces, calibrate_traces_main, CalibrationConfig

end # module CLOUD18_calibration

