"""
    CLOUD18_calibration
A Julia package for calibrating CLOUD18 data.
"""
module CLOUD18_calibration

include("calibrate_traces.jl")
export CalibrateTraces, calibrate_traces_main

end
