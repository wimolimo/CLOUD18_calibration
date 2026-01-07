"""
    CLOUD18_calibration
A Julia package for calibrating CLOUD18 data.
"""
module CLOUD18_calibration

include("humidity_dependence_calibration.jl")
include("calibrate_traces.jl")
using .HumidityDependence, .CalibrateTraces
export HumidityDependence, CalibrateTraces, calibrate_traces_main

end
