"""
    CLOUD18_calibration
A Julia package for calibrating CLOUD18 data.
"""
module CLOUD18_calibration

include("humidity_dependence_calibration.jl")
using .analyze_humidity_dependence
export analyze_humidity_dependence

end
