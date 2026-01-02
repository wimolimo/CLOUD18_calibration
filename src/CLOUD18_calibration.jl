"""
    CLOUD18_calibration
A Julia package for calibrating CLOUD18 data.
"""
module CLOUD18_calibration

include("humidity_dependence_calibration.jl")
using .humidity_dependence_analysis
export humidity_dependence_analysis

end
