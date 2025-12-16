module analyze_humidity_dependence

using HDF5
#import PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PythonPlot
using Dates
using CSV
using DataFrames
import LsqFit
using TOFTracer2
import TOFTracer2.InterpolationFunctions as IntpF
import TOFTracer2.CalibrationFunctions as CalF
import TOFTracer2.ExportFunctions as ExpF
import TOFTracer2.ImportFunctions as ImpF

export analyze_humidity_dependence

dirname = joinpath(@__DIR__, "..", "..", "CLOUD18_data", "Licor")
filename = joinpath(dirname, "2025-11-21.txt")

bg_dirname = joinpath(@__DIR__, "..", "..", "CLOUD18_data", "Calibration", "humidity_dependent_BG")
bg_filename = joinpath(bg_dirname, "")



function analyze_humidity_dependence()

end

end # module analyze_humidity_dependence