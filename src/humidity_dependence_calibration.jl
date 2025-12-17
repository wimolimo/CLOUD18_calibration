"""
Analyze humidity dependence of calibration data from CLOUD18 Licor measurements.
"""
module analyze_humidity_dependence

using HDF5
# import PyCall
# pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using Plots
using Dates
using CSV
using DataFrames
# import LsqFit
using TOFTracer2
import TOFTracer2.InterpolationFunctions as IntpF
import TOFTracer2.CalibrationFunctions as CalF
import TOFTracer2.ExportFunctions as ExpF
import TOFTracer2.ImportFunctions as ImpF
import TOFTracer2.massLibrary as massLibrary

export analyze_humidity_dependence

plotStart = DateTime(2000, 1, 1, 0, 0, 0)
plotEnd = DateTime(3000, 1, 1, 0, 0, 0)

dirname = joinpath(@__DIR__, "..", "..", "CLOUD18_data", "Licor")
filename = joinpath(dirname, "2025-11-21.txt")

bg_dirname = joinpath(@__DIR__, "..", "..", "CLOUD18_data", "Calibration", "humidity_dependent_BG")
bg_filename = joinpath(bg_dirname, "")

ions2plot = "NH4+" # "NH4+" # "all", "NH4+", "H+"
STD_masses_dict = massLibrary.CLOUD_brownSTD_masses

####################################
# select masses and ions to analyze
####################################

massesToPlot = []
keysToPlot = []
if ions2plot == "NH4+"
    for key in keys(STD_masses_dict)
        append!(massesToPlot, STD_masses_dict[key][1][2])
        push!(keysToPlot, key)
    end
    ion = "NH4+"
elseif ions2plot == "H+"
    for key in keys(STD_masses_dict)
        append!(massesToPlot, STD_masses_dict[key][1][1])
        push!(keysToPlot, key)
    end
    ion = "NH4+"
elseif ions2plot == "all"
    for key in ["TMB"] # you choose, which
        append!(massesToPlot, STD_masses_dict[key][1])
    end
    ion = "NH4+"
end





function humidity_dependence_analysis()

end

end # module analyze_humidity_dependence