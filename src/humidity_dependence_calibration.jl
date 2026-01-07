"""
Analyze humidity dependence of calibration data from CLOUD18 Licor measurements.
(partially adapted from analyze_humidityDependentCalibrations.jl)
"""
module HumidityDependence

using HDF5
# import PyCall
# pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx

# Force matplotlib to use a headless backend before PyCall / Plots load
# ENV["MPLBACKEND"] = "agg"
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

function analyze_humidity_dependence()

plotStart = DateTime(2000, 1, 1, 0, 0, 0)
plotEnd = DateTime(3000, 1, 1, 0, 0, 0)

# humidity dependent std files
std_dirname = joinpath(@__DIR__, "..", "..", "CLOUD18_data", "Calibration", "Humidity-dependent_std", "results")
std_filename = joinpath(std_dirname, "_result.hdf5")

# humidity measurements (need correct time frame!!!!!)
licor_dirname = joinpath(@__DIR__, "..", "..", "CLOUD18_data", "Licor")
licor_filename = joinpath(licor_dirname, "2025-11-21.txt")

# background measurements
bg_dirname = joinpath(@__DIR__, "..", "..", "CLOUD18_data", "Calibration", "humidity_dependent_BG")
bg_filename = joinpath(bg_dirname, "2025-11-24 01h00m59_1.5ppth_1ppb_std_brown.h5")

println("load data")
ions2plot = "NH4+" # "NH4+" # "all", "NH4+", "H+"
STD_masses_dict = massLibrary.CLOUD_brownSTD_masses

####################################
# select masses and ions to analyze
####################################

println("select masses and ions to analyze")

massesToPlot = []
keysToPlot = []
print(keys(STD_masses_dict))

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
    ion = "H+"
elseif ions2plot == "all"
    for key in ["TMB"] # you choose, which
        append!(massesToPlot, STD_masses_dict[key][1])
    end
    ion = "NH4+"
end

##################################
# plot raw data and select filters
##################################

println("plot raw data and select filters")

tracesFig, tracesAx, measResult = PlotFunctions.plotTracesFromHDF5(std_filename, massesToPlot;
    plotHighTimeRes=false,
    smoothing=10,
    backgroundSubstractionMode=0,
    timedelay=Dates.Hour(0),
    isobarToPlot=0,
    plotsymbol=".-",
    timeFrame2plot=(plotStart, plotEnd),
    ion=ion
)
savefig("$(std_filename)Traces_$(ions2plot).png")

humDat = PlotFunctions.load_plotLicorData(licor_filename; ax=tracesAx)
tracesFig.tight_layout()

# if bgTimes, humidityLimits and signalTimes is not globally defined, define them interactively here
if !isdefined(Main, :bgTimes)
    println("do the next part interactively:\n")
    IFIG = PlotFunctions.InteractivePlot(std_filename, tracesAx)
    (humidityLimits, bgTimes, signalTimes) = CalF.humCal_getDatalimitsFromPlot(IFIG)
end

if (length(IFIG.deleteXlim) > 0) && (length(IFIG.deleteXlim) % 2 == 0)
    for i in collect(1:2:length(IFIG.deleteXlim))
        measResult.Traces[IFIG.deleteXlim[i].<=measResult.Times.<=IFIG.deleteXlim[i+1], :] .= NaN
    end
end

#######################################
# calculate and plot calibration points
#######################################

println("calculate and plot calibration points")

(calibData, calibData_std, humidities) = CalF.humcal_getHumidityDependentSensitivity(measResult, humDat;
    hums=humidityLimits,
    bgtimes=bgTimes,
    signaltimes=signalTimes,
    pptInInlet=1000)

fig = figure(figsize=(14, 9))
(calibFig, calibAx) = PlotFunctions.scatter_errorbar(fig, measResult, humidities, calibData, calibData_std; ion=ion)
xlabel("absolute humidity [mmol mol⁻¹]")
ylabel("sensitivity [dcps ppt⁻¹]")

################################
# plot and export fit parameters
################################

println("plot and export fit parameters")

hum4plot = collect(0:0.2:12)
fitParams = []
fitParamErrors = []
colornames = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan", "tab:blue", "tab:orange", "tab:green", "tab:red"]

if length(keysToPlot) == length(measResult.MasslistMasses)
    for (i, m) in enumerate(measResult.MasslistMasses)
        (param, stderror, fitlabel) = CalF.fitParameters_DoubleExponential(humidities, calibData[:, i])
        push!(fitParams, param)
        push!(fitParamErrors, stderror)
        plot(hum4plot, CalF.DoubleExponential(hum4plot, param),
            color=colornames[i],
            label=string(round(measResult.MasslistMasses[i], digits=3), " - ", keysToPlot[i], ".", ions2plot, " -- sens(AH) = $(round(param[1],sigdigits=3)) * exp(-$(round(param[2],sigdigits=3))*AH) + $(round(param[3],sigdigits=3))*exp(-$(round(param[4],sigdigits=3))*AH) + $(round(param[5],sigdigits=3))"))
    end
else
    for (i, m) in enumerate(measResult.MasslistMasses)
        (param, stderror, fitlabel) = CalF.fitParameters_DoubleExponential(humidities, calibData[:, i])
        push!(fitParams, param)
        push!(fitParamErrors, stderror)
        plot(hum4plot, CalF.DoubleExponential(hum4plot, param),
            color=colornames[i],
            label="m/z $(round(m,digits=3)), $(MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i])) -- sens(AH) = $(round(param[1],sigdigits=3)) * exp(-$(round(param[2],sigdigits=3))*AH) + $(round(param[3],sigdigits=3))*exp(-$(round(param[4],sigdigits=3))*AH) + $(round(param[5],sigdigits=3))")
    end
end


legend()
savefig("$(std_filename)calibration_lin_$(ions2plot).png")
yscale("log")
savefig("$(std_filename)calibration_log_$(ions2plot).png")

fitParams2Export = hvcat(length(measResult.MasslistMasses), (fitParams[a][j] for a in 1:length(measResult.MasslistMasses), j in 1:length(fitParams[1]))...)
fitParamErrors2Export = hvcat(length(measResult.MasslistMasses), (fitParamErrors[a][j] for a in 1:length(measResult.MasslistMasses), j in 1:length(fitParamErrors[1]))...)

ExpF.exportFitParameters("$(std_filename)fitParameters.txt", fitParams2Export, fitParamErrors2Export,
    measResult.MasslistMasses, measResult.MasslistCompositions;
    fitfunction="sensitivity(AH) = p1 * exp.(-p2*AH) .+ p3*exp.(-p4*AH) .+ p5")


#######################################
# calculate and plot calibration points
#######################################

println("re-calculate and plot calibration points with bg corrections")

(calibData, calibData_std, humidities) = CalF.humcal_getHumidityDependentSensitivity(measResult, humDat;
    hums=humidityLimits,
    bgtimes=bgTimes,
    signaltimes=signalTimes,
    pptInInlet=1000)

fig = figure(figsize=(14, 9))
(calibFig, calibAx) = PlotFunctions.scatter_errorbar(fig, measResult, humidities, calibData, calibData_std; ion=ion)
xlabel("absolute humidity [mmol mol⁻¹]")
ylabel("sensitivity [dcps ppt⁻¹]")

################################
# plot and export fit parameters
################################

println("plot and export fit parameters")

hum4plot = collect(0:0.2:12)
fitParams = []
fitParamErrors = []
colornames = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan", "tab:blue", "tab:orange", "tab:green", "tab:red"]

if length(keysToPlot) == length(measResult.MasslistMasses)
    for (i, m) in enumerate(measResult.MasslistMasses)
        (param, stderror, fitlabel) = CalF.fitParameters_DoubleExponential(humidities, calibData[:, i])
        push!(fitParams, param)
        push!(fitParamErrors, stderror)
        plot(hum4plot, CalF.DoubleExponential(hum4plot, param),
            color=colornames[i],
            label=string(round(measResult.MasslistMasses[i], digits=3), " - ", keysToPlot[i], ".", ions2plot, " -- sens(AH) = $(round(param[1],sigdigits=3)) * exp(-$(round(param[2],sigdigits=3))*AH) + $(round(param[3],sigdigits=3))*exp(-$(round(param[4],sigdigits=3))*AH) + $(round(param[5],sigdigits=3))"))
    end
else
    for (i, m) in enumerate(measResult.MasslistMasses)
        (param, stderror, fitlabel) = CalF.fitParameters_DoubleExponential(humidities, calibData[:, i])
        push!(fitParams, param)
        push!(fitParamErrors, stderror)
        plot(hum4plot, CalF.DoubleExponential(hum4plot, param),
            color=colornames[i],
            label="m/z $(round(m,digits=3)), $(MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i])) -- sens(AH) = $(round(param[1],sigdigits=3)) * exp(-$(round(param[2],sigdigits=3))*AH) + $(round(param[3],sigdigits=3))*exp(-$(round(param[4],sigdigits=3))*AH) + $(round(param[5],sigdigits=3))")
    end
end

legend()
savefig("$(std_filename)calibration_lin_$(ions2plot).png")
yscale("log")
savefig("$(std_filename)calibration_log_$(ions2plot).png")

fitParams2Export = hvcat(length(measResult.MasslistMasses), (fitParams[a][j] for a in 1:length(measResult.MasslistMasses), j in 1:length(fitParams[1]))...)
fitParamErrors2Export = hvcat(length(measResult.MasslistMasses), (fitParamErrors[a][j] for a in 1:length(measResult.MasslistMasses), j in 1:length(fitParamErrors[1]))...)

ExpF.exportFitParameters("$(std_filename)fitParameters.txt", fitParams2Export, fitParamErrors2Export,
    measResult.MasslistMasses, measResult.MasslistCompositions;
    fitfunction="sensitivity(AH) = p1 * exp.(-p2*AH) .+ p3*exp.(-p4*AH) .+ p5")

##########################################################################
# correct fit parameters relative to Hexanone and save these also to file
##########################################################################

println("calculate and export fit parameters relative to Hexanone")

if round(massLibrary.HEXANONE_nh4[1],digits=3) in round.(measResult.MasslistMasses,digits=3)
    fitParamsHex = fitParams2Export[:, isapprox.(measResult.MasslistMasses, massLibrary.HEXANONE_nh4[1], atol=0.0001)]
    maxHexanone = maximum(CalF.DoubleExponential(hum4plot, vec(fitParamsHex)))

    fitParams2Export_rel = fitParams2Export ./ [maxHexanone, 1, maxHexanone, 1, maxHexanone]
    fitParamErrors2Export_rel = fitParamErrors2Export ./ [maxHexanone, 1, maxHexanone, 1, maxHexanone]

    ExpF.exportFitParameters("$(std_filename)fitParameters_relative.txt", fitParams2Export_rel, fitParamErrors2Export_rel,
        measResult.MasslistMasses, measResult.MasslistCompositions;
        fitfunction="relative_sensitivity_to_Hexanone(AH) = p1 * exp.(-p2*AH) .+ p3*exp.(-p4*AH) .+ p5")
end

end # function analyze_humidity_dependence


end # module humidity_dependence

if abspath(PROGRAM_FILE) == @__FILE__
    analyze_humidity_dependence()
end