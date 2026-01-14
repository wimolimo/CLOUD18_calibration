using HDF5
#import PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx

using PyCall
using PyPlot
using Dates
using CSV
using DataFrames
import LsqFit
using TOFTracer2
import TOFTracer2.InterpolationFunctions as IntpF
import TOFTracer2.CalibrationFunctions as CalF
import TOFTracer2.ExportFunctions as ExpF
import TOFTracer2.ImportFunctions as ImpF

fp = joinpath(@__DIR__, "..", "..", "CLOUD18_data", "Calibration", "Humidity-dependent_std", "results")
file = joinpath(fp, "_result.hdf5")

humfile = joinpath(@__DIR__, "..", "..", "CLOUD18_data", "Licor", "2025-11-21.txt")

bg_fp = joinpath(@__DIR__, "..", "..", "CLOUD18_data", "Calibration", "humidity_dependent_BG", "results")
bg_file = joinpath(bg_fp, "_result.hdf5")

plotStart = DateTime(2000, 1, 1, 0, 0, 0)
plotEnd = DateTime(3000, 1, 1, 0, 0, 0)

println("Measurement time range: ", plotStart, " — ", plotEnd)

ions2plot = "NH4+" # "NH4+" # "all", "NH4+", "H+"
#STD_masses_dict = massLibrary.CLOUD_greenSTD_masses # STD1
STD_masses_dict = massLibrary.CLOUD_brownSTD_masses # 
# ["Acetic Acid", "Hexanone", "Acetaldehyde", "Apinene", "Acetonitrile", "Benzene", "Octanone",
# "Xylene", "Hexenal", "MVK", "Toluene", "DMS", "Acetone"]

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
    ion = "H+"
elseif ions2plot == "all"
    for key in ["TMB"] # you choose, which
        append!(massesToPlot, STD_masses_dict[key][1])
    end
    ion = "H+"
end

# massesToPlot = massLibrary.FullPrimaryionslist_NH4soft

##################################
# plot raw data and select filters
##################################

(tracesFig, tracesAx, measResult) = PlotFunctions.plotTracesFromHDF5(file, massesToPlot;
    plotHighTimeRes = false,
    smoothing = 1,
    timeFrame2plot = (plotStart, plotEnd)
    )


# show bg data as well
(tracesFig_bg, tracesAx_bg, measResult_bg) = PlotFunctions.plotTracesFromHDF5(bg_file, massesToPlot;
    plotHighTimeRes = false,
    smoothing = 1,
    timeFrame2plot = (plotStart, plotEnd)

    )

# plot licor data into same figure
humDat = PlotFunctions.load_plotLicorData(humfile; ax=tracesAx, header=2)
tracesFig.tight_layout()

humDat_bg = PlotFunctions.load_plotLicorData(humfile; ax=tracesAx_bg, header=2)
tracesFig_bg.tight_layout()

# bgTimes=[plotStart, plotEnd]
# signalTimes=[plotStart, plotEnd]
# humidityLimits=[0.0, 20] # mmol mol⁻¹ 

# not doing it interactively for now
if !isdefined(Main, :bgTimes)
    println("\ndo the next part interactively:")
    IFIG = PlotFunctions.InteractivePlot(file, tracesAx)
    (humidityLimits, bgTimes, signalTimes) = CalF.humCal_getDatalimitsFromPlot(IFIG)
end



#######################################
# calculate and plot calibration points
#######################################

println("bgTimes: ", bgTimes[1], " — ", bgTimes[2])
println("signalTimes: ", signalTimes[1], " — ", signalTimes[2])
println("humidityLimits: ", humidityLimits)

(calibData, calibData_std, humidities) = CalF.humcal_getHumidityDependentSensitivity(measResult, humDat;
    hums=humidityLimits,
    bgtimes=bgTimes,
    signaltimes=signalTimes,
    pptInInlet=1000)

println("calibData: ", calibData)
println("calibData_std: ", calibData_std)
println("humidities: ", humidities)

fig = figure(figsize=(10, 6))
(calibFig, calibAx) = PlotFunctions.scatter_errorbar(fig, measResult, humidities, calibData, calibData_std; ion=ion)
xlabel("absolute humidity [mmol mol⁻¹]")
ylabel("sensitivity [dcps ppt⁻¹]")

################################
# plot and export fit parameters
################################

hum4plot = collect(0:0.2:12)
fitParams = []
fitParamErrors = []
colornames = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan", "tab:blue", "tab:orange", "tab:green", "tab:red"]

#println(keysToPlot)
println("n_x = ", length(humidities))
println("n_y = ", size(calibData))
println("size of keystoplot: ", size(keysToPlot))
println("size of Masslist: ", size(measResult.MasslistMasses))

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
        println("Fitting mass: ", m)
        (param, stderror, fitlabel) = CalF.fitParameters_DoubleExponential(humidities, calibData[:, i])
        push!(fitParams, param)
        push!(fitParamErrors, stderror)
        plot(hum4plot, CalF.DoubleExponential(hum4plot, param),
            color=colornames[i],
            label="m/z $(round(m,digits=3)), $(MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i])) -- sens(AH) = $(round(param[1],sigdigits=3)) * exp(-$(round(param[2],sigdigits=3))*AH) + $(round(param[3],sigdigits=3))*exp(-$(round(param[4],sigdigits=3))*AH) + $(round(param[5],sigdigits=3))")
    end
end

legend()
savefig("$(fp)calibration_lin_$(ions2plot).png")
yscale("log")
savefig("$(fp)calibration_log_$(ions2plot).png")

fitParams2Export = hvcat(length(measResult.MasslistMasses), (fitParams[a][j] for a in 1:length(measResult.MasslistMasses), j in 1:length(fitParams[1]))...)
fitParamErrors2Export = hvcat(length(measResult.MasslistMasses), (fitParamErrors[a][j] for a in 1:length(measResult.MasslistMasses), j in 1:length(fitParamErrors[1]))...)

ExpF.exportFitParameters("$(fp)fitParameters.txt", fitParams2Export, fitParamErrors2Export,
    measResult.MasslistMasses, measResult.MasslistCompositions;
    fitfunction="sensitivity(AH) = p1 * exp.(-p2*AH) .+ p3*exp.(-p4*AH) .+ p5") #fit function as value?

##########################################################################
# correct fit parameters relative to Hexanone and save these also to file
##########################################################################

#=

if round(massLibrary.HEXANONE_nh4[1],digits=3) in round.(measResult.MasslistMasses,digits=3)
    fitParamsHex = fitParams2Export[:, isapprox.(measResult.MasslistMasses, massLibrary.HEXANONE_nh4[1], atol=0.0001)]
    maxHexanone = maximum(CalF.DoubleExponential(hum4plot, vec(fitParamsHex)))

    fitParams2Export_rel = fitParams2Export ./ [maxHexanone, 1, maxHexanone, 1, maxHexanone]
    fitParamErrors2Export_rel = fitParamErrors2Export ./ [maxHexanone, 1, maxHexanone, 1, maxHexanone]

    ExpF.exportFitParameters("$(fp)fitParameters_relative.txt", fitParams2Export_rel, fitParamErrors2Export_rel,
        measResult.MasslistMasses, measResult.MasslistCompositions;
        fitfunction="relative_sensitivity_to_Hexanone(AH) = p1 * exp.(-p2*AH) .+ p3*exp.(-p4*AH) .+ p5")
end

=#