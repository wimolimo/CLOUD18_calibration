using HDF5
using PythonCall
using PythonPlot
using Dates
using CSV
using DataFrames
using Statistics
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

random_file = joinpath("C://Users//c7441399//Documents//Atemluft", "2026-01-22-beginn-der-aufzeichnungen.txt")

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

#######################################
# calculate and plot calibration points
#######################################

"""
    humcal_getHumidityDependentSensitivityOfBG(mRes, humdat; signaltimes=[DateTime(0),DateTime(3000)], pptInInlet=1.0)

    Calculate the humidity-dependent sensitivity of the measurement result `mRes`
    using the humidity data `humdat`. The sensitivity is calculated as the ratio
    of the measured traces (corrected for mass-dependent transmission) to the
    known ppt in the inlet, for the specified signal times.

    # Arguments
    - `mRes::MeasurementResult`: The measurement result containing traces and times.
    - `humdat::DataFrame`: The humidity data containing DateTime and humidity values.
    - `signaltimes::Vector{DateTime}`: A vector with two DateTime values specifying
      the start and end times for signal averaging.
    - `pptInInlet::Float64`: The known ppt concentration in the inlet.

    # Returns
    - `calibData_noNaN::Array{Float64,2}`: The calibration data (sensitivity) without NaN values.
    - `hums_noNaN::Vector{Float64}`: The corresponding humidity values without NaN entries.
"""
function humcal_getHumidityDependentSensitivityOfBG(mRes,humdat; mRes_bg=[],humdat_bg=[],signaltimes=[DateTime(0),DateTime(3000)],pptInInlet=1.0)
	
    # if bg is given, get bg humidity dependency and subtract bg
    if !isempty(mRes_bg) && !isempty(humdat_bg)
        hums_bg = IntpF.interpolateSelect(mRes_bg.Times,humdat_bg.DateTime,humdat_bg[!,"H₂O_(mmol_mol⁻¹)"];selTimes=signaltimes)
        Traces_dcps_bg = mRes_bg.Traces .* transpose(sqrt.(100 ./mRes_bg.MasslistMasses)) # duty cycle correction
        calibData_bg = (Traces_dcps_bg)./(pptInInlet)     # umrechnung in ppb

        # delete for Nans
        calibData_bg_noNaN = calibData_bg[.!(vec(all(isnan.(calibData_bg),dims=2))),:]
        hums_bg_noNaN = hums_bg[.!(vec(all(isnan.(calibData_bg),dims=2))),:]

        for (i, m) in enumerate(measResult.MasslistMasses)
            println("Fitting mass: ", m)
            (param, stderror, fitlabel) = CalF.fitParameters_DoubleExponential(hums_bg_noNaN, calibData_bg_noNaN[:, i])
            push!(fitParams, param)
            push!(fitParamErrors, stderror)
            plot(hum4plot, CalF.DoubleExponential(hum4plot, param),
            color=colornames[i],
            label="m/z $(round(m,digits=3)), $(MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i])) -- sens(AH) = $(round(param[1],sigdigits=3)) * exp(-$(round(param[2],sigdigits=3))*AH) + $(round(param[3],sigdigits=3))*exp(-$(round(param[4],sigdigits=3))*AH) + $(round(param[5],sigdigits=3))")
        end
    end


    hums = IntpF.interpolateSelect(mRes.Times,humdat.DateTime,humdat[!,"H₂O_(mmol_mol⁻¹)"];selTimes=signaltimes)

	Traces_dcps = mRes.Traces .* transpose(sqrt.(100 ./mRes.MasslistMasses)) # duty cycle correction

    calibData = (Traces_dcps)./(pptInInlet)     # umrechnung in ppb
    #calibData_std_i = (signalVShum_std)./(pptInInlet)

	# delete for Nans
	calibData_noNaN = calibData[.!(vec(all(isnan.(calibData),dims=2))),:]
	#calibData_std_noNaN = calibData_std[.!(vec(all(isnan.(calibData),dims=2))),:]
	hums_noNaN = hums[.!(vec(all(isnan.(calibData),dims=2))),:]

	return (calibData_noNaN,vec(hums_noNaN))
end

calibData_bg, humidities_bg = humcal_getHumidityDependentSensitivityOfBG(measResult_bg, humDat_bg;
    pptInInlet=1000)

#TODO: get humidity dependent bg and subtract from measResult, but for now, just bg averaged over whole time
#calibData, humidities = humcal_getHumidityDependentSensitivityOfBG(measResult, humDat;
#    pptInInlet=1000)

# put std to zeros for now
calibData_std = zeros(size(calibData))

println("calibData: ", size(calibData))
println("humidities: ", size(humidities))

#=
(calibData, calibData_std, humidities) = CalF.humcal_getHumidityDependentSensitivity(measResult, humDat;
    hums=avgs,
    signaltimes=signalTimes,
    pptInInlet=1000)

println("calibData: ", calibData)
println("calibData_std: ", calibData_std)
println("humidities: ", humidities)
=#

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
    println("Masses not found: ", setdiff(massesToPlot, measResult.MasslistMasses))
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


if round(massLibrary.HEXANONE_nh4[1],digits=3) in round.(measResult.MasslistMasses,digits=3)
    fitParamsHex = fitParams2Export[:, isapprox.(measResult.MasslistMasses, massLibrary.HEXANONE_nh4[1], atol=0.0001)]
    maxHexanone = maximum(CalF.DoubleExponential(hum4plot, vec(fitParamsHex)))

    fitParams2Export_rel = fitParams2Export ./ [maxHexanone, 1, maxHexanone, 1, maxHexanone]
    fitParamErrors2Export_rel = fitParamErrors2Export ./ [maxHexanone, 1, maxHexanone, 1, maxHexanone]

    println("Exporting fit parameters relative to Hexanone to file.")
    ExpF.exportFitParameters("$(fp)fitParameters_relative.txt", fitParams2Export_rel, fitParamErrors2Export_rel,
        measResult.MasslistMasses, measResult.MasslistCompositions;
        fitfunction="relative_sensitivity_to_Hexanone(AH) = p1 * exp.(-p2*AH) .+ p3*exp.(-p4*AH) .+ p5")

    # TODO: plot relative sensitivities
end

