using HDF5
#import PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx

using PyCall
using PyPlot
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

bgTimes=[plotStart, plotEnd]
signalTimes=[plotStart, plotEnd]
humidityLimits=[0.0, 20] # mmol mol⁻¹ 

# not doing it interactively for now
if !isdefined(Main, :bgTimes)
    println("\ndo the next part interactively:")
    IFIG = PlotFunctions.InteractivePlot(file, tracesAx)
    (humidityLimits, bgTimes, signalTimes) = CalF.humCal_getDatalimitsFromPlot(IFIG)
end

##############################################
# get humidity ranges around measurement times
##############################################
struct output
        avg::Float64
        start_time::DateTime
        stop_time::DateTime
        indices::Vector{Vector{Int64}}
        values::Vector{Float64}
    end 

out = Vector{output}()

"""
    getHumidityRange(humDat, query_times; max_abs_dev=0.2, max_rel_dev=0.05, min_points=3)

    Given a DataFrame `humDat` with licor humidity data and a vector of `query_times` (DateTime),
    this function finds, for each query time, a range of humidity values around the closest measurement
    time that form a plateau within specified tolerances. The plateau is defined as the largest contiguous set of humidity
    values around the query time that do not deviate from the center value (at the query time) by more than
    `max_abs_dev` (absolute) or `max_rel_dev` (relative). If the plateau contains fewer than `min_points` points,
    the window is expanded symmetrically until at least `min_points` are included.
Returns a vector of NamedTuples, each containing:
- `avg`: the average humidity in the identified range
- `start_time`: the DateTime of the first point in the range
- `stop_time`: the DateTime of the last point in the range
- `indices`: the indices of the original DataFrame corresponding to the points in the range
"""    
function getHumidityRange(humDat::DataFrame, query_times::Vector{DateTime};
        max_abs_dev = 0.2,        # absolute mmol mol⁻¹ tolerance for plateau
        max_rel_dev = 0.05,       # relative tolerance (fraction of center)
        min_points = 3)           # ensure at least this many points in window

    empty!(out)

    date_col_str = "System_Date_(Y-M-D)"
    time_col_str = "System_Time_(h:m:s)"
    hum_col_str  = "H₂O_(mmol_mol⁻¹)"

    names_str = String.(names(humDat))
    !(date_col_str in names_str && time_col_str in names_str && hum_col_str in names_str) &&
        throw(ArgumentError("Expected licor columns not found. Available: $(names(humDat))"))

    # build DateTime vector from date+time columns (index by string column names)
    n = nrow(humDat)
    licor_times = Vector{DateTime}(undef, n)
    fmt1 = DateFormat("yyyy-mm-dd H:M:S")
    fmt2 = DateFormat("yyyy-mm-dd H:M")
    for i in 1:n
        s = string(humDat[i, date_col_str]) * " " * string(humDat[i, time_col_str])
        licor_times[i] = try
            DateTime(s, fmt1)
        catch
            try
                DateTime(s, fmt2)
            catch
                DateTime(s)   # fallback generic parser
            end
        end
    end

    # numeric humidity (NaN for missing / non-convertible)
    raw = humDat[!, hum_col_str]
    hum = Array{Float64}(undef, length(raw))
    for i in eachindex(raw)
        v = raw[i]
        hum[i] = v === missing ? NaN : try Float64(v) catch NaN end
    end

    length(licor_times) == length(hum) || throw(ArgumentError("licor_times length mismatch"))

    valid_mask = .!map(isnan, hum) .& .!map(ismissing, licor_times)
    if count(valid_mask) == 0
        throw(ArgumentError("No valid humidity/time rows found in humDat"))
    end

    times_v = licor_times[valid_mask]
    hum_v = hum[valid_mask]
    orig_indices = findall(valid_mask)

    for q in query_times

        # get closest humidity measurement to query time
        dvals = abs.(Dates.value.(q .- times_v))
        idx = argmin(dvals)
        center = hum_v[idx]

        within_tol(val) = abs(val - center) <= max_abs_dev || abs(val - center) <= max_rel_dev * max(abs(center), eps())

        l = idx
        while l > 1 && within_tol(hum_v[l-1])
            l -= 1
        end

        r = idx
        while r < length(hum_v) && within_tol(hum_v[r+1])
            r += 1
        end

        if (r - l + 1) < min_points
            extra = min_points - (r - l + 1)
            addl = min(extra ÷ 2 + extra % 2, l - 1)
            addr = min(extra ÷ 2, length(hum_v) - r)
            l -= addl
            r += addr
        end

        values = hum_v[l:r]
        avg = mean(values)
        inds_in_original = orig_indices[l:r]
        # push!(out, (avg = avg, start_time = times_v[l], stop_time = times_v[r], indices = inds_in_original, values = values))
        push!(out, output(avg, times_v[l], times_v[r], [inds_in_original], values))
    end

    return out
end

# get humidity averages for each measurement time of calibration data
# humDat_time = humDat[!,"System_Time_(h:m:s)"]
#timeFilter = PlotFunctions.matplotlib2datetime.(tracesAx.get_xlim()[1]) .< humDat_time .< PlotFunctions.matplotlib2datetime.(tracesAx.get_xlim()[2])
#humtime = humDat_time[timeFilter]
results = getHumidityRange(humDat, measResult.Times)
avgs = [r.avg for r in results]

#starttimes = [r.start_time for r in results]
#endtimes = [r.stop_time for r in results]
println("Averages: ", avgs)


# plot humidity points in both figures (std and bg)
fig1 = tracesAx[:figure]
ax_humidity = filter(ax -> ax[:get_ylabel]() == "Humidity [mmol mol⁻¹]", fig1[:axes])[1]
#comparison values are interpolated humidity at measurement times, not averages
comparison_value = IntpF.interpolateSelect(measResult.Times,humDat.DateTime,humDat[!,"H₂O_(mmol_mol⁻¹)"];selTimes=[DateTime(0),DateTime(3000)])

# use center time of each selected window and plot as black markers
centers = measResult.Times
try
    ax_humidity[:scatter](centers, avgs; s=40, c="k", zorder=5)
    ax_humidity[:scatter](centers, comparison_value; s=40, c="r", zorder=5)
    # dashed horizontal segments for each averaging window
    for r in results
        ax_humidity[:plot](
            [r.start_time, r.stop_time],
            [r.avg, r.avg];
            linestyle="--",
            color="r",
            linewidth=1.2,
            zorder=4
        )
    end
    # optional: label each point with its value (rounded)
    for (c, a) in zip(centers, avgs)
        ax_humidity[:text](c, a-1.0, string(round(a, digits=3)); fontsize=8, va="bottom", ha="center", color="k")
        ax_humidity[:text](c, a+0.5, string(round(comparison_value[findfirst(==(c), measResult.Times)], digits=3)); fontsize=8, va="bottom", ha="center", color="r")
    end
catch err
    @warn "Could not plot averages on tracesAx_bg: $err"
end

# now for bg data
results_bg = getHumidityRange(humDat_bg, measResult_bg.Times)
avgs_bg = [r.avg for r in results_bg]
starttimes_bg = [r.start_time for r in results_bg]
endtimes_bg = [r.stop_time for r in results_bg]

fig2 = tracesAx_bg[:figure]
ax_humidity_bg = filter(ax -> ax[:get_ylabel]() == "Humidity [mmol mol⁻¹]", fig2[:axes])[1]
# use center time of each selected window and plot as black markers
centers_bg = measResult_bg.Times

ax_humidity_bg[:scatter](centers_bg, avgs_bg; s=40, c="k", zorder=5)
# dashed horizontal segments for each averaging window
for r in results_bg
    ax_humidity_bg[:plot](
        [r.start_time, r.stop_time],
        [r.avg, r.avg];
        linestyle="--",
        color="r",
        linewidth=1.2,
        zorder=4
    )
end
# label each point with its value (rounded)
for (c, a) in zip(centers_bg, avgs_bg)
    ax_humidity_bg[:text](c, a, string(round(a, digits=3)); fontsize=8, va="bottom", ha="center", color="k")
end

tracesFig.tight_layout()
tracesFig_bg.tight_layout()

#######################################
# calculate and plot calibration points
#######################################

"""
	humcal_getHumidityDependentSensitivityOfBG(mRes,humdat;hums=collect(0,0.1,1),bgtimes=[],signaltimes=[DateTime(0),DateTime(3000)],pptInInlet=1.0)

calculates from a processed dataset of a humidity-dependent calibration and an output file of a LiCOR the humidity dependent calibration factors
"""
function humcal_getHumidityDependentSensitivityOfBG(mRes,humdat; signaltimes=[DateTime(0),DateTime(3000)],pptInInlet=1.0)
	
    hums = IntpF.interpolateSelect(mRes.Times,humdat.DateTime,humdat[!,"H₂O_(mmol_mol⁻¹)"];selTimes=signaltimes)

	Traces_dcps = mRes.Traces .* transpose(sqrt.(100 ./mRes.MasslistMasses))

    humRanges = getHumidityRange(humdat, mRes.Times[signaltimes[1] .< mRes.Times .< signaltimes[2]])

    println("hums: ", hums)
    println("Traces_dcps size: ", size(Traces_dcps))

    calibData = (Traces_dcps)./(pptInInlet)
    #calibData_std_i = (signalVShum_std)./(pptInInlet)

	# delete for Nans
	calibData_noNaN = calibData[.!(vec(all(isnan.(calibData),dims=2))),:]
	#calibData_std_noNaN = calibData_std[.!(vec(all(isnan.(calibData),dims=2))),:]
	hums_noNaN = hums[.!(vec(all(isnan.(calibData),dims=2))),:]

	return (calibData_noNaN,vec(hums_noNaN))
end


println("bgTimes: ", bgTimes[1], " — ", bgTimes[2])
println("signalTimes: ", signalTimes[1], " — ", signalTimes[2])
println("humidityLimits: ", humidityLimits)

calibData, humidities = humcal_getHumidityDependentSensitivityOfBG(measResult, humDat;
    signaltimes=signalTimes,
    pptInInlet=1000)

calibData_std = zeros(size(calibData))

println("calibData: ", calibData)
println("humidities: ", humidities)

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

    ExpF.exportFitParameters("$(fp)fitParameters_relative.txt", fitParams2Export_rel, fitParamErrors2Export_rel,
        measResult.MasslistMasses, measResult.MasslistCompositions;
        fitfunction="relative_sensitivity_to_Hexanone(AH) = p1 * exp.(-p2*AH) .+ p3*exp.(-p4*AH) .+ p5")
end

