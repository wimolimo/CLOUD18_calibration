#module HumidityDependenceCalibration

using HDF5
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
function get_ion(ions2plot, STD_masses_dict)
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
    return massesToPlot, keysToPlot, ion
end


function plot_bg_avg!(tracesAx, measResult_bg)
    # Calculate and plot horizontal lines for background averages
    bg_means = vec(mean(measResult_bg.Traces, dims=1))
    for i in 1:length(bg_means)
        # Use "C$(i-1)" to match the default color cycle of the traces
        tracesAx.axhline(bg_means[i], color="C$(i-1)", linestyle="--", alpha=0.6)
    end
    tracesFig.canvas.draw()
    tracesFig.tight_layout()

end



#######################################
# calculate and plot calibration points
#######################################

"""
    plot_humidity_dependent_calibration(humcalibfile, ionization)
    
Plot humidity-dependent calibration results and return calibration DataFrame.

# Arguments
- `humcalibfile::String`: Path to humidity calibration file (txt), relative to hexanone.
- `ionization::String`: Ionization method used (e.g., "NH4+").

# Returns
- `calibDF::DataFrame`: DataFrame containing calibration results.

# Saves
- Humidity-dependent relative sensitivity to hexanone plot in the directory of `humcalibfile`.
"""
function plot_humidity_dependent_calibration(humcalibfile, ionization)
    calibDF = CSV.read(humcalibfile, DataFrame; delim='\t', header=2)

    #instead of CalF.plot_humdep_fromCalibParameters:
    hum4plot=collect(0:0.2:18)
    ionization=ionization

    fig, ax = subplots(figsize=(10,6))
    for (name, mass) in zip(calibDF[!, "Sumformula"], calibDF[!, "Mass"])
        f = findfirst(calibDF[!, "Sumformula"] .== name)
        if isnothing(f)
            println("Warning: Could not find formula $name in calibration DataFrame.")
            continue
        end
        # get all params:
        params = calibDF[f, [:p1, :p2, :p3, :p4, :p5]] #fitparameters for this compound
        humdep = CalF.applyFunction(hum4plot,params;functiontype=humdepcalibRelationship="double exponential")
        ax.plot(hum4plot, humdep, label=string(round(mass, digits=3), " - ", name))
    end
    ax.set_xlabel("absolute humidity [mmol mol⁻¹]")
    ax.set_ylabel("relative sensitivity to Hexanone []")
    ax.set_title("Humidity-dependent calibration - Ionization: $(ionization)")
    ax.legend()
    fig.savefig("$(dirname(humcalibfile))/calibration_relHexanone_lin_$(ionization).png")
    ax.set_yscale("log")
    fig.savefig("$(dirname(humcalibfile))/calibration_relHexanone_log_$(ionization).png")

    return calibDF
end

"""
    scatter_errorbar_xy(fig,measResult::ResultFileFunctions.MeasurementResult,xdata::Vector,ydata::Matrix,xerr::Matrix,yerr::Matrix;ion="NH4+")

plots traces as averaged datapoints with their repective given errors
"""
function scatter_errorbar_xy(fig,measResult::ResultFileFunctions.MeasurementResult,xdata::Vector,ydata::Matrix,xerr::Vector,yerr::Matrix;ion="NH4+")
     ax = subplot(111)
     legStrings = []
     if ion in ["all","H+","H3O+"]
         for i = 1:length(measResult.MasslistMasses)
            errorbar(xdata, ydata[:,i], yerr=yerr[:,i], xerr=xerr, marker="o", linestyle="None", capsize=3)
             push!(legStrings,"m/z $(round(measResult.MasslistMasses[i],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i])).H+")
         end
     elseif ion=="NH4+"
         for i = 1:length(measResult.MasslistMasses)
            errorbar(xdata, ydata[:,i], yerr=yerr[:,i], xerr=xerr, marker="o", linestyle="None", capsize=3)
             push!(legStrings,"m/z $(round(measResult.MasslistMasses[i],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray((measResult.MasslistCompositions .- [0,0,3,0,1,0,0,0])[:,i])).NH4+")
         end
     end
    legend(legStrings)
    return fig,ax
end


"""

"""
function humcal_getHumidityDependentSensitivity(mRes, humdat; mRes_bg=[], signaltimes=[DateTime(0),DateTime(3000)], pptInInlet=1.0)
    
    # Ask for user input inside the function
    print("Enter number of minutes to average humidity after each measurement (default 4): ")
    input_str = readline()
    avg_minutes = isempty(strip(input_str)) ? 4.0 : parse(Float64, input_str)

    # Filter mRes rows by signaltimes
    mask_signal = (signaltimes[1] .< mRes.Times) .& (mRes.Times .< signaltimes[2])
    sel_times = mRes.Times[mask_signal]
    Traces_dcps = mRes.Traces[mask_signal, :] .* transpose(sqrt.(100 ./ mRes.MasslistMasses))   # duty cycle correction

    # For each measurement time, compute mean/std of humidity in the next X minutes
    hums_avg = Float64[]
    hums_std = Float64[]
    
    for t_start in sel_times
        # Convert fractional minutes to seconds for duration calculation
        t_end = t_start + Second(round(Int, avg_minutes * 60))
        
        # Filter LiCOR data for this time window
        mask_h2o = (humdat.DateTime .>= t_start) .& (humdat.DateTime .<= t_end)
        h2o_vals = humdat[mask_h2o, "H₂O_(mmol_mol⁻¹)"]
        
        if isempty(h2o_vals)
            push!(hums_avg, NaN)
            push!(hums_std, NaN)
        else
            push!(hums_avg, mean(h2o_vals))
            push!(hums_std, std(h2o_vals))
        end
    end

    # Background subtraction (if mRes_bg provided)
    if mRes_bg == []
        println("No background subtraction applied.")
        calibData = Traces_dcps ./ pptInInlet
        calibData_std = zeros(size(calibData))
    else
        println("Background subtraction applied.")
        Traces_bg = mRes_bg.Traces .* transpose(sqrt.(100 ./ mRes_bg.MasslistMasses))   # duty cycle correction

        # simplify by taking global mean of background for now
        bg_avg = mean(Traces_bg, dims=1)
        bg_std = std(Traces_bg, dims=1)
        
        calibData = (Traces_dcps .- bg_avg) ./ pptInInlet
        # standard deviation propagation (simplified)
        calibData_std = (ones(size(calibData, 1)) * bg_std) ./ pptInInlet
    end

    # Clean up results (remove indices where calibData or humidity is NaN)
    valid_mask = .!(vec(all(isnan.(calibData), dims=2))) .& .!isnan.(hums_avg)
    
    calibData_final = calibData[valid_mask, :]
    calibData_std_final = calibData_std[valid_mask, :]
    hums_avg_final = hums_avg[valid_mask]
    hums_std_final = hums_std[valid_mask]

    return (calibData_final, calibData_std_final, hums_avg_final, hums_std_final, avg_minutes)
 end

 # Removed the readline logic from here and updated the function call





 function relative_err_plot(humidities, hum_stds, calibData, calibData_std, measResult, ions)

     # Calculate relative errors (standard deviation / mean)
    rel_err_x = hum_stds ./ humidities
    rel_err_y = calibData_std ./ calibData

    # Print relative error summary to console
    println("\n--- Error Analysis ---")
    println("Maximum Relative Error X (Humidity): ", round(maximum(filter(!isnan, rel_err_x)) * 100, digits=4), "%")
    for i in 1:size(rel_err_y, 2)
        max_rel_y = maximum(filter(!isnan, rel_err_y[:, i])) * 100
        println("Maximum Relative Error Y (Mass $(round(measResult.MasslistMasses[i], digits=2))): ", round(max_rel_y, digits=4), "%")
    end
    println("----------------------\n")

    fig = figure(figsize=(10, 6))
    (calibFig, calibAx) = scatter_errorbar_xy(fig, measResult, humidities, calibData, hum_stds, calibData_std; ion=ions)
    xlabel("absolute humidity [mmol mol⁻¹]")
    ylabel("sensitivity [dcps ppt⁻¹]")

end

################################
# plot and export fit parameters
################################

function fit_humidity_dependence(humidities, calibData, measResult, ions2plot, keysToPlot, fp)

    hum4plot = collect(0.0:0.1:maximum(humidities)+0.5)
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

end

##########################################################################
# correct fit parameters relative to Hexanone and save these also to file
##########################################################################


function relative_to_hexanone(fitParams2Export, fitParamErrors2Export, measResult)

    # find index of Hexanone
    hex_idx = findfirst(isapprox.(measResult.MasslistMasses, massLibrary.HEXANONE_nh4[1], atol=0.0001))
    fitParamsHex = fitParams2Export[:, hex_idx]
    
    # maxHexanone is usually p1+p3+p5 at AH=0
    maxHexanone = maximum(CalF.DoubleExponential(hum4plot, vec(fitParamsHex)))
    
    fitParams2Export_rel = fitParams2Export ./ [maxHexanone, 1, maxHexanone, 1, maxHexanone]
    fitParamErrors2Export_rel = fitParamErrors2Export ./ [maxHexanone, 1, maxHexanone, 1, maxHexanone]

    println("Exporting fit parameters relative to Hexanone to file.")
    ExpF.exportFitParameters("$(fp)fitParameters_relative.txt", fitParams2Export_rel, fitParamErrors2Export_rel,
        measResult.MasslistMasses, measResult.MasslistCompositions;
        fitfunction="relative_sensitivity_to_Hexanone(AH) = p1 * exp.(-p2*AH) .+ p3*exp.(-p4*AH) .+ p5")

    # Plot sensitivities relative to Hexanone
    fig_rel, ax_rel = subplots(figsize=(10, 6))
    for (i, m) in enumerate(measResult.MasslistMasses)
        # Calculate relative curve using normalized parameters
        rel_params = fitParams2Export_rel[:, i]
        rel_curve = CalF.DoubleExponential(hum4plot, rel_params)
        
        # Calculate Propagated Error for data points
        # σ_rel = y_rel * sqrt((σ_y/y)^2
        y_rel = calibData[:, i] ./ maxHexanone
        y_rel_err = y_rel .* sqrt.((calibData_std[:, i] ./ calibData[:, i]).^2)
        
        label_str = length(keysToPlot) == length(measResult.MasslistMasses) ? 
                    "$(round(m, digits=3)) - $(keysToPlot[i])" : 
                    "m/z $(round(m, digits=3))"
        
        ax_rel.plot(hum4plot, rel_curve, color=colornames[i], label=label_str)
        # Plot markers with propagated error bars in both directions
        ax_rel.errorbar(humidities, y_rel, yerr=y_rel_err, xerr=hum_stds, 
                        marker="o", linestyle="None", color=colornames[i], capsize=3)
    end
    
    ax_rel.set_xlabel("absolute humidity [mmol mol⁻¹]")
    ax_rel.set_ylabel("relative sensitivity to Hexanone []")
    ax_rel.set_title("Humidity-dependent Calibration (Relative to Hexanone) - Ionization: $(ions2plot)")
    ax_rel.legend()
    
    fig_rel.savefig("$(fp)calibration_relHexanone_lin_$(ions2plot).png")
    ax_rel.set_yscale("log")
    fig_rel.savefig("$(fp)calibration_relHexanone_log_$(ions2plot).png")

end

function run_humidity_dependence_calibration()

    massesToPlot, keysToPlot, ion = get_ion(ions2plot, STD_masses_dict)

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

    plot_bg_avg!(tracesAx_bg, measResult_bg)

    # plot licor data into same figure
    humDat = PlotFunctions.load_plotLicorData(humfile; ax=tracesAx, header=2)
    humDat_bg = PlotFunctions.load_plotLicorData(humfile; ax=tracesAx_bg, header=2)

    calibData, calibData_std, humidities, hum_stds, used_avg_mins = humcal_getHumidityDependentSensitivity(measResult, humDat;
    pptInInlet=1000,
    mRes_bg=measResult_bg)

    relative_err_plot(humidities, hum_stds, calibData, calibData_std, measResult, ion)

    # Shade the averaging window in the first plot (tracesAx)
    for t_start in measResult.Times[ (plotStart .< measResult.Times) .& (measResult.Times .< plotEnd) ]
        t_end = t_start + Second(round(Int, used_avg_mins * 60))
        tracesAx.axvspan(t_start, t_end, color="gray", alpha=0.2, label=t_start == measResult.Times[ (plotStart .< measResult.Times) .& (measResult.Times .< plotEnd) ][1] ? "avg window" : "")
    end
    tracesFig.canvas.draw()

    fit_humidity_dependence(humidities, calibData, measResult, ion, keysToPlot, fp)
    
    #get fit parameters to export
    fitParams2Export = hvcat(length(measResult.MasslistMasses), (fitParams[a][j] for a in 1:length(measResult.MasslistMasses), j in 1:length(fitParams[1]))...)
    fitParamErrors2Export = hvcat(length(measResult.MasslistMasses), (fitParamErrors[a][j] for a in 1:length(measResult.MasslistMasses), j in 1:length(fitParamErrors[1]))...)

    ExpF.exportFitParameters("$(fp)fitParameters.txt", fitParams2Export, fitParamErrors2Export,
        measResult.MasslistMasses, measResult.MasslistCompositions;
        fitfunction="sensitivity(AH) = p1 * exp.(-p2*AH) .+ p3*exp.(-p4*AH) .+ p5") 

    if round(massLibrary.HEXANONE_nh4[1],digits=3) in round.(measResult.MasslistMasses,digits=3)
        relative_to_hexanone(fitParams2Export, fitParamErrors2Export, measResult)
    end

end

run_humidity_dependence_calibration()

#end # module HumidityDependenceCalibration
