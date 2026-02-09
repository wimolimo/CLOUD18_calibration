module HumidityDependenceCalibration

export run_humidity_dependence_calibration, get_ion_metadata, get_calibData, print_relative_error_summary, DoubleExponential_and_fit, plot_relative_normalization, export_sensitivities

using HDF5, PyCall, PyPlot, Dates, CSV, DataFrames, Statistics
import LsqFit
using TOFTracer2
import TOFTracer2.MasslistFunctions # FIX: Add missing import
import TOFTracer2.InterpolationFunctions as IntpF
import TOFTracer2.CalibrationFunctions as CalF
import TOFTracer2.ExportFunctions as ExpF

const COLORS = ["blue", "orange", "green", "red", "purple", "brown", "pink", "gray", "olive", "cyan"]

# --- Helper Functions ---

"""
    get_ion_metadata(ions_type::String, std_dict::AbstractDict)

Extracts target m/z values and compound labels from a standard dictionary based on the ionization mode.

# Arguments
- `ions_type`: The ionization mode, typically "NH4+" or "H+".
- `std_dict`: A dictionary (e.g., from `massLibrary`) containing compound names as keys and 
  m/z values as the first element of the value tuple.

# Returns
- A tuple containing: `(Vector{Float64} masses, Vector{String} labels, String ion)`.
"""
function get_ion_metadata(ions_type::String, std_dict::AbstractDict)
    masses, keys_list = Float64[], String[]
    idx = (ions_type == "NH4+") ? 2 : 1
    for key in keys(std_dict)
        push!(masses, std_dict[key][1][idx])
        push!(keys_list, key)
    end
    return masses, keys_list, ions_type
end

"""
    get_calibData(mRes, humdat, mRes_bg; ppt=1000.0, signaltimes=[])

Averages humidity and signal data for each measurement point and calculates sensitivities.

This function interactively prompts the user for an averaging window (in minutes) to be 
applied after each measurement timestamp. It performs duty cycle correction on the TOF 
signals and subtracts background levels if `mRes_bg` is provided.

# Arguments
- `mRes`: Measurement result object for the calibration standards.
- `humdat`: DataFrame containing humidity data from a Licor sensor.
- `mRes_bg`: (Optional) Measurement result object for the background periods.
- `ppt`: The mixing ratio concentration of the standards in ppt (default 1000.0).
- `signaltimes`: A vector of `[start, end]` DateTime to filter the measurement periods.

# Returns
- A tuple: `(calibData, calibData_std, hums_avg, hums_std, window_minutes)`.
"""
function get_calibData(
    mRes::TOFTracer2.ResultFileFunctions.MeasurementResult, 
    humdat::DataFrame, 
    mRes_bg::Union{TOFTracer2.ResultFileFunctions.MeasurementResult, Vector{Any}}=[]; 
    ppt::Float64=1000.0, 
    signaltimes::Vector{DateTime}=[DateTime(0), DateTime(3000)]
)

    print("Enter number of minutes to average humidity after each measurement (press Enter for 4): ")

    input_str = readline()
    window_val = isempty(strip(input_str)) ? 4.0 : parse(Float64, input_str)

    mask = (signaltimes[1] .< mRes.Times) .& (mRes.Times .< signaltimes[2])
    sel_times = mRes.Times[mask]
    traces = mRes.Traces[mask, :] .* transpose(sqrt.(100 ./ mRes.MasslistMasses))   # duty cycle correction

    h_avg, h_std = Float64[], Float64[]
    for t in sel_times
        t_end = t + Second(round(Int, window_val * 60)) # window in seconds
        h_vals = humdat[(humdat.DateTime .>= t) .& (humdat.DateTime .<= t_end), "H₂O_(mmol_mol⁻¹)"]
        push!(h_avg, isempty(h_vals) ? NaN : mean(h_vals))
        push!(h_std, isempty(h_vals) ? NaN : std(h_vals))
    end

    if mRes_bg == []
        calibData, calibData_std = traces ./ ppt, zeros(size(traces))                           # assumption: std compounds all have 1 ppb
    else
        bg_traces = mRes_bg.Traces .* transpose(sqrt.(100 ./ mRes_bg.MasslistMasses))           # duty cycle correction
        calibData = (traces .- mean(bg_traces, dims=1)) ./ ppt                                  # assumption: std compounds all have 1 ppb
        calibData_std = (ones(size(traces, 1)) * std(bg_traces, dims=1)) ./ ppt
    end

    valid = .!(vec(all(isnan.(calibData), dims=2))) .& .!isnan.(h_avg)
    return calibData[valid, :], calibData_std[valid, :], h_avg[valid], h_std[valid], window_val
end

"""
    print_relative_error_summary(h_avg, h_std, calibData, calibData_std, mRes)

Calculates and prints the maximum relative errors (standard deviation / mean) for both 
humidity (X-axis) and sensitivity (Y-axis) for each ion to the terminal.
"""
function print_relative_error_summary(
    h_avg::Vector{Float64}, 
    h_std::Vector{Float64}, 
    calibData::Matrix{Float64}, 
    calibData_std::Matrix{Float64}, 
    mRes::TOFTracer2.ResultFileFunctions.MeasurementResult
)
    rel_err_x = h_std ./ h_avg
    println("\n--- Maximum Relative Error Summary ---")
    println("X (Humidity): ", round(maximum(filter(!isnan, rel_err_x)) * 100, digits=2), "%")
    for i in 1:size(calibData, 2)
        rel_err_y = calibData_std[:, i] ./ calibData[:, i]
        avg_rel_y = maximum(filter(!isnan, rel_err_y)) * 100
        println("Y (m/z $(round(mRes.MasslistMasses[i], digits=3))): ", round(avg_rel_y, digits=3), "%")
    end
    println("--------------------------------------\n")
end

# ---  Plotting & Fitting ---

"""
    DoubleExponential_and_fit(hums, hums_stds, calibData, calibData_std, mRes, ion, out_fp; yscale="log")

Fits a double exponential curve to sensitivity data and generates a visualization.

The fit uses the function: `sensitivity(AH) = p1*exp(-p2*AH) + p3*exp(-p4*AH) + p5`.

# Arguments
- `hums`: Averaged humidity values (mmol/mol or corresponding units).
- `calibData`: Calculated sensitivities (dcps/ppt).
- `mRes`: Metadata object containing mass information for labeling.
- `out_fp`: Directory path where the resulting plot will be saved.
- `yscale`: Scaling for the plot y-axis (e.g., "linear", "log").

# Returns
- A tuple: `(params_matrix, errors_matrix, hum_plot_range)`.
"""
function DoubleExponential_and_fit(
    hums::Vector{Float64}, 
    hums_stds::Vector{Float64}, 
    calibData::Matrix{Float64}, 
    calibData_std::Matrix{Float64}, 
    mRes::TOFTracer2.ResultFileFunctions.MeasurementResult, 
    ion::String, 
    out_fp::String; 
    yscale::String="log"
)

    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan"]

    hums_plot_range = collect(0.0:0.1:maximum(hums)+0.5)
    params_all, params_err_all = [], []
    
    fig, ax = subplots(figsize=(10, 6))
    for i in 1:length(mRes.MasslistMasses)
        p, err, _ = CalF.fitParameters_DoubleExponential(hums, calibData[:, i])
        push!(params_all, p); push!(params_err_all, err)
        
        # Plot data points with error bars and fit curves in matching colors
        ax.errorbar(hums, calibData[:, i], xerr=hums_stds, yerr=calibData_std[:, i], marker="o", linestyle="None", color=colors[mod1(i, 10)], capsize=3)
        ax.plot(hums_plot_range, CalF.DoubleExponential(hums_plot_range, p), color=colors[mod1(i, 10)], 
        label="m/z $(round(mRes.MasslistMasses[i],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray(mRes.MasslistCompositions[:,i], ion = ion))")
    end
    
    ax.set_xlabel("Absolute Humidity [mmol mol⁻¹]")
    ax.set_ylabel("Sensitivity [dcps/ppt]")
    ax.set_yscale(yscale)
    ax.legend(); ax.grid(true, alpha=0.3)
    
    savefig(joinpath(out_fp, "calibration_$(yscale)_$ion.png"))
    
    # Construct matrices for export (5 rows: p1..p5, N columns: masses)
    params_mat = convert(Matrix{Float64}, hvcat(length(params_all), (params_all[a][j] for a in 1:length(params_all), j in 1:5)...))
    params_err_mat = convert(Matrix{Float64}, hvcat(length(params_err_all), (params_err_all[a][j] for a in 1:length(params_err_all), j in 1:5)...))
    return params_mat, params_err_mat, hums_plot_range
end

"""
    export_sensitivities(params_mat, params_err_mat, mRes, filename; out_fp=getcwd())

Exports the calculated fit parameters and their associated uncertainties to a formatted text file.
"""
function export_sensitivities(
    params_mat::Matrix{Float64}, 
    params_err_mat::Matrix{Float64}, 
    mRes::TOFTracer2.ResultFileFunctions.MeasurementResult, 
    filename::String; 
    out_fp::String = getcwd()
)
    full_path = joinpath(out_fp, filename)
    ExpF.exportFitParameters(full_path, params_mat, params_err_mat, mRes.MasslistMasses, mRes.MasslistCompositions;
        fitfunction="sensitivity(AH) = p1*exp(-p2*AH) + p3*exp(-p4*AH) + p5")
    println("Parameters exported to: $full_path")
end

"""
    plot_relative_normalization(h_plot, hums, hums_stds, calibData, calibData_std, params_mat, params_err_mat, mRes, ion, out_fp; yscale="linear")

Normalizes all calibration curves to the maximum sensitivity of Hexanone and plots the result.

This function performs error propagation for the relative sensitivities and generates a 
plot showing how the sensitivity of each ion changes relative to the reference standard 
(Hexanone) over the humidity range.
"""
function plot_relative_normalization(
    h_plot::Vector{Float64}, 
    hums::Vector{Float64}, 
    hums_stds::Vector{Float64}, 
    calibData::Matrix{Float64}, 
    calibData_std::Matrix{Float64}, 
    params_mat::Matrix{Float64}, 
    params_err_mat::Matrix{Float64}, 
    mRes::TOFTracer2.ResultFileFunctions.MeasurementResult, 
    ion::String, 
    out_fp::String; 
    yscale::String="log"
)

    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan"]

    hex_idx = findfirst(isapprox.(mRes.MasslistMasses, massLibrary.HEXANONE_nh4[1], atol=0.0001))
    isnothing(hex_idx) && (println("Hexanone not found, skipping relative plot."); return)      # hexanone must be in the mass list
    
    max_hex = maximum(CalF.DoubleExponential(h_plot, vec(params_mat[:, hex_idx])))
    max_hex_rel_err = params_err_mat[5, hex_idx] / params_mat[5, hex_idx] # plateau error proxy

    p_rel = params_mat ./ [max_hex, 1, max_hex, 1, max_hex]
    e_rel = params_err_mat ./ [max_hex, 1, max_hex, 1, max_hex]

    fig, ax = subplots(figsize=(10, 6))
    for i in 1:length(mRes.MasslistMasses)

        y_rel = calibData[:, i] ./ max_hex
        # σ_rel = y_rel * sqrt((σy/y)^2 + (σhex/hex)^2)
        y_rel_err = y_rel .* sqrt.((calibData_std[:, i] ./ calibData[:, i]).^2 )
        
        ax.plot(h_plot, CalF.DoubleExponential(h_plot, p_rel[:, i]), color=colors[mod1(i, 10)])
        ax.errorbar(hums, y_rel, yerr=y_rel_err, xerr=hums_stds, marker="o", linestyle="None", color=colors[mod1(i, 10)], capsize=3, 
        label="m/z $(round(mRes.MasslistMasses[i],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray(mRes.MasslistCompositions[:,i], ion = "NH3H+"))")
    end
    ax.set_ylabel("Rel. Sensitivity to Hexanone"); ax.set_xlabel("Absolute Humidity"); ax.set_title("Relative Humidity-dependent Calibration")
    ax.set_yscale(yscale); ax.legend(); ax.grid(true, alpha=0.3)
    savefig(joinpath(out_fp, "calibration_relHexanone_$(yscale)_$ion.png"))
    
    export_sensitivities(p_rel, e_rel, mRes, "fitParameters_relative.txt", out_fp=out_fp)
end

# --- main function --------------------------------------------------

"""
    run_humidity_dependence_calibration(fp, calibFile_std, calibFile_bg, humFile; plotTime=...)

Main function for the humidity-dependent calibration workflow.

Loads the HDF5 results and Licor files, processes the averages, performs the 
Double-Exponential fitting, generates all plots, and exports the final parameters.

# Arguments
- `fp`: The root directory for exporting results.
- `calibFile_std`: Path to the standard calibration HDF5 result file.
- `calibFile_bg`: Path to the background HDF5 result file.
- `humFile`: Path to the Licor humidity text file.
- `plotTime`: Vector containing start and end `DateTime` for analysis.
"""
function run_humidity_dependence_calibration(fp::String, calibFile_std::String, calibFile_bg::String, humFile::String; plotTime::Vector{DateTime} = [DateTime(2000, 1, 1), DateTime(3000, 1, 1)])

    plottime_start = plotTime[1]
    plottime_end = plotTime[2]

    # 1. Metadata Selection
    masses_plot, keys_list, ion = get_ion_metadata("NH4+", massLibrary.CLOUD_brownSTD_masses)

    # 2. Raw Loading & Plotting
    (fig, axTraces, legstrings, measResult) = PlotFunctions.plotTracesFromHDF5(calibFile_std, masses_plot; plotHighTimeRes=false, timeFrame2plot=(plottime_start, plottime_end), ion=ion)
    (fig_bg, axTraces_bg, legstrings_bg, measResult_bg) = PlotFunctions.plotTracesFromHDF5(calibFile_bg, masses_plot; plotHighTimeRes=false, timeFrame2plot=(plottime_start, plottime_end), ion=ion)
    
    # Plot BG averages line
    bg_m = vec(mean(measResult_bg.Traces, dims=1))
    for i in 1:length(bg_m) 
        axTraces_bg.axhline(bg_m[i], color="C$(i-1)", linestyle="--", alpha=0.6) 
    end

    # plot Licordata into raw plot
    humdat = PlotFunctions.load_plotLicorData(humFile; ax=axTraces, header=2)
    humdat_bg = PlotFunctions.load_plotLicorData(humFile; ax=axTraces_bg, header=2)

    # 3. Time Averaging Logic
    calibData, calibData_std, hums_avg, hums_stds, win = get_calibData(measResult, humdat, measResult_bg; signaltimes=plotTime)

    # 4. Shade Averaging Windows on data Plot
    ax_hum = axTraces.figure.get_axes()[2]  # twin axis is the last axis added to the figure
    for (i, t) in enumerate(measResult.Times[(plottime_start .< measResult.Times) .& (measResult.Times .< plottime_end)])
        axTraces.axvspan(t, t + Second(round(Int, win * 60)), color="gray", alpha=0.15)
        ax_hum.axhline(hums_avg[i], color = "black", linestyle="--", alpha=0.5)
    end
    
    # 5. Result Summaries, Fitting & Export
    print_relative_error_summary(hums_avg, hums_stds, calibData, calibData_std, measResult)
    
    # Absolute Sensitivities
    params, params_err, fig_sensitivity = DoubleExponential_and_fit(hums_avg, hums_stds, calibData, calibData_std, measResult, ion, fp, yscale="log")
    println("Fitted parameters type: ", typeof(params), typeof(params_err))
    export_sensitivities(params, params_err, measResult, "fitParameters.txt", out_fp=fp)
    
    # Relative Sensitivities
    plot_relative_normalization(fig_sensitivity, hums_avg, hums_stds, calibData, calibData_std, params, params_err, measResult, ion, fp, yscale="log")
    println("Humidity dependence calibration completed.")
end

end # module HumidityDependenceCalibration

#run_humidity_dependence_calibration()