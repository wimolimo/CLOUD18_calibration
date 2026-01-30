#module HumidityDependenceCalibration

#export run_humidity_dependence_calibration, get_ion_metadata, compute_window_averages, print_relative_error_summary, fit_and_export_sensitivities, plot_relative_normalization

using HDF5, PyCall, PyPlot, Dates, CSV, DataFrames, Statistics
import LsqFit
using TOFTracer2
import TOFTracer2.InterpolationFunctions as IntpF
import TOFTracer2.CalibrationFunctions as CalF
import TOFTracer2.ExportFunctions as ExpF
import TOFTracer2.ImportFunctions as ImpF

# --- Global Configurations ---
const FP = joinpath(@__DIR__, "..", "..", "CLOUD18_data", "Calibration", "Humidity-dependent_std", "results")
const FILE = joinpath(FP, "_result.hdf5")
const HUMFILE = joinpath(@__DIR__, "..", "..", "CLOUD18_data", "Licor", "2025-11-21.txt")
const BG_FILE = joinpath(@__DIR__, "..", "..", "CLOUD18_data", "Calibration", "humidity_dependent_BG", "results", "_result.hdf5")

const PLOT_START = DateTime(2000, 1, 1)
const PLOT_END = DateTime(3000, 1, 1)
const COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan"]




# --- Helper Functions ---

"""
    get_ion_metadata(ions_type, std_dict)

This function extracts the target m/z values and compound labels from a provided standard 
dictionary, filtering by the selected ionization mode (e.g., "NH4+" or "H+").
"""
function get_ion_metadata(ions_type, std_dict)
    masses, keys_list = Float64[], String[]
    idx = (ions_type == "NH4+") ? 2 : 1
    for key in keys(std_dict)
        push!(masses, std_dict[key][1][idx])
        push!(keys_list, key)
    end
    return masses, keys_list, ions_type
end

"""
    compute_window_averages(mRes, humdat, mRes_bg; ppt=1000.0, signaltimes=[])

This function handles the time-averaging logic. It interactively prompts the user for 
an averaging window, calculates humidity and signal statistics (mean/std) for each 
measurement point, and performs background subtraction if a background dataset is provided.
"""
function compute_window_averages(mRes, humdat, mRes_bg; ppt=1000.0, signaltimes=[DateTime(0), DateTime(3000)])
    print("Enter number of minutes to average humidity after each measurement (default 4): ")
    input_str = readline()
    window_val = isempty(strip(input_str)) ? 4.0 : parse(Float64, input_str)

    mask = (signaltimes[1] .< mRes.Times) .& (mRes.Times .< signaltimes[2])
    sel_times = mRes.Times[mask]
    traces = mRes.Traces[mask, :] .* transpose(sqrt.(100 ./ mRes.MasslistMasses))

    h_avg, h_std = Float64[], Float64[]
    for t in sel_times
        t_end = t + Second(round(Int, window_val * 60)) # window in seconds
        h_vals = humdat[(humdat.DateTime .>= t) .& (humdat.DateTime .<= t_end), "H₂O_(mmol_mol⁻¹)"]
        push!(h_avg, isempty(h_vals) ? NaN : mean(h_vals))
        push!(h_std, isempty(h_vals) ? NaN : std(h_vals))
    end

    if mRes_bg == []
        cal_data, cal_std = traces ./ ppt, zeros(size(traces))
    else
        bg_traces = mRes_bg.Traces .* transpose(sqrt.(100 ./ mRes_bg.MasslistMasses))
        cal_data = (traces .- mean(bg_traces, dims=1)) ./ ppt
        cal_std = (ones(size(traces, 1)) * std(bg_traces, dims=1)) ./ ppt
    end

    valid = .!(vec(all(isnan.(cal_data), dims=2))) .& .!isnan.(h_avg)
    return cal_data[valid, :], cal_std[valid, :], h_avg[valid], h_std[valid], window_val
end

"""
    print_relative_error_summary(h_avg, h_std, c_data, c_std, mRes)

This function calculates and prints a summary of the average relative errors (std/mean) 
for the humidity (X-axis) and the sensitivities of each ion (Y-axis) to the console.
"""
function print_relative_error_summary(h_avg, h_std, c_data, c_std, mRes)
    rel_err_x = h_std ./ h_avg
    println("\n--- Maximum Relative Error Summary ---")
    println("X (Humidity): ", round(maximum(filter(!isnan, rel_err_x)) * 100, digits=2), "%")
    for i in 1:size(c_data, 2)
        rel_err_y = c_std[:, i] ./ c_data[:, i]
        avg_rel_y = maximum(filter(!isnan, rel_err_y)) * 100
        println("Y (m/z $(round(mRes.MasslistMasses[i], digits=3))): ", round(avg_rel_y, digits=2), "%")
    end
    println("--------------------------------------\n")
end

# --- main Plotting & Fitting ---

"""
    fit_and_plot_sensitivities(hums, hums_stds, c_data, c_std, mRes, keys_list, ion, out_fp; yscale="linear")

This function performs a double-exponential fit for each ion in the mass list. It generates 
a visualization containing the averaged data points with error bars and the resulting 
fit curves. The final plot is saved to the specified output directory.

Returns `(fitParamsMatrix, fitErrorsMatrix, humAxisForPlot)`.
"""
function fit_and_plot_sensitivities(hums, hums_stds, c_data, c_std, mRes, keys_list, ion, out_fp; yscale="linear")
    h_plot = collect(0.0:0.1:maximum(hums)+0.5)
    p_all, e_all = [], []
    
    fig, ax = subplots(figsize=(10, 6))
    for i in 1:length(mRes.MasslistMasses)
        p, err, _ = CalF.fitParameters_DoubleExponential(hums, c_data[:, i])
        push!(p_all, p); push!(e_all, err)
        
        lbl = i <= length(keys_list) ? "$(round(mRes.MasslistMasses[i], digits=3)) - $(keys_list[i])" : "m/z $(round(mRes.MasslistMasses[i], digits=3))"
        # Plot data points with error bars and fit curves in matching colors
        ax.errorbar(hums, c_data[:, i], xerr=hums_stds, yerr=c_std[:, i], marker="o", linestyle="None", color=COLORS[mod1(i, 10)], capsize=3)
        ax.plot(h_plot, CalF.DoubleExponential(h_plot, p), color=COLORS[mod1(i, 10)], label=lbl)
    end
    
    ax.set_xlabel("Absolute Humidity [mmol mol⁻¹]")
    ax.set_ylabel("Sensitivity [dcps/ppt]")
    ax.set_yscale(yscale)
    ax.legend(); ax.grid(true, alpha=0.3)
    
    savefig(joinpath(out_fp, "calibration_$(yscale)_$ion.png"))
    
    # Construct matrices for export (5 rows: p1..p5, N columns: masses)
    p_mat = hvcat(length(p_all), (p_all[a][j] for a in 1:length(p_all), j in 1:5)...)
    e_mat = hvcat(length(e_all), (e_all[a][j] for a in 1:length(e_all), j in 1:5)...)
    
    return p_mat, e_mat, h_plot
end

"""
    export_sensitivities(p_mat, e_mat, mRes, out_fp, filename)

This function handles the file I/O operations for the calibration parameters. It formats 
the fitting coefficients and their associated uncertainties into a structured text file 
using the TOFTracer2 export utility.
"""
function export_sensitivities(p_mat::Matrix{Float64}, e_mat::Matrix{Float64}, mRes, filename; out_fp = getcwd())

    full_path = joinpath(out_fp, filename)
    ExpF.exportFitParameters(full_path, p_mat, e_mat, mRes.MasslistMasses, mRes.MasslistCompositions;
        fitfunction="sensitivity(AH) = p1*exp(-p2*AH) + p3*exp(-p4*AH) + p5")
    println("Parameters exported to: $full_path")
end

"""
    plot_relative_normalization(h_plot, hums, hums_stds, c_data, c_std, p_mat, e_mat, mRes, keys_list, ion, out_fp; yscale="linear")

This function normalizes the sensitivities to the maximum value of Hexanone. It 
performs error propagation for the relative data, generates normalized plots with 
scatter points and error bars, and exports the relative fit parameters.
"""
function plot_relative_normalization(h_plot, hums, hums_stds, calibData, calibData_std, p_mat, e_mat, mRes, keys_list, ion, out_fp; yscale="linear")

    hex_idx = findfirst(isapprox.(mRes.MasslistMasses, massLibrary.HEXANONE_nh4[1], atol=0.0001))
    isnothing(hex_idx) && (println("Hexanone not found, skipping relative plot."); return)      # hexanone must be in the mass list
    
    max_hex = maximum(CalF.DoubleExponential(h_plot, vec(p_mat[:, hex_idx])))
    max_hex_rel_err = e_mat[5, hex_idx] / p_mat[5, hex_idx] # plateau error proxy

    p_rel = p_mat ./ [max_hex, 1, max_hex, 1, max_hex]
    e_rel = e_mat ./ [max_hex, 1, max_hex, 1, max_hex]

    fig, ax = subplots(figsize=(10, 6))
    for i in 1:length(mRes.MasslistMasses)

        y_rel = calibData[:, i] ./ max_hex
        # σ_rel = y_rel * sqrt((σy/y)^2 + (σhex/hex)^2)
        y_rel_err = y_rel .* sqrt.((calibData_std[:, i] ./ calibData[:, i]).^2 )
        
        ax.plot(h_plot, CalF.DoubleExponential(h_plot, p_rel[:, i]), color=COLORS[mod1(i, 10)])
        ax.errorbar(hums, y_rel, yerr=y_rel_err, xerr=hums_stds, marker="o", linestyle="None", color=COLORS[mod1(i, 10)], capsize=3, label=keys_list[i])
    end
    ax.set_ylabel("Rel. Sensitivity to Hexanone"); ax.set_xlabel("Absolute Humidity"); ax.set_title("Relative Humidity-dependent Calibration")
    ax.set_yscale(yscale); ax.legend(); ax.grid(true, alpha=0.3)
    savefig(joinpath(out_fp, "calibration_relHexanone_$(yscale)_$ion.png"))
    
    export_sensitivities(p_rel, e_rel, mRes, out_fp, "fitParameters_relative.txt")
end

# --- Orchestrator ---

function run_humidity_dependence_calibration()

    # 1. Metadata Selection
    m_plot, keys_list, ion = get_ion_metadata("NH4+", massLibrary.CLOUD_brownSTD_masses)

    println("Output get_ion_metadata: ", m_plot, keys_list, ion)

    # 2. Raw Loading & Plotting
    (fig_raw, ax_raw, measResult) = PlotFunctions.plotTracesFromHDF5(FILE, m_plot; plotHighTimeRes=false, timeFrame2plot=(PLOT_START, PLOT_END))
    (_, ax_bg, measResult_bg) = PlotFunctions.plotTracesFromHDF5(BG_FILE, m_plot; plotHighTimeRes=false, timeFrame2plot=(PLOT_START, PLOT_END))
    
    # Plot BG averages line
    bg_m = vec(mean(measResult_bg.Traces, dims=1))
    for i in 1:length(bg_m) 
        ax_bg.axhline(bg_m[i], color="C$(i-1)", linestyle="--", alpha=0.6) 
    end

    # plot Licordata into raw plot
    humdat = PlotFunctions.load_plotLicorData(HUMFILE; ax=ax_raw, header=2)
    _ = PlotFunctions.load_plotLicorData(HUMFILE; ax=ax_bg, header=2)

    # 3. Time Averaging Logic
    c_data, c_std, hums_avg, hums_stds, win = compute_window_averages(measResult, humdat, measResult_bg; signaltimes=[PLOT_START, PLOT_END])
    println("c_data, c_std, hums_avg, hums_stds sizes, win: ", size(c_data), size(c_std), length(hums_avg), length(hums_stds), win)

    # 4. Shade Averaging Windows on Raw Plot
    for t in measResult.Times[(PLOT_START .< measResult.Times) .& (measResult.Times .< PLOT_END)]
        ax_raw.axvspan(t, t + Second(round(Int, win * 60)), color="gray", alpha=0.15)
    end
    
    # 5. Result Summaries, Fitting & Export
    print_relative_error_summary(hums_avg, hums_stds, c_data, c_std, measResult)
    
    # Absolute Sensitivities
    p_exp, e_exp, h_plot = fit_and_plot_sensitivities(hums_avg, hums_stds, c_data, c_std, measResult, keys_list, ion, FP)
    export_sensitivities(p_exp, e_exp, measResult, FP, "fitParameters.txt")
    
    # Relative Sensitivities
    plot_relative_normalization(h_plot, hums_avg, hums_stds, c_data, c_std, p_exp, e_exp, measResult, keys_list, ion, FP, yscale="log")
    println("Humidity dependence calibration completed.")
end

#end # module HumidityDependenceCalibration

#using .HumidityDependenceCalibration

run_humidity_dependence_calibration()

