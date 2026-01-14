module CalibrateTraces

export calibrate_traces_main, CalibrationConfig

#using ... from ...
using HDF5
using CSV
using DataFrames
using Dates
using Plots
using DelimitedFiles
#import Statistics

#check !
using TOFTracer2
import TOFTracer2.CalibrationFunctions as CalF
import TOFTracer2.ImportFunctions as ImpF
import TOFTracer2.InterpolationFunctions as IntpF
import TOFTracer2.ExportFunctions as ExpF
import TOFTracer2.ResultFileFunctions
import TOFTracer2.PlotFunctions
import TOFTracer2.MasslistFunctions
import TOFTracer2.massLibrary

### STRUCTS ###

"""
    CalibrationConfig

Struct to hold configuration parameters for trace calibration.

# Fields
- `dir_licor_data::String`: Directory path containing Licor data files.
- `humcalibfile::String`: Path to humidity calibration file (CSV).
- `drycalibsfile::String`: Path to dry calibration file (HDF5 or CSV).
- `resultfp::String`: File path to save the results.
- `resultfiles::Vector{String}`: List of paths to result files (HDF5).
- `ionization::String`: Ionization method used (e.g., "NH4+").
- `primaryionslist::Vector{Float64}`: List of primary ion masses to load.
- `refMass::Float64`: Mass of the reference compound.
- `refName::String`: Name of the reference compound.
- `exportTraces::Bool`: Flag to indicate whether to export calibrated traces.
- `HeaderForExportDict::Dict{String,Any}`: Dictionary containing header information for export.
"""
struct CalibrationConfig
    dir_licor_data::String
    humcalibfile::String
    drycalibsfile::String
    resultfp::String
    resultfiles::Vector{String}
    ionization::String
    primaryionslist::Vector{Float64}
    refMass::Float64
    refName::String
    exportTraces::Bool
    HeaderForExportDict::Dict{String,Any}
end


### FUNCTIONS ###
"""
    load_hexVSpis_params(drycalibsfile::String)

Load or create and save hexanone vs. primary ion calibration parameters, using TOFTracer2 with CalibrationFunctions.dryCal_selectPIandRefDataFromIFIG.

# Arguments
- `drycalibsfile::String`: Path to dry calibration file (HDF5 or CSV)

# Returns
- `hexVSpis_params`: Tuple containing parameters, errors, and metadata

# Saves
- CSV file with hexanone vs. primary ion calibration parameters if loaded from HDF5, in the same directory as `drycalibsfile`.
"""
function load_hexVSpis_params(drycalibsfile::String)
    if HDF5.ishdf5(drycalibsfile)
        hexVSpis_params = CalF.dryCal_selectPIandRefDataFromIFIG(drycalibsfile)
        hexVSpis_params2export = vcat(hexVSpis_params[3], ["parameters" "errors"], hcat(hexVSpis_params[1], hexVSpis_params[2]))
        CSV.write("$(dirname(drycalibsfile))Hexanone_VS_PIs_params.csv", DataFrame(hexVSpis_params2export, :auto)) #save parameters for later use?
    else
        a = CSV.read(drycalibsfile, DataFrame, header=[2])
        b = CSV.read(drycalibsfile, DataFrame; footerskip=3, header=false)
        hexVSpis_params = (a.parameters, a.errors, [values(b[1, :])[1], values(b[1, :])[2]])
    end
    return hexVSpis_params
end


"""
    load_licor_data(dir_licor_data::String)
    
Load Licor humidity data from text files in the specified directory.

# Arguments
- `dir_licor_data::String`: Directory path containing Licor data files.

# Returns
- `licorDat::DataFrame`: DataFrame containing loaded Licor humidity data (datetime, H2O_mmolpermol).
"""
function load_licor_data(dir_licor_data::String)
    return ImpF.createLicorData_fromFiles(dir_licor_data;
        filefilter=r"licor_.*\.txt", #rename file to use with licor_restofthename.txt
        headerrow=2,
        columnNameOfInterest="H₂O_(mmol_mol⁻¹)",
        type_columnOfInterest=Float64)
end


"""
    plot_humidity_dependent_calibration(humcalibfile, ionization)
    
Plot humidity-dependent calibration results and return calibration DataFrame.

# Arguments
- `humcalibfile::String`: Path to humidity calibration file (txt).
- `ionization::String`: Ionization method used (e.g., "NH4+").

# Returns
- `calibDF::DataFrame`: DataFrame containing calibration results.

# Saves
- Humidity-dependent relative sensitivity to hexanone plot in the directory of `humcalibfile`.
"""
function plot_humidity_dependent_calibration(humcalibfile, ionization)
    calibDF = CSV.read(humcalibfile, DataFrame; header=2)

    #instead of CalF.plot_humdep_fromCalibParameters:
    hum4plot=collect(0:0.2:12)
    ionization=ionization

    plt = plot()
    for (name, mass) in zip(calibDF[!, "Sumformula"], calibDF[!, "Mass"])
        f = findfirst(calibDF[!, "Sumformula"] .== name)
        # get all params:
        params = calibDF[f, [:p1, :p2, :p3, :p4, :p5]]
        humdep = applyFunction(hum4plot,params;functiontype=humdepcalibRelationship="double exponential")
        plot(hum4plot, humdep, label=string(round(mass, digits=3), " - ", name))
    end
    xlabel!(plt, "absolute humidity [mmol mol⁻¹]")
    ylabel!(plt, "relative sensitivity to Hexanone []")
    # linear y-scale
    plot!(plt, legend = :best)
    savefig(plt, "$(dirname(humcalibfile))calibration_relHexanone_lin_$(ionization).png")
    # logarithmic y-scale
    plot!(plt, yscale = :log10)
    savefig(plt, "$(dirname(humcalibfile))calibration_relHexanone_log_$(ionization).png")

    return calibDF
end
#should it also return the plot? or just save it?


"""
    load_and_merge_results(resultfiles, primaryionslist)

Load and merge measurement results from multiple result files and return the merged results for all masses and for the primary ions only.

# Arguments
- `resultfiles::Vector{String}`: List of paths to result files (HDF5).
- `primaryionslist::Vector{Float64}`: List of primary ion masses to load.

# Returns
- `mResfinal::struct`: Merged measurement results for all masses (selectedTimes, selectedMasslistMasses, masslistElements, masslistElementsMasses, selectedMassesCompositions, traces).
- `mResfinal_PIs::struct`: Merged measurement results for primary ions only.
"""
function load_and_merge_results(resultfiles, primaryionslist)
    mResfinal_PIs = ResultFileFunctions.loadResults(resultfiles[1]; useAveragesOnly = true, massesToLoad = primaryionslist)
    mResfinal = ResultFileFunctions.loadResults(resultfiles[1]; useAveragesOnly = true)

    if length(resultfiles) > 1
        for i in 2:length(resultfiles)
            mResfinal_PIs = ResultFileFunctions.joinResultsTime(mResfinal_PIs, ResultFileFunctions.loadResults(resultfiles[i]; useAveragesOnly = true, massesToLoad = primaryionslist))
            mResfinal = ResultFileFunctions.joinResultsTime(mResfinal, ResultFileFunctions.loadResults(resultfiles[i]; useAveragesOnly = true))
        end
    end

    return mResfinal, mResfinal_PIs
end


"""
    compute_summed_primary_ions(mResfinal_PIs)

Returns the summed primary ions from the measurement results.

# Arguments
- `mResfinal_PIs::struct`: Measurement results for primary ion intensities.

# Returns
- `summedPIs::Vector`: Summed primary ion intensities.
"""
function compute_summed_primary_ions(mResfinal_PIs)
    summedPIs = mResfinal_PIs.Traces .* sqrt.(100 ./ mResfinal_PIs.MasslistMasses) #Traces is matrix times vs PI
    summedPIs[summedPIs .<= 0] .= 0
    return summedPIs
end


"""
    interpolate_licor_to_ptr_time(mResfinal, licorDat)

Return interpolated Licor humidity data to match the time points of the measurement results.

# Arguments
- `mResfinal::struct`: Measurement results containing time points.
- `licorDat::DataFrame`: DataFrame containing Licor humidity data (datetime, H2O_mmolpermol).

# Returns
- `licor_final::Vector`: Interpolated Licor humidity data aligned with measurement results time points. 
"""
function interpolate_licor_to_ptr_time(mResfinal, licorDat)
    return IntpF.sortSelectAverageSmoothInterpolate(
        mResfinal.Times,
        licorDat.datetime,
        licorDat.H2O_mmolpermol,
        returnSTdev = false,
        selectY = [-Inf, Inf] # option to get rid of outliers via the kwarg 'selectY' (use wisely) #[-21, 5]
    )
end



"""
    build_calibration_traces(mResfinal, summedPIs, licor_final, calibDF, hexVSpis_params, refName)

Build calibration traces for all compounds based on humidity-dependent and dry calibration.

# Arguments
- `mResfinal::struct`: Measurement results containing mass list and time points.
- `summedPIs::Vector`: Summed primary ion intensities.
- `licor_final::Vector`: Interpolated Licor humidity data aligned with measurement results time points.
- `calibDF::DataFrame`: Calibration data frame.
- `hexVSpis_params::Vector`: Hexanone vs primary ion parameters.
- `refName::String`: Reference compound name.

# Returns
- `dcps_per_ppb::Matrix`: Calibration traces in dcps per ppb for all compounds.
- `indices::Vector{Int}`: Indices of calibrated compounds.
"""
function build_calibration_traces(mResfinal, summedPIs, licor_final, calibDF, hexVSpis_params, refName)

    #initialize calibration traces matrix with zeros
    dcps_per_ppb = zeros(length(mResfinal.Times), length(mResfinal.MasslistMasses))

    fref = findfirst(calibDF[!, "Sumformula"] .== refName)
    refparams = calibDF[fref, [:p1, :p2, :p3, :p4, :p5]]

    #composition based filters
    undeffilter = BitVector(sum(mResfinal.MasslistCompositions; dims=1)[1, :] .== 0)
    oxygen_number = findfirst(==("O"), mResfinal.MasslistElements)
    oneoxygenfilter = BitVector(mResfinal.MasslistCompositions[oxygen_number, :] .== 1)
    twoplusoxygenfilter = BitVector(mResfinal.MasslistCompositions[oxygen_number, :] .>= 2)

    println("found $(sum(undeffilter)) undefined masses, $(sum(oneoxygenfilter)) masses with 1 oxygen atom, and $(sum(twoplusoxygenfilter)) masses with >=2 oxygen atoms")
    
    #Hexanone dry sensitivity [cps/ppb] vs primary ion cps, calculated for all time points
    f_hex = CalF.applyFunction(summedPIs, hexVSpis_params[1]; functiontype = hexVSpis_params[3][1])

    #multiply f_hex with wet sensitivity of different masses, relative to dry Hexanone
    #dry / kinetic limit calibration for undef and >=2 O
    println("calibrating all compounds with >2 oxygen atoms and undefined ones with reference $(refName) dry.")
    dcps_per_ppb[:, (undeffilter .| twoplusoxygenfilter)] .= #vector of length Times 
        f_hex .*
        CalF.applyFunction(zeros(length(mResfinal.Times)), refparams; functiontype = "double exponential") #vector containing the zero humidity point of the humidity dependent calibration

    #humid / equilibrium calibration for 1 O
    println("calibrating all compounds with 1 oxygen atom humidity-dependent with reference $(refName).")
    dcps_per_ppb[:, oneoxygenfilter] .=
        f_hex .*
        CalF.applyFunction(licor_final, refparams; functiontype = "double exponential") #use licor_final instead of CalF.applyFunction(fpfinal, humparams[1]; functiontype = humparams[3][1]) only if icor data is complete

    indices = Int[]

 #= for (name, mass) in zip(calibDF[!, "Sumformula"], calibDF[!, "Mass"])
        params = calibDF[findfirst(calibDF[!, "Sumformula"] .== name), [:p1, :p2, :p3, :p4, :p5]]
        index = findfirst(isapprox.(mResfinal.MasslistMasses, mass, atol=0.0001))
 =#
    for row in eachrow(calibDF)
        params = row[[:p1, :p2, :p3, :p4, :p5]]
        index = findfirst(isapprox.(mResfinal.MasslistMasses, row.Mass; atol=0.0001))

        if index isa Int
            dcps_per_ppb[:, index] =
                f_hex .*
                CalF.applyFunction(licor_final, params; functiontype = "double exponential") 

            push!(indices, index)
        end
    end

    return dcps_per_ppb, indices
end
#what I changed: remove frostpoint dependency, use only licor instead. use humidity dependant calibration for 1 oxygen and kinetic limit calibration for 2+ oxygen.
#compounds with 0 oxygen are ignored (zero sensitivity)


"""
    plot_calibration_traces(mResfinal, dcps_per_ppb, summedPIs, indices, resultfp)

Plot calibration traces for selected compounds and save the plots. 

# Arguments
- `mResfinal::struct`: Measurement results containing mass list and time points.
- `dcps_per_ppb::Matrix`: Calibration traces in dcps per ppb for all compounds.
- `summedPIs::Vector`: Summed primary ion intensities.
- `indices::Vector{Int}`: Indices of calibrated compounds.
- `resultfp::String`: File path to save the plots.

# Saves
- Calibration traces plot as PNG and PDF files.
"""
function plot_calibration_traces(mResfinal, dcps_per_ppb, summedPIs, indices, resultfp)
    p1 = plot(size=(1000,600),
        xlabel="time [UTC]",
        ylabel="calibration factor [dcps / ppb]")

    for i in indices
        plot!(p1,
            mResfinal.Times,
            dcps_per_ppb[:, i],
            label = "index $i - $(round(mResfinal.MasslistMasses[i],digits=2)), " * MasslistFunctions.sumFormulaStringFromCompositionArray(mResfinal.MasslistCompositions[:,i])
        )
    end

    plot!(p1, mResfinal.Times, summedPIs, label="summed primary ions", lw=2)

    savefig(p1, "$(resultfp)CalibrationTraces.png")
    savefig(p1, "$(resultfp)CalibrationTraces.pdf")
end


"""
    export_calibrated_traces(mResfinal, dcps_per_ppb, indices, ionization, HeaderForExportDict, resultfp)

Export calibrated traces for masses with exactly one nitrogen and more than one oxygen to a CSV file with specified headers.

# Arguments
- `mResfinal::struct`: Measurement results containing mass list and time points.
- `dcps_per_ppb::Matrix`: Calibration traces in dcps per ppb for all compounds.
- `ionization::String`: Ionization method used (e.g., "NH4+").
- `HeaderForExportDict::Dict{String,Any}`: Dictionary containing header information for export.
- `resultfp::String`: File path to save the exported CSV file.

# Saves
- Exported calibrated traces as a CSV file.
"""
function export_calibrated_traces(mResfinal, dcps_per_ppb, ionization, HeaderForExportDict, resultfp)

    #filter for masses with at least one carbon
    filterCnr = mResfinal.MasslistCompositions[findfirst(mResfinal.MasslistElements .== "C"), :] .>= 1

    #filter for masses with calibration data
    filterNoCalib = vec(sum(dcps_per_ppb; dims=1) .> 0)

    #filter for masses with exactly one nitrogen
    filterNnr = mResfinal.MasslistCompositions[findfirst(mResfinal.MasslistElements .== "N"), :] .== 1 #adjust if want to export other masses

    finalfilter = filterCnr .& filterNoCalib .& filterNnr

    calibResult =
        ResultFileFunctions.MeasurementResult(
            mResfinal.Times,
            mResfinal.MasslistMasses[finalfilter],
            mResfinal.MasslistElements,
            mResfinal.MasslistElementsMasses,
            mResfinal.MasslistCompositions[:,finalfilter],
            1000.0 .* mResfinal.Traces[:,finalfilter] ./ dcps_per_ppb[:,finalfilter]
        )

    iifig = PlotFunctions.InteractivePlot(calibResult)
    PlotFunctions.scrollAddTraces(iifig)
    IndOfinterest = unique(iifig.activeIndices)

    HeaderForExport = TOFTracer2.ExportFunctions.CLOUDheader(calibResult.Times;
            title = HeaderForExportDict["title"],
            level=HeaderForExportDict["level"],
            version=HeaderForExportDict["version"],
            authorname_mail=HeaderForExportDict["authorname_mail"],
            units=HeaderForExportDict["units"],
            addcomment=HeaderForExportDict["addcomment"],
            threshold=HeaderForExportDict["threshold"],
            nrrows_addcomment = HeaderForExportDict["nrrows_addcomment"]
    )

    ExpF.exportTracesCSV_CLOUD(
        resultfp,
        calibResult.MasslistElements,
        calibResult.MasslistMasses[IndOfinterest],
        calibResult.MasslistCompositions[:,IndOfinterest],
        calibResult.Times,
        calibResult.Traces[:,IndOfinterest];
        transmission = 0,
        headers = HeaderForExport,
        ion = ionization,
        average = 0
    )
end


################################
# Main Script
################################

"""
    calibrate_traces_main(config::CalibrationConfig)

Main function to calibrate measurement traces based on humidity-dependent and dry calibration.

# Arguments
- `config::CalibrationConfig`: Struct containing paths and parameters for calibration.

# Saves
- Humidity-dependent calibration plot in the directory of `humcalibfile`.
- Calibration traces plot as PNG and PDF files.
- Exported calibrated traces as CSV file.
"""
function calibrate_traces_main(config::CalibrationConfig)

    hexVSpis_params =
#        !isdefined(hexVSpis_params) ?
            load_hexVSpis_params(config.drycalibsfile)
#            :hexVSpis_params

    licorDat =
#        !isdefined(licorDat) ?
            load_licor_data(config.dir_licor_data)
#            :licorDat


    calibDF = plot_humidity_dependent_calibration(config.humcalibfile, config.ionization)

    primaryionslist =
        isempty(config.primaryionslist) ?
            massLibrary.FullPrimaryionslist_NH4soft :
            config.primaryionslist

    mResfinal, mResfinal_PIs = load_and_merge_results(config.resultfiles, primaryionslist)
    
    summedPIs = compute_summed_primary_ions(mResfinal_PIs)

    licor_final = interpolate_licor_to_ptr_time(mResfinal, licorDat)

    dcps_per_ppb, indices = build_calibration_traces(mResfinal, summedPIs, licor_final, calibDF, hexVSpis_params, config.refName)

    plot_calibration_traces(mResfinal, dcps_per_ppb, summedPIs, indices, config.resultfp)

    if config.exportTraces
        export_calibrated_traces(mResfinal, dcps_per_ppb, config.ionization, config.HeaderForExportDict, config.resultfp)
    end
end

#error estimation still missing

end # module CalibrateTraces


