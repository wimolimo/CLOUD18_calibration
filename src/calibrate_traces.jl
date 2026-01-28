module CalibrateTraces

export calibrate_traces_main, CalibrationConfig

#using ... from ...
using HDF5
using CSV
using DataFrames
using Dates
using DelimitedFiles
using PyPlot
using PyCall
import Statistics

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

### STRUCT ###
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
### helper functions ###
"""
    parse_formula_to_composition(formula::String)

Parse a chemical formula string and return a dictionary of element counts.

# Arguments
- `formula::String`: Chemical formula string (e.g., "C6H15NO.H+" or "C6H15H+NO")

# Returns
- `Dict{String,Int}`: Dictionary mapping element symbols to their counts, ignoring charges.
"""
function parse_formula_to_composition(formula::String)
    
    formula_clean = replace(formula, r"\d*[+-]" => "") # handles multiple charges, parentheses not handled
    composition = Dict{String,Int}()
    
    # Match element symbol (uppercase letter optionally followed by lowercase) -> group [1]
    # followed by optional number -> group [2]
    # ignores everything that is not an element symbol + number (e.g. +, -, ., etc)
    pattern = r"([A-Z][a-z]?)(\d*)" #regrex pattern to match elements and their counts
    
    for match in eachmatch(pattern, formula_clean) # iterate over all formulas in the structure of pattern, ignoring charges
        element = match.captures[1] 
        count_str = match.captures[2]
        count = isempty(count_str) ? 1 : parse(Int, count_str) # default count is 1 if not specified, otherwise parse the number
        composition[element] = get(composition, element, 0) + count # accumulate counts for each element: get element count. if element does not exist asume count=0; add current count
    end
    
    return composition
end


"""
    find_formula_index(formulas, target_formula::String)

Find the index of a formula in a vector by comparing atomic composition.
Returns `nothing` if no match is found.

# Arguments
- `formulas`: Vector of formula strings from DataFrame (any string-like type)
- `target_formula::String`: Target formula to search for

# Returns
- `Int` or `Nothing`: Index of matching formula or nothing if not found
"""
function find_formula_index(formulas, target_formula::String)
    target_composition = parse_formula_to_composition(target_formula)
    
    for (i, formula) in enumerate(formulas)
        formula_composition = parse_formula_to_composition(String(formula))
        if formula_composition == target_composition
            return i #if more than 1 would match, return first
        end
    end
    return nothing # if no match was found
end


"""
    scatterDryCalibs2(drycalibsfile::String; referenceMasses=[TOFTracer2.massLibrary.HEXANONE_nh4[1]],primaryions=[])

Plot dry calibration data for primary ions and reference masses from the specified dry calibration file.

# Arguments
- `drycalibsfile::String`: Path to dry calibration file (HDF5).
- `referenceMasses::Vector{Float64}`: List of reference masses to plot (default: hexanone + NH4+).
- `primaryions::Vector{Float64}`: List of primary ion masses to plot (default: all water and ammonium clusters).

# Returns
- `dryCalibFig`: Figure object of the dry calibration plot.
- `dryCalibAx`: Axes object of the dry calibration plot.
- `mResDryCalibs`: Measurement results loaded from the dry calibration file.
- `primaryiontraces`: Traces of primary ions.
- `referencetraces`: Traces of reference masses.
"""
function scatterDryCalibs2(drycalibsfile::String; referenceMasses=[TOFTracer2.massLibrary.HEXANONE_nh4[1]],primaryions=[]) #modified from PlotFunctions.scatterDryCalibs
    if isempty(primaryions)
        primaryions = [
            MasslistFunctions.massFromComposition(H=2, O=1)
            MasslistFunctions.massFromComposition(H=4, O=2)
            MasslistFunctions.massFromComposition(H=6, O=3)
            MasslistFunctions.massFromComposition(H=8, O=4)
            MasslistFunctions.massFromComposition(H=3, N=1)
            MasslistFunctions.massFromComposition(H=5, N=1, O=1)
            MasslistFunctions.massFromComposition(H=7, N=1, O=2)
            MasslistFunctions.massFromComposition(H=9, N=1, O=3)
            MasslistFunctions.massFromComposition(H=6, N=2)
            MasslistFunctions.massFromComposition(H=8, N=2, O=1)
            MasslistFunctions.massFromComposition(H=10, N=2, O=2)
            MasslistFunctions.massFromComposition(H=9, N=3)
            MasslistFunctions.massFromComposition(H=11, N=3, O=1)
        ]
    end
    massesDryCalibToPlot = vcat(primaryions, referenceMasses)
    
    mResDryCalibs = ResultFileFunctions.loadResults(drycalibsfile; useAveragesOnly=true, massesToLoad=massesDryCalibToPlot)# load data
    dryCalibFig = figure()
    dryCalibAx = subplot(111)
    
    filterarray = falses(length(mResDryCalibs.MasslistMasses)) # create bitarray based on occurence of primaryions in massesDryCalibToPlot
    for m in primaryions
        filterarray .|= isapprox.(mResDryCalibs.MasslistMasses,m;atol=0.00001)
    end

    primaryionmasses = mResDryCalibs.MasslistMasses[filterarray]
    primaryiontraces = mResDryCalibs.Traces[:,filterarray] * sqrt.(100 ./ primaryionmasses)
    referencetraces = mResDryCalibs.Traces[:,(!).(filterarray)] * sqrt.(100 ./ referenceMasses)
    
    scatter(mResDryCalibs.Times, primaryiontraces, label="sum of primary ions") #plot primary ions summed dcps trace
    scatter(mResDryCalibs.Times, referencetraces, label="sum of reference ions - m/z $(round.(referenceMasses;digits=3))") # plot reference dcps trace
    xlabel("Time")
    ylabel("signals [dcps]")
    title("Dry Calibrations")
    legend()
    yscale("log")
    grid()
    tight_layout()
    savefig("$(dirname(drycalibsfile))dryCalibs.png")
    return dryCalibFig, dryCalibAx, mResDryCalibs, primaryiontraces, referencetraces
end


"""
    dryCal_selectPIandRefDataInteractive(drycalibsfile::String)

Interactively select dry calibration data points to exclude from fit and return fitted hexanone vs primary ion parameters.

# Arguments
- `drycalibsfile::String`: Path to dry calibration file (hdf5).

# Returns
- `hexVSpis_params::Tuple`: Tuple containing fit parameters, errors, and functiontype.
"""
function dryCal_selectPIandRefDataInteractive(drycalibsfile::String) #modified from CalibrationFunctions.dryCal_selectPIandRefDataFromIFIG

    dryCalibFig, dryCalibAx, mResDryCalibs, primaryiontraces, referencetraces = scatterDryCalibs2(drycalibsfile; referenceMasses=[TOFTracer2.massLibrary.HEXANONE_nh4[1]], primaryions=[])
    println("please give the minimum y-value to show")
    dryCalibAx.set_ylim(bottom=parse(Int, readline()))
    println("How many dry calibration data points do you want to exclude from the fit?")
    nrOfExcludeCalibs = parse(Int, readline())
    IFIG = PlotFunctions.InteractivePlot(drycalibsfile, dryCalibAx)
    println("Select primary ion calibration coordinates to exclude by moving the mouse to the respective primary ion data point and press 'c'. Repeat until you have selected $nrOfExcludeCalibs point(s).")
    PlotFunctions.getMouseCoords(IFIG; datetime_x=true)
    while (length(IFIG.coords) < nrOfExcludeCalibs)
        sleep(0.1)
    end
    
    #Find closest datapoints
    exclude_idx = Int[]
    for i in 1:nrOfExcludeCalibs
        t_click  = IFIG.coords[i][1] # ISO 8601 datetime string like "2025-10-11T09:31:27.405"
        y_click  = IFIG.coords[i][2]
        t_click_unix = Dates.datetime2unix(DateTime(t_click))
        t_data_unix = Dates.datetime2unix.(mResDryCalibs.Times)
        dists = (t_data_unix .- t_click_unix).^2 .+ (primaryiontraces .- y_click).^2 # Distance^2 in (Time, PrimaryIonsSum) space ######check datatype of Times: 15-element Vector{Dates.DateTime}: 2025-10-02T13:07:34.863, 2025-10-07T10:00:04.705, 2025-10-08T14:19:40.794; use matplotlib2datetime or similar?
        push!(exclude_idx, argmin(dists))
    end
    exclude_idx = unique(exclude_idx) # in case user clicked very close points, double selections are excluded
    df = DataFrame(
        Time = mResDryCalibs.Times[Not(exclude_idx)], #exclude selected time points via their Indices
        PrimaryIonsSum = vec(sum(primaryiontraces[Not(exclude_idx), :]; dims=2)), #sum across masses and convert to vector
        ReferenceSignal = vec(sum(referencetraces[Not(exclude_idx), :]; dims=2)) #sum across masses and convert to vector
    )

    # Fit
    hexVSpis_params = CalF.fitParameters(df.PrimaryIonsSum, df.ReferenceSignal; functiontype="power")
    nrOfCalibs = nrow(df)

    figure()
    scatter(df.PrimaryIonsSum, df.ReferenceSignal, label="data")
    xforfit = collect(floor(minimum(df.PrimaryIonsSum); sigdigits=1):1000: ceil(maximum(df.PrimaryIonsSum); sigdigits=1))
    fill_between(xforfit,
        CalF.PowerFunction(xforfit, hexVSpis_params[1] .- hexVSpis_params[2] / sqrt(nrOfCalibs)),
        CalF.PowerFunction(xforfit, hexVSpis_params[1] .+ hexVSpis_params[2] / sqrt(nrOfCalibs)),
        label="uncertainty", 
        alpha=0.25)
    plot(xforfit, CalF.PowerFunction(xforfit, hexVSpis_params[1]), label=hexVSpis_params[3])
    legend()
    xlabel("sum of primary ions [dcps]")
    ylabel("signal on reference mass [dcps/ppb]")
    yscale("log")
    savefig("$(dirname(drycalibsfile))Hexanone_VS_PIs.png")
    savefig("$(dirname(drycalibsfile))Hexanone_VS_PIs.pdf")
    return hexVSpis_params
end


#the rest is modified from calibration script
### calibration steps functions ###
"""
    load_hexVSpis_params(drycalibsfile::String)

Load or create and save hexanone vs. primary ion calibration parameters, using TOFTracer2 with CalibrationFunctions.dryCal_selectPIandRefDataFromIFIG.

# Arguments
- `drycalibsfile::String`: Path to dry calibration file (HDF5 or CSV)

# Returns
- `hexVSpis_params::Tuple`: Tuple containing parameters, errors, and metadata

# Saves
- CSV file with hexanone vs. primary ion calibration parameters if loaded from HDF5, in the same directory as `drycalibsfile`.
"""
function load_hexVSpis_params(drycalibsfile::String)
    if HDF5.ishdf5(drycalibsfile)
        hexVSpis_params = dryCal_selectPIandRefDataInteractive(drycalibsfile)
        hexVSpis_params2export = vcat(hexVSpis_params[3], ["parameters" "errors"], hcat(hexVSpis_params[1], hexVSpis_params[2]))
        CSV.write("$(dirname(drycalibsfile))Hexanone_VS_PIs_params.csv", DataFrame(hexVSpis_params2export, :auto)) #save parameters for later use?
    else
        a = CSV.read(drycalibsfile, DataFrame, header=[2])
        b = CSV.read(drycalibsfile, DataFrame; footerskip=3, header=false)
        hexVSpis_params = (a.parameters, a.errors, [values(b[1, :])[1], values(b[1, :])[2]])
    end
    return hexVSpis_params
end
# change to load from CSV if exists?

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
    calibDF = CSV.read(humcalibfile, DataFrame; delim='\t', header=2)

    ################ Plot humidity-dependent calibration results: move this to timo???#############################
    #instead of CalF.plot_humdep_fromCalibParameters:
    hum4plot=collect(0:0.2:12)
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

    #mResfinal_PIs.MasslistCompositions = #########TO DO: rewrite name from C10H15H+.NH4+ to C10H16.NH4+

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
#    summedPIs = mResfinal_PIs.Traces .* sqrt.(100 ./ mResfinal_PIs.MasslistMasses) #Traces is matrix times vs PI
    # Traces is matrix (time × masses), MasslistMasses is vector (masses,)
    # Multiply each column by sqrt(100/mass), then sum across masses to get one value per time point
    summedPIs = vec(sum(mResfinal_PIs.Traces .* reshape(sqrt.(100 ./ mResfinal_PIs.MasslistMasses), 1, :); dims=2))
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
        licorDat.H2O_mmolpermol;
        returnSTdev = false,
        selectY = [-Inf, Inf] # option to get rid of outliers via the kwarg 'selectY' (use wisely) #[-21, 5]
    )
end


"""
    build_calibration_traces(mResfinal, summedPIs, licor_final, calibDF, hexVSpis_params, refName)

Build calibration traces for all compounds based on humidity-dependent and dry calibration. User selects wheather to apply humidity-dependent calibration.

# Arguments
- `mResfinal::MeasurementResult`: Measurement results containing mass list and time points.
- `summedPIs::Vector`: Summed primary ion intensities.
- `licor_final::Vector`: Interpolated Licor humidity data aligned with measurement results time points.
- `calibDF::DataFrame`: Humidity dependent calibration data.
- `hexVSpis_params::Tuple`: Hexanone vs primary ion parameters.
- `refName::String`: Reference compound name.

# Returns
- `dcps_per_ppb::Matrix`: Calibration traces in dcps per ppb for all compounds, for each humidity and PI sum.
- `dcps_per_ppb_err::Matrix`: Uncertainties of calibration traces in dcps per ppb.
- `indices::Vector{Int}`: Indices of calibrated compounds.
"""
function build_calibration_traces(mResfinal, summedPIs, licor_final, calibDF, hexVSpis_params, refName)

    #initialize calibration traces matrix with zeros
    dcps_per_ppb = zeros(length(mResfinal.Times), length(mResfinal.MasslistMasses))
    dcps_per_ppb_err = zeros(size(dcps_per_ppb))

    fref = find_formula_index(calibDF[!, "Sumformula"], refName) # Match reference compound by atomic composition instead of exact string matching
    if fref === nothing
        error("Reference compound '$refName' not found in calibration data. Available formulas: $(calibDF[!, "Sumformula"])")
    end
    
    ref_params = calibDF[fref, [:p1, :p2, :p3, :p4, :p5]]

    #composition based filters
    undeffilter = BitVector(sum(mResfinal.MasslistCompositions; dims=1)[1, :] .== 0)
    oxygen_number = findfirst(==("O"), mResfinal.MasslistElements)
    oneoxygenfilter = BitVector(mResfinal.MasslistCompositions[oxygen_number, :] .== 1)
    twoplusoxygenfilter = BitVector(mResfinal.MasslistCompositions[oxygen_number, :] .>= 2)

    println("found $(sum(undeffilter)) undefined masses, $(sum(oneoxygenfilter)) masses with 1 oxygen atom, and $(sum(twoplusoxygenfilter)) masses with >=2 oxygen atoms")
    
    #Hexanone dry sensitivity [cps/ppb] vs primary ion cps, calculated for all time points
    f_hex = CalF.applyFunction(summedPIs, hexVSpis_params[1]; functiontype = hexVSpis_params[3][1]) #summedPI has length times???
    f_hex_err = sqrt.( 
        (hexVSpis_params[2][1]./(hexVSpis_params[1][1]))^2 #relative variance of parameter a (slope), (delta a / a)^2
        .+ (log.(summedPIs[summedPIs .> 0]).*hexVSpis_params[2][2]).^2 #relative variance of parameter b (exponent), (ln x *delta b)^2
        .+ 2*(hexVSpis_params[2][1]./(hexVSpis_params[1][1])*log.(summedPIs[summedPIs .> 0]).*hexVSpis_params[2][2]) #covariance term (2 delta a/a * ln x * delta b); 2 cov(a,b) = delta a * delta b is assumed
    )
    f_hex_err[summedPIs .<= 0] .= 0

    # wet sensitivity of different masses
    f_hum0 = CalF.applyFunction(zeros(length(mResfinal.Times)), ref_params; functiontype = "double exponential") #vector containing the zero humidity point of the humidity dependent calibration
    f_hum_err = zeros(length(mResfinal.Times)) #placeholder for later addition of relative humidity dependent calibration error. Note that this part holds ONLY true, if the errors from the humidity-dependent calibration fit are negligible compared to the errors of the fit to the primary ions

    println("Do you want to apply the humidity-dependent calibration for 1 oxygen compounds (recommended for T>0°C)? (y/n)")
    userinput = readline()
    if userinput == "y"
        f_hum = CalF.applyFunction(licor_final, ref_params; functiontype = "double exponential") # licor_final has length of times #use licor_final instead of CalF.applyFunction(fpfinal, humparams[1]; functiontype = humparams[3][1]) only if icor data is complete
        
        #dry / kinetic limit calibration for undef and >=2 O
        println("calibrating all compounds with >=2 oxygen atoms and undefined ones with reference $(refName) dry.")
        dcps_per_ppb[:, (undeffilter .| twoplusoxygenfilter)] .= f_hex .* f_hum0
        dcps_per_ppb_err[:, (undeffilter .| twoplusoxygenfilter)] .= dcps_per_ppb[:, (undeffilter .| twoplusoxygenfilter)] .* sqrt.(f_hex_err.^2 .+ 0.0.^2) # no humidity dependent calibration error for dry calibration

        #humid / equilibrium calibration for 1 O
        println("calibrating all compounds with 1 oxygen atom humidity-dependent with reference $(refName).")
        dcps_per_ppb[:, oneoxygenfilter] .= f_hex .* f_hum
        dcps_per_ppb_err[:, oneoxygenfilter] .= dcps_per_ppb[:, oneoxygenfilter] .* sqrt.(f_hex_err.^2 .+ f_hum_err.^2) # combine errors from dry and humid calibration

        indices = Int[]
    
        for row in eachrow(calibDF)
            params = row[[:p1, :p2, :p3, :p4, :p5]]
            index = findfirst(isapprox.(mResfinal.MasslistMasses, row.Mass; atol=0.0001)) #find masses in masslist matching compounds in calibDF (humidity dependent calibration data)
            f_hum_row = CalF.applyFunction(licor_final, params; functiontype = "double exponential") 

            if index isa Int
                dcps_per_ppb[:, index] = f_hex .* f_hum_row
                dcps_per_ppb_err[:, index] = dcps_per_ppb[:, index] .* sqrt.(f_hex_err.^2 .+ f_hum_err.^2) # combine errors from dry and humid calibration
                push!(indices, index)
            end
        end
        
    elseif userinput == "n"
        f_hum = f_hum0
        println("calibrating all compounds (except 0 oxygen compounds) dry with reference $(refName).")
        dcps_per_ppb[:, (undeffilter .| twoplusoxygenfilter .| oneoxygenfilter)] .= f_hex .* f_hum0
        dcps_per_ppb_err[:, (undeffilter .| twoplusoxygenfilter .| oneoxygenfilter)] .= dcps_per_ppb[:, (undeffilter .| twoplusoxygenfilter .| oneoxygenfilter)] .* sqrt.(f_hex_err.^2 .+ 0.0.^2) # no humidity dependent calibration error for dry calibration

    else
        error("Invalid input for humidity-dependent calibration choice. Please enter 'y' or 'n'.")
    end


    total_relative_error = sqrt.(f_hex_err .^2 .+ f_hum_err .^2)
    mean_relative_error = Statistics.mean(total_relative_error)
    std_of_mean_relative_error = Statistics.std(total_relative_error)
    println("The relative standard error of the calibration factor for this method alone is a factor ",
        round(mean_relative_error,digits=3), " ± ",round(std_of_mean_relative_error,digits=3),
        " due to the error of the normalization to the primary ions.") #change print when f_hum_err != 0

    return dcps_per_ppb, dcps_per_ppb_err, indices
end
#what I changed: remove frostpoint dependency, use only licor instead. use humidity dependant calibration for 1 oxygen and kinetic limit calibration for 2+ oxygen.
#added error here, but f_hum_err is still zero
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
    fig, ax = subplots(figsize = (10, 6))

    ax.set_xlabel("time [UTC]")
    ax.set_ylabel("calibration factor [dcps / ppb]")
    ax.set_yscale("log")

    for i in indices
        ax.plot(mResfinal.Times, dcps_per_ppb[:, i],
            label = "$(round(mResfinal.MasslistMasses[i],digits=2)) * $(MasslistFunctions.sumFormulaStringFromCompositionArray(mResfinal.MasslistCompositions[:,i]))"
        )
    end

    ax.plot(mResfinal.Times, summedPIs, label="summed primary ions", linewidth=2)
    ax.legend()
    ax.set_title("Calibration Traces")

    fig.savefig(joinpath(resultfp, "CalibrationTraces.png"))
    fig.savefig(joinpath(resultfp, "CalibrationTraces.pdf"))
end


"""
    plot_and_filter_directly_calibrated_traces!(mResfinal, dcps_per_ppb, indices, summedPIs, resultfp)

Plot directly calibrated traces, allow interactive selection of time periods to delete, and set affected traces to NaN in-place.

# Arguments
- `mResfinal`: MeasurementResult (modified in-place)
- `dcps_per_ppb::Matrix`: calibration factors [dcps / ppb]
- `indices::Vector{Int}`: indices of calibrated masses
- `summedPIs::Vector`: summed primary ions
- `resultfp::String`: output file path prefix

# Returns
- `mResfinal` (modified in-place)
"""
function plot_and_filter_directly_calibrated_traces(mResfinal, dcps_per_ppb, indices, summedPIs, resultfp)
    fig = figure(figsize=(10,6))
    ax = subplot(111)
    plot(mResfinal.Times, 1000.0 .* mResfinal.Traces[:, indices] ./ dcps_per_ppb[:, indices])
    plot(mResfinal.Times, summedPIs)
    legStrings = Array{String,1}(undef, length(indices) + 1)
    for (i, idx) in enumerate(indices)
        legStrings[i] = "m/z = " * string(round(mResfinal.MasslistMasses[idx], digits=2)) * " - " * MasslistFunctions.sumFormulaStringFromCompositionArray(mResfinal.MasslistCompositions[:, idx])
    end
    legStrings[end] = "summed primary ions"

    legend(legStrings)
    ylabel("concentration [ppt]")
    xlabel("time [UTC]")
    title("Directly Calibrated Traces")
    yscale("log")
    ylim(1e-2, maximum(summedPIs) * 2)

    savefig(joinpath(resultfp, "DirectlyCalibratedTraces.png"))
    savefig(joinpath(resultfp, "DirectlyCalibratedTraces.pdf"))

    return mResfinal
end
################################# add filter or rename


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

    filterCnr = mResfinal.MasslistCompositions[findfirst(mResfinal.MasslistElements .== "C"), :] .>= 1 #filter for masses with at least one carbon
    filterNoCalib = vec(sum(dcps_per_ppb; dims=1) .> 0) #filter for masses with calibration data
    filterNnr = mResfinal.MasslistCompositions[findfirst(mResfinal.MasslistElements .== "N"), :] .== 1 #filter for masses with exactly one nitrogen; adjust if want to export other masses
    finalfilter = filterCnr .& filterNoCalib .& filterNnr

    calibResult = #contains only filtered masses
        ResultFileFunctions.MeasurementResult(
            mResfinal.Times,
            mResfinal.MasslistMasses[finalfilter],
            mResfinal.MasslistElements,
            mResfinal.MasslistElementsMasses,
            mResfinal.MasslistCompositions[:,finalfilter],
            1000.0 .* mResfinal.Traces[:,finalfilter] ./ dcps_per_ppb[:,finalfilter] #converts traces from cps to ppt
        )

    # filter for interesting traces:
    # either with
    # - findVaryingMasses (often filters too harsh!!!) only for long measurements (10 days) or high Res data
    # - findChangingMasses (if BG and signal times clear!)
    
    c = ResultFileFunctions.findVaryingMasses(calibResult.MasslistMasses,
        calibResult.MasslistCompositions,
        calibResult.Traces;
        sigmaThreshold=3,
        noNitrogen = false,
        onlySaneMasses = false,
        filterCrosstalkMasses=false,
        pointsForSmoothing = 5)
    IndOfinterest = c[1]
    
    #or:
 #= 
    #interactive selection of traces to export
    iifig = PlotFunctions.InteractivePlot(calibResult)
    println("Please select the traces you want to export by scrolling and pressing 'a'.") #how can you terminate the selection?
    PlotFunctions.scrollAddTraces(iifig) ############specify how to end it!
    IndOfinterest = unique(iifig.activeIndices) #indices of selected traces #is empty??????????????
=#   
   
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

    primaryionslist =
        isempty(config.primaryionslist) ?
            massLibrary.FullPrimaryionslist_NH4soft :
            config.primaryionslist

    mResfinal, mResfinal_PIs = load_and_merge_results(config.resultfiles, primaryionslist)
    
    summedPIs = compute_summed_primary_ions(mResfinal_PIs)

    licor_final = interpolate_licor_to_ptr_time(mResfinal, licorDat)

    #mResfinal_PIs TOFTracer2.ResultFileFunctions.MeasurementResult([Dates.DateTime("2025-09-19T20:30:08.188"), Dates.DateTime("2025-09-19T21:00:13.986"), Dates.DateTime("2025-09-19T23:12:30.472"), Dates.DateTime("2025-09-19T23:42:36.350"), Dates.DateTime("2025-09-20T00:12:42.240"), Dates.DateTime("2025-09-20T00:42:48.158"), Dates.DateTime("2025-09-20T01:12:54.075"), Dates.DateTime("2025-09-20T01:42:59.975"), Dates.DateTime("2025-09-20T02:13:05.886"), Dates.DateTime("2025-09-20T02:43:11.798"), Dates.DateTime("2025-09-20T03:13:17.704"), Dates.DateTime("2025-09-20T03:43:23.570"), Dates.DateTime("2025-09-20T04:13:29.466"), Dates.DateTime("2025-09-20T04:43:35.367"), Dates.DateTime("2025-09-20T05:13:41.246"), Dates.DateTime("2025-09-20T05:43:47.149"), Dates.DateTime("2025-09-20T06:13:53.041"), Dates.DateTime("2025-09-20T06:43:58.915"), Dates.DateTime("2025-09-20T07:14:04.811"), Dates.DateTime("2025-09-20T09:36:41.358")], [18.033836, 19.017856000000002, 35.060396, 36.044416, 37.028436, 52.086956, 53.070975999999995, 54.054996, 55.039016000000004, 73.04959600000001], ["C", "C(13)", "H", "H+", "N", "O", "O(18)", "S"], [12.0, 13.00335, 1.00783, 1.007276, 14.00307, 15.99492, 17.99916, 31.97207], [0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0; 3 2 6 5 4 9 8 7 6 8; 1 1 1 1 1 1 1 1 1 1; 1 0 2 1 0 3 2 1 0 0; 0 1 0 1 2 0 1 2 3 4; 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0], [98286.93225551702 20276.200144167084 13760.359393732337 106576.21711008846 206353.05833118974 5.540041578681565 53.20697621822088 16305.664411277161 565756.2346037049 196576.72270328057; 123790.0196691867 60723.86846133821 3245.15757838127 75478.04375116661 467539.09592942306 3.5391954438538464 12.491981154978415 7888.368748199251 678790.6643132545 216272.49430382805; 18015.827960227954 2865.1008557494915 2741.480127539125 13991.721590956107 23711.09529230624 0.4809877639673661 9.316950224671483 1585.8902325967665 64545.153704710494 9388.222005960792; 21387.511024937452 3689.875324073351 2857.491556044269 16989.74439827576 26882.955827285125 0.5243960056863659 9.148252121751717 1764.4113470528175 68545.98069296248 9936.76153393458; 30833.29664193825 3888.599572171985 370.3476643959922 43084.584927050906 20141.763812362813 0.32753753548888315 4.276416564278284 8547.439258190005 86341.2464147376 20784.294766616535; 33052.03927443273 4007.639497826516 234.17739546399454 48230.3250001063 20005.93722532665 0.3118794567298591 3.6639251709734473 10191.343598728705 90700.74953247135 23624.560261809427; 31350.812505128888 3842.8501944890695 195.233580699104 45250.34908992221 19030.108528182143 0.33339537734081914 3.346300388238356 10285.154185245083 90876.96460566131 24206.916246091354; 30158.080163261253 3700.3217575764747 172.9635219693943 42600.258279777176 17670.461266547645 0.3167763786809431 3.0558525023187313 10069.893120442652 87598.55888492538 23900.404172291208; 30885.89695181431 3758.471045890598 168.69299754981205 43822.49767869541 17903.836754627464 0.33639503057535336 3.0679878807494374 10408.814100533273 89505.65519719047 24647.985895331487; 32480.181612844048 3924.233155506358 178.4883698787874 48376.04252127736 19362.368900578444 0.3209545892925484 3.3567818250237256 11429.034784013533 96363.88718063105 26315.508015528474; 33781.817910622914 4079.297255932955 188.46133995908662 52525.32899756504 20799.49514029783 0.3714109918044376 3.40101769703494 12281.29643550142 101674.85711051225 27537.456686626043; 35406.19952983123 4285.610193517204 192.04556570709508 54058.041947783015 21261.191751665734 0.3049318975236086 3.5542341662836607 12531.805564726039 102419.17979075183 27836.991506279275; 35596.730539509685 4355.5191694429495 195.9708059808146 55187.88812618638 21895.639749520404 0.3627022796117051 3.6684833235617287 12898.611235253718 105527.30426664757 28516.576353987377; 36314.37741625003 4449.395340392695 194.37772382139897 54372.763253833014 21410.234741030992 0.31949959414343065 3.5518491647281003 12723.52514893747 102453.15712040327 28022.228814251677; 35119.25185604586 4322.3468048386485 189.1562532418771 52675.14906571881 20838.762493193924 0.3378743160232251 3.6285215075468877 12721.115033377047 102308.5918342868 28199.715586403378; 33760.16272330062 4124.417251783643 172.17976604944275 49084.294533535125 18488.23859199092 0.334594166511454 3.3718331945248647 12186.908661744737 94986.62748252667 26982.600898444853; 37568.103037706765 4604.439812520627 209.3493946867622 57578.73553934395 22187.0856361281 0.3590139652839054 3.786963951750081 13468.662557290369 106591.64134878838 29021.178633988435; 36269.034490577396 4389.131790067037 195.13864132633324 52929.75357799637 20209.856005321213 0.37395029952913583 3.67244309857531 12542.524421983686 99685.18390111456 27616.055907651073; 29408.728433629603 3642.314468626307 155.17007306050883 40602.09876262366 16882.748909363094 0.27143598350746484 3.3781772093847033 11119.775513754348 83559.66190240037 24711.043699274684; 83960.75252792466 4144.455474648014 1116.2603633567353 404892.24862105725 32788.62814256558 3.25818646705125 132.98859987717898 350455.98910356587 657981.7874060586 821971.3850543515])
    #summedPIs [1.8326621251434386e6, 2.509030704542346e6, 216051.35942147192, 242561.60481569121, 339343.57919438276, 364393.51928160747, 354422.6832581542, 339530.2374917781, 347695.88636327593, 374411.6256519539, 396869.869093455, 406180.3984985117, 415204.31532281777, 409990.9501424034, 403119.4102949167, 377572.22567518946, 427701.6923022052, 400918.9019061557, 329912.5875915326, 3.262869031055997e6]
    #licor_final [6.897544056033412, 6.851115702563275, 6.804692668629231, 6.7582678108335354, 6.71184234513227, ...]
    
    calibDF = plot_humidity_dependent_calibration(config.humcalibfile, config.ionization)

    dcps_per_ppb, dcps_per_ppb_err, indices = build_calibration_traces(mResfinal, summedPIs, licor_final, calibDF, hexVSpis_params, config.refName)

    plot_calibration_traces(mResfinal, dcps_per_ppb, summedPIs, indices, config.resultfp)

    plot_and_filter_directly_calibrated_traces(mResfinal, dcps_per_ppb, indices, summedPIs, config.resultfp)

    if config.exportTraces
        export_calibrated_traces(mResfinal, dcps_per_ppb, config.ionization, config.HeaderForExportDict, config.resultfp)
    end
end

#error estimation still missing

end # module CalibrateTraces


