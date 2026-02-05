"""
    module CalibrateTraces

A module for calibrating traces using dry and humidity-dependent calibration data.
"""
module CalibrateTraces

export calibrate_traces_main

import HDF5: ishdf5
using CSV
using DataFrames
using Dates
using PyPlot
using DelimitedFiles
import Statistics: mean, std

using TOFTracer2
import TOFTracer2.CalibrationFunctions as CalF
import TOFTracer2.ImportFunctions as ImpF
import TOFTracer2.InterpolationFunctions as IntpF
import TOFTracer2.ExportFunctions as ExpF
import TOFTracer2.ResultFileFunctions
import TOFTracer2.PlotFunctions
import TOFTracer2.MasslistFunctions
import TOFTracer2.massLibrary

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
    
    # Handle multiple charges -> remove charge notation in parentheses like (4+), (2-) where the number represents charge magnitude
    formula_clean = replace(formula, r"\(\d*[+-]\)" => "")
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
    plot_fit(functiontype)

Return the fitting function corresponding to the given function type.

# Arguments
- `functiontype::String`: Type of fitting function ("power", "linear", "exponential", "double exponential").

# Returns
- `Function`: Corresponding fitting function from TOFTracer2.CalibrationFunctions module.
"""
function plot_fit(functiontype)
        if functiontype == "power"
            return CalF.PowerFunction
        elseif functiontype == "linear"
            return CalF.LinearFunction
        elseif functiontype == "exponential"
            return CalF.Exponential
        #elseif functiontype == "double exponential"
        #    return CalF.DoubleExponential
        else
            error("Unsupported function type: $functiontype")
        end
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
function scatterDryCalibs2(drycalibsfile::String; referenceMasses=refMass, primaryions=primaryionslist) #modified from PlotFunctions.scatterDryCalibs
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
    primaryiontraces = mResDryCalibs.Traces[:,filterarray] * sqrt.(100 ./ primaryionmasses) #duty cycle corrected primary ion traces
    referencetraces = mResDryCalibs.Traces[:,(!).(filterarray)] * sqrt.(100 ./ referenceMasses) #duty cycle corrected sum of reference ion traces
    
    scatter(mResDryCalibs.Times, primaryiontraces, label="sum of primary ions") #plot primary ions summed dcps trace
    scatter(mResDryCalibs.Times, referencetraces, label="reference ion - m/z $(round.(referenceMasses;digits=3))") #this is hexanone
    xlabel("Time")
    ylabel("signals [dcps]")
    title("Dry Calibrations")
    legend()
    yscale("log")
    grid()
    tight_layout()
    savefig("$(dirname(drycalibsfile))\\dryCalibs.png")
    return dryCalibFig, dryCalibAx, mResDryCalibs, primaryiontraces, referencetraces
end


"""
    dryCal_selectPIandRefDataInteractive(drycalibsfile::String, refMass, primaryionslist)

Interactively select dry calibration data points to exclude from fit and return fitted hexanone vs primary ion parameters.

# Arguments
- `drycalibsfile::String`: Path to dry calibration file (hdf5).
- `refMass`: Reference mass for calibration (e.g., hexanone + NH4+).
- `primaryionslist`: List of primary ion masses for calibration.

# Returns
- `hexVSpis_params::Tuple`: Tuple containing fit parameters, errors, and functiontype.
"""
function dryCal_selectPIandRefDataInteractive(drycalibsfile::String, refMass, primaryionslist) #modified from CalibrationFunctions.dryCal_selectPIandRefDataFromIFIG

    dryCalibFig, dryCalibAx, mResDryCalibs, primaryiontraces, referencetraces = scatterDryCalibs2(drycalibsfile; referenceMasses=[refMass], primaryions=primaryionslist)
    # println("please give the minimum y-value to show")
    # dryCalibAx.set_ylim(bottom=parse(Int, readline()))
    println("How many dry calibration data points do you want to include in the fit?")
    nrOfIncludeCalibs = parse(Int, readline())
    while !(nrOfIncludeCalibs in 1:length(mResDryCalibs.Times))
        println("Please enter a valid number between 1 and $(length(mResDryCalibs.Times)).")
        nrOfIncludeCalibs = parse(Int, readline())
    end

    IFIG = PlotFunctions.InteractivePlot(drycalibsfile, dryCalibAx)
    println("Select primary ion calibration coordinates to include by moving the mouse to the respective primary ion data point and press 'c'. Repeat until you have selected $nrOfIncludeCalibs point(s).")
    PlotFunctions.getMouseCoords(IFIG; datetime_x=true)
    while (length(IFIG.coords) < nrOfIncludeCalibs)
        sleep(0.1)
    end
    println("primaryionslist: ", primaryionslist)
    
    #Find closest datapoints
    include_idx = Int[]
    for i in 1:nrOfIncludeCalibs
        t_click  = IFIG.coords[i][1] # ISO 8601 datetime string like "2025-10-11T09:31:27.405"
        y_click  = IFIG.coords[i][2]
        t_click_unix = Dates.datetime2unix(DateTime(t_click))
        t_data_unix = Dates.datetime2unix.(mResDryCalibs.Times)
        dists = (t_data_unix .- t_click_unix).^2 .+ (primaryiontraces .- y_click).^2 # Distance^2 in (Time, PrimaryIonsSum) space
        push!(include_idx, argmin(dists))
    end

    include_idx = unique(include_idx) # in case user clicked very close points, double selections are excluded
    df = DataFrame(
        Time = mResDryCalibs.Times[include_idx], #include selected time points via their Indices
        PrimaryIonsSum = vec(sum(primaryiontraces[include_idx, :]; dims=2)), #sum across masses and convert to vector
        ReferenceSignal = vec(sum(referencetraces[include_idx, :]; dims=2)) #sum across masses and convert to vector
    )

    println("Origin point (0,0) is added to the data for the fit.")
    df = vcat(DataFrame(PrimaryIonsSum=0, ReferenceSignal=0, Time=NaN), df)
    
    figure()
    scatter(df.PrimaryIonsSum, df.ReferenceSignal, label="data", color="gray")
    legend()
    xlabel("sum of primary ions [dcps]")
    ylabel("signal on reference mass [dcps/ppb]") # =dcps since std contains 1ppb of hexanone
    title("Dry Calibration Data Points Selected")

    if length(df.PrimaryIonsSum) < 4 # 4 not 3 because (0,0) is added already
        # If less than 3 calibration points are selected, only a linear fit can be applied, since for more complex fits the covariance matrix cannot be calculated and thus no fit parameters can be obtained. The linear fit is applied automatically in this case, but the user is informed about this.
        userinput = "linear"
        println("Less than 3 calibration points selected, linear fit will be applied.")
        
    else
        # Fit user input
        println("Please choose fit function type: 'linear' or 'power'.") #'double exponential' doesnt work yet, due to covariance matrix issues # 'exponential' fit looks weird
        userinput = readline()
        while !(userinput in ["linear", "power"])
            println("This function is not implemented yet. Please choose fit function type: 'linear' or 'power'.")
            userinput = readline()
        end
    end

    

    functiontype = userinput
    println(df)

    # Fit selected data points
    hexVSpis_params = CalF.fitParameters(df.PrimaryIonsSum, df.ReferenceSignal; functiontype=functiontype)
    nrOfCalibs = nrow(df)
    xforfit = collect(floor(minimum(df.PrimaryIonsSum); sigdigits=1):1000: ceil(maximum(df.PrimaryIonsSum); sigdigits=1))
    fit_func = plot_fit(functiontype)
    fill_between(xforfit,
        fit_func(xforfit, hexVSpis_params[1] .- hexVSpis_params[2] / sqrt(nrOfCalibs)),
        fit_func(xforfit, hexVSpis_params[1] .+ hexVSpis_params[2] / sqrt(nrOfCalibs)),
        label="uncertainty", 
        alpha=0.25)
    plot(xforfit, fit_func(xforfit, hexVSpis_params[1]), label=hexVSpis_params[3])
    legend()
    #yscale("log")
    savefig("$(dirname(drycalibsfile))\\Hexanone_VS_PIs.png")
    savefig("$(dirname(drycalibsfile))\\Hexanone_VS_PIs.pdf")
    return hexVSpis_params
end


"""
    calc_fhex(summedPIs, hexVSpis_params)

Calculate hexanone dry sensitivity vs primary ions and its uncertainty.

# Arguments
- `summedPIs::Vector`: Summed primary ion intensities in dcps.
- `hexVSpis_params::Tuple`: Hexanone vs primary ion parameters.

# Returns
- `f_hex::Vector`: Hexanone dry sensitivity [cps/ppb] vs primary ion cps.
- `f_hex_err::Vector`: Uncertainty of hexanone dry sensitivity.
"""
function calc_fhex(summedPIs, hexVSpis_params)
    #Hexanone dry sensitivity [cps/ppb] vs primary ion cps, calculated for all time points
    f_hex = CalF.applyFunction(summedPIs, hexVSpis_params[1]; functiontype = hexVSpis_params[3][1]) #summedPIs has length times
    f_hex_err = sqrt.( 
        (hexVSpis_params[2][1]./(hexVSpis_params[1][1]))^2 #relative variance of parameter a (slope), (delta a / a)^2
        .+ (log.(summedPIs[summedPIs .> 0]).*hexVSpis_params[2][2]).^2 #relative variance of parameter b (exponent), (ln x *delta b)^2
        .+ 2*(hexVSpis_params[2][1]./(hexVSpis_params[1][1])*log.(summedPIs[summedPIs .> 0]).*hexVSpis_params[2][2]) #covariance term (2 delta a/a * ln x * delta b); 2 cov(a,b) = delta a * delta b is assumed
    )
    f_hex_err[summedPIs .<= 0] .= 0

    total_relative_error = f_hex_err # no humidity dependent calibration error included yet, f_hum_err is zero
    mean_relative_error = mean(total_relative_error)
    std_of_mean_relative_error = std(total_relative_error)
    println("The relative standard error of the calibration factor for this method alone is a factor ",
        round(mean_relative_error,digits=3), " ± ",round(std_of_mean_relative_error,digits=3),
        " due to the error of the normalization to the primary ions.") #change print when f_hum_err != 0

    return f_hex, f_hex_err
end

function ask_for_PIList()
    println("Which primary ions will be used for the calibration?")
    println("press 'f': full list (all water and ammonium clusters)")
    println("press 'w': H3O+ and H2O.H3O+")
    println("press 'a': NH4+ and NH3.NH4+")
    println("press 'o': NH4+ only")

    userinput = readline()
    while !(userinput in ["f", "w", "a", "o"])
        println("Invalid input. Please enter 'f' for full list, 'a' for NH4+ and NH3.NH4+, 'w' for H3O+ and H2O.H3O+, or 'o' for NH4+ only.")
        userinput = readline()
    end
    println("massLibrary: ", massLibrary)
    if userinput == "f"
        primaryionslist = massLibrary.FullPrimaryionslist_NH4soft #adding all possible water and ammonium clusters: [18.033836, 19.017856000000002, 35.060396, 36.044416, 37.028436, 52.086956, 53.070975999999995, 54.054996, 72.065576]
    elseif userinput == "o"
        primaryionslist = MasslistFunctions.massFromComposition(H=3, N=1) #H+ is added automatically in massFromComposition
    elseif userinput == "a"
        primaryionslist = [MasslistFunctions.massFromComposition(H=3, N=1), 
                           MasslistFunctions.massFromComposition(H=6, N=2)] #NH4+ and NH3.NH4+
    elseif userinput == "w"
        primaryionslist = [MasslistFunctions.massFromComposition(H=2, O=1), 
                           MasslistFunctions.massFromComposition(H=4, O=2)] #H2O and H2O.H3O+ clusters only
    end
    return primaryionslist
end

#the rest is modified from calibration script from Wiebke
#################### calibration steps functions #############################
"""
    load_hexVSpis_params(drycalibsfile::String)

Load or create and save hexanone vs. primary ion calibration parameters, using TOFTracer2 with CalibrationFunctions.dryCal_selectPIandRefDataFromIFIG.

# Arguments
- `drycalibsfile::String`: Path to dry calibration file (HDF5 or CSV)
- `refMass`: Reference mass for calibration (e.g., hexanone + NH4+).
- `primaryionslist`: List of primary ion masses for calibration.

# Returns
- `hexVSpis_params::Tuple`: Tuple containing parameters, errors, and metadata

# Saves
- CSV file with hexanone vs. primary ion calibration parameters if loaded from HDF5, in the same directory as `drycalibsfile`.
"""
function load_hexVSpis_params(drycalibsfile::String, refMass, primaryionslist)
    if ishdf5(drycalibsfile)
        hexVSpis_params = dryCal_selectPIandRefDataInteractive(drycalibsfile, refMass, primaryionslist)
        hexVSpis_params2export = vcat(hexVSpis_params[3], ["parameters" "errors"], hcat(hexVSpis_params[1], hexVSpis_params[2]))
        CSV.write("$(dirname(drycalibsfile))\\Hexanone_VS_PIs_params.csv", DataFrame(hexVSpis_params2export, :auto)) 
    else
        a = CSV.read(drycalibsfile, DataFrame, header=[2])
        b = CSV.read(drycalibsfile, DataFrame; footerskip=3, header=false)
        hexVSpis_params = (a.parameters, a.errors, [values(b[1, :])[1], values(b[1, :])[2]])
    end
    return hexVSpis_params
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
        for i in 2:eachindex(resultfiles)
            mResfinal_PIs = ResultFileFunctions.joinResultsTime(mResfinal_PIs, ResultFileFunctions.loadResults(resultfiles[i]; useAveragesOnly = true, massesToLoad = primaryionslist))
            mResfinal_PIs.Traces .= mResfinal_PIs.Traces .* transpose(sqrt.(100 ./ mResfinal_PIs.MasslistMasses)) #duty cycle correction for primary ions
            mResfinal = ResultFileFunctions.joinResultsTime(mResfinal, ResultFileFunctions.loadResults(resultfiles[i]; useAveragesOnly = true))
            mResfinal.Traces .= mResfinal.Traces .* transpose(sqrt.(100 ./ mResfinal.MasslistMasses)) #duty cycle correction for all masses
        end
    end

    #Loaded (58, 1195) traces
    #Loaded and merged 1 result file with a total of 58 time points and 1195 masses.
    
    return mResfinal, mResfinal_PIs
end


"""
    compute_summed_primary_ions(mResfinal_PIs)

Returns the summed primary ions from the measurement results.

# Arguments
- `mResfinal_PIs::struct`: Measurement results for primary ion intensities.

# Returns
- `summedPIs::Vector`: Summed primary ion intensities in dcps.
"""
function compute_summed_primary_ions(mResfinal_PIs)
    summedPIs = vec(sum(mResfinal_PIs.Traces, dims=2)) # sum across masses to get one value per time point
    summedPIs[summedPIs .<= 0] .= 0
    return summedPIs #in dcps
end

################# humidity dependent calibration functions ####################
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
    interpolate_licor_to_ptr_time(mResfinal, licorDat)

Return interpolated Licor humidity data to match the time points of the measurement results.

# Arguments
- `mResfinal::struct`: Measurement results containing time points.
- `licorDat::DataFrame`: DataFrame containing Licor humidity data (datetime, H2O_mmolpermol).

# Returns
- `licor_final::Vector`: Interpolated Licor humidity data aligned with measurement results time points. 
"""
function interpolate_licor_to_ptr_time(mResfinal, licorDat)
    #= return IntpF.sortSelectAverageSmoothInterpolate(
        mResfinal.Times,
        licorDat.datetime,
        licorDat.H2O_mmolpermol;
        returnSTdev = false,
        selectY = [-Inf, Inf] # option to get rid of outliers via the kwarg 'selectY' (use wisely) #[-21, 5]
    )=#
    return IntpF.interpolate(
        mResfinal.Times,
        licorDat.datetime,
        licorDat.H2O_mmolpermol
    )
end


##################################################################

"""
    build_calibration_traces(mResfinal, summedPIs, hexVSpis_params, refName, licor_final, calibDF)
    build_calibration_traces(mResfinal, summedPIs, hexVSpis_params, refName)

Build calibration traces for all compounds. Two methods are available:
- With humidity calibration: Uses humidity-dependent calibration for compounds with 1 oxygen atom, dry calibration for >=2 oxygen and undefined compounds, and compounds from the gas standard individually humidity dependent.
- Dry calibration only: Applies dry calibration to all compounds except 0 oxygen compounds.

# Arguments
- `mResfinal::MeasurementResult`: Measurement results containing mass list and time points.
- `summedPIs::Vector`: Summed primary ion intensities in dcps.
- `hexVSpis_params::Tuple`: Hexanone vs primary ion parameters.
- `refName::String`: Reference compound name.
- `licor_final::Vector`: (Humidity dependent method only) Interpolated Licor humidity data aligned with measurement results time points.
- `calibDF::DataFrame`: (Humidity dependent method only) Humidity dependent calibration data.

# Returns
- `dcps_per_ppb::Matrix`: Calibration traces in dcps per ppb for all compounds.
- `indices::Vector{Int}`: (Humidity dependent method only) Indices of individually calibrated compounds.
"""
function build_calibration_traces(mResfinal, summedPIs, hexVSpis_params, refName, licor_final, calibDF)

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
    
    #initialize calibration traces matrix with zeros
    dcps_per_ppb = zeros(length(mResfinal.Times), length(mResfinal.MasslistMasses))
    dcps_per_ppb_err = zeros(size(dcps_per_ppb))

    # dry sensitivity of hexanone vs PIs
    f_hex, f_hex_err = calc_fhex(summedPIs, hexVSpis_params)

    # wet sensitivity of hexanone vs AH
    f_hum = CalF.applyFunction(licor_final, ref_params; functiontype = "double exponential") #unitless # licor_final has length of times #use licor_final instead of CalF.applyFunction(fpfinal, humparams[1]; functiontype = humparams[3][1]) only if icor data is complete
    println("ref_params: ", ref_params)
    println("f_hum: ", f_hum)
    println("licor_final: ", licor_final)
    #1) dry / kinetic limit calibration for undef and >=2 O
    println("calibrating all compounds with >=2 oxygen atoms and undefined ones with reference $(refName) dry.")
    dcps_per_ppb[:, (undeffilter .| twoplusoxygenfilter)] .= f_hex # f_hum0 not needed since f_hum0 = ones(length(mResfinal.Times)) # 1 because normalized to dry point of humidity dependent calibration of hexanone
    dcps_per_ppb_err[:, (undeffilter .| twoplusoxygenfilter)] .= dcps_per_ppb[:, (undeffilter .| twoplusoxygenfilter)] .* f_hex_err # no humidity dependent calibration error for dry calibration

    #2) humid / equilibrium calibration for 1 O
    println("calibrating all compounds with 1 oxygen atom humidity-dependent with reference $(refName).")
    dcps_per_ppb[:, oneoxygenfilter] .= f_hex .* f_hum
    dcps_per_ppb_err[:, oneoxygenfilter] .= dcps_per_ppb[:, oneoxygenfilter] .* sqrt.(f_hex_err.^2 .+ 0.0^2) # combine errors from dry and humid calibration

    #3) individual calib (hum dep) for the compounds in the gas standard -> hum dep calibration file.
    indices = Int[]
    for row in eachrow(calibDF) 
        params = row[[:p1, :p2, :p3, :p4, :p5]]
        index = findfirst(isapprox.(mResfinal.MasslistMasses, row.Mass; atol=0.0001)) #find masses in masslist matching compounds in calibDF (humidity dependent calibration data)
        f_hum_row = CalF.applyFunction(licor_final, params; functiontype = "double exponential") 
        
        if index isa Int64
            dcps_per_ppb[:, index] = f_hex .* f_hum_row
            dcps_per_ppb_err[:, index] = dcps_per_ppb[:, index] .* sqrt.(f_hex_err.^2 .+ 0.0^2) # combine errors from dry and humid calibration
            push!(indices, index)
        end
    end    

    return dcps_per_ppb, indices #indices of gas standard compounds calibrated individually
end
#what I changed: remove frostpoint dependency, use only licor instead. use humidity dependant calibration for 1 oxygen and kinetic limit calibration for 2+ oxygen.
#added error here, but f_hum_err is still zero
#compounds with 0 oxygen are ignored (zero sensitivity)

function build_calibration_traces(mResfinal, summedPIs, hexVSpis_params, refName)

    undeffilter = BitVector(sum(mResfinal.MasslistCompositions; dims=1)[1, :] .== 0)
    oxygen_number = findfirst(==("O"), mResfinal.MasslistElements)
    oneplusoxygenfilter = BitVector(mResfinal.MasslistCompositions[oxygen_number, :] .>= 1)
    println("found $(sum(undeffilter)) undefined masses, and $(sum(oneplusoxygenfilter)) masses with >=1 oxygen atom")

    #initialize calibration traces matrix with zeros
    dcps_per_ppb = zeros(length(mResfinal.Times), length(mResfinal.MasslistMasses))
    dcps_per_ppb_err = zeros(size(dcps_per_ppb))

    f_hex, f_hex_err = calc_fhex(summedPIs, hexVSpis_params)

    println("calibrating all compounds (except 0 oxygen compounds) dry with reference $(refName).")
    dcps_per_ppb[:, (undeffilter .| oneplusoxygenfilter)] .= f_hex
    dcps_per_ppb_err[:, (undeffilter .| oneplusoxygenfilter)] .= dcps_per_ppb[:, (undeffilter .| oneplusoxygenfilter)] .* f_hex_err # no humidity dependent calibration error for dry calibration

    return dcps_per_ppb
end


"""
    plot_calibration_traces(mResfinal, dcps_per_ppb, summedPIs, indices, resultfp)

Plot calibration traces for selected compounds and save the plots. 

# Arguments
- `mResfinal::struct`: Measurement results containing mass list and time points.
- `dcps_per_ppb::Matrix`: Calibration traces in dcps per ppb for all compounds.
- `summedPIs::Vector`: Summed primary ion intensities in dcps.
- `indices::Vector{Int}`: Indices of calibrated compounds.
- `resultfp::String`: File path to save the plots.

# Saves
- Calibration traces plot as PNG and PDF files.
"""
function plot_calibration_traces(mResfinal, dcps_per_ppb, summedPIs, indices, resultfp, ionization)
    fig, ax = subplots(figsize = (10, 6))

    ax.set_xlabel("time [UTC]")
    ax.set_ylabel("calibration factor [dcps / ppb]")
    ax.set_yscale("log")
    
    for i in indices
        ax.plot(mResfinal.Times, dcps_per_ppb[:, i],
            label = "$(round(mResfinal.MasslistMasses[i],digits=2)) * $(MasslistFunctions.sumFormulaStringFromCompositionArray(mResfinal.MasslistCompositions[:,i], ion=ionization))"
        )
    end

    ax.plot(mResfinal.Times, summedPIs, label="summed primary ions [dcps]", color="black", linestyle="--") 

    ax.set_title("Calibration Traces")
    ax.legend(loc="upper right")
   
    fig.savefig(joinpath(resultfp, "CalibrationTraces.png"))
    fig.savefig(joinpath(resultfp, "CalibrationTraces.pdf"))
end


"""
    plot_directly_calibrated_traces(mResfinal, dcps_per_ppb, summedPIs, indices, resultfp)

Plot directly calibrated traces.

# Arguments
- `mResfinal`: MeasurementResult (modified in-place)
- `dcps_per_ppb::Matrix`: calibration factors [dcps / ppb]
- `indices::Vector{Int}`: indices of calibrated masses
- `summedPIs::Vector`: summed primary ions
- `resultfp::String`: output file path prefix

# Saves
- Directly calibrated traces plot as PNG and PDF files.
"""
function plot_directly_calibrated_traces(mResfinal, dcps_per_ppb, summedPIs, indices, resultfp, ionization)
    fig = figure(figsize=(10,6))
    ax = subplot(111)
    ax.set_yscale("log")
    ax.plot(mResfinal.Times, 1000.0 .* mResfinal.Traces[:, indices] ./ dcps_per_ppb[:, indices])
    ax.plot(mResfinal.Times, summedPIs, color="black", linestyle="--") 
    
    legStrings = Array{String,1}(undef, length(indices)+1)
    for (i, idx) in enumerate(indices)
        legStrings[i] = "m/z = " * string(round(mResfinal.MasslistMasses[idx], digits=2)) * " - " * MasslistFunctions.sumFormulaStringFromCompositionArray(mResfinal.MasslistCompositions[:, idx], ion=ionization)
    end
    legStrings[end] = "summed primary ions [dcps]"
    
    ax.set_title("Directly Calibrated Traces")
    ax.legend(legStrings, loc="upper right")
    ax.set_ylabel("concentration [ppt]")
    ax.set_xlabel("time [UTC]")
    #ax.set_ylim(1e-2, maximum(summedPIs) * 2)

    savefig(joinpath(resultfp, "DirectlyCalibratedTraces.png"))
    savefig(joinpath(resultfp, "DirectlyCalibratedTraces.pdf"))
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


###############################################################################
# Main Script
###############################################################################

"""
    calibrate_traces_main(dir_licor_data, humcalibfile, drycalibsfile, resultfp, resultfiles, ionization, refName, exportTraces, HeaderForExportDict)
    calibrate_traces_main(drycalibsfile, resultfp, resultfiles, ionization, refName, exportTraces, HeaderForExportDict)

Main function to calibrate measurement traces. Two methods are available:
- With humidity-dependent calibration: Applies humidity-dependent calibration for 1-oxygen compounds and individually calibrated gas standard compounds. Requires Licor humidity data and humidity calibration parameters.
- Dry calibration only: Applies dry calibration to all compounds (except 0-oxygen compounds). Suitable for temperatures >0°C or when humidity data is unavailable.

# Arguments
- `dir_licor_data::String`: (Humidity method only) Directory path containing Licor data files.
- `humcalibfile::String`: (Humidity method only) File path for humidity-dependent calibration parameters (txt format).
- `drycalibsfile::String`: File path for dry calibration data (HDF5) or hexanone vs. primary ion parameters (CSV).
- `resultfp::String`: Directory path for saving calibration results and plots.
- `resultfiles::Vector{String}`: List of measurement result file paths to be calibrated (HDF5).
- `ionization::String`: Ionization type used in measurements (e.g., "NH4+", "H+").
- `refMass`: Reference compound mass for calibration (e.g., hexanone + NH4+).
- `refName::String`: Reference compound chemical formula (e.g., "C6H12O.NH4+").
- `exportTraces::Bool`: Flag to export calibrated traces to CLOUD format CSV.
- `HeaderForExportDict::Dict{String,Any}`: Dictionary containing metadata for exported CSV headers.

# Interactive Prompts
During execution, the function will prompt for:
- Primary ion selection (full list or NH4+ only)
- Number of dry calibration points to exclude
- Fit function type (linear, power)

# Saves
- Dry calibration plots: `dryCalibs.png`, `Hexanone_VS_PIs.png/pdf` in dry calibration directory
- (Humidity-dependent only) Calibration traces plot for gas standard compounds: `CalibrationTraces.png/pdf` in `resultfp`
- (Humidity-dependent only) Directly calibrated traces for gas standard compounds: `DirectlyCalibratedTraces.png/pdf` in `resultfp`
- Hexanone vs. PI parameters: `Hexanone_VS_PIs_params.csv` (if processing HDF5 dry calibration file)
- Exported traces: `CLOUD_PTR_<ionization>_ambient_vXX.txt` in `resultfp` (if `exportTraces=true`)
"""
function calibrate_traces_main(dir_licor_data, humcalibfile, drycalibsfile, resultfp, resultfiles, ionization, refMass, refName, exportTraces, HeaderForExportDict)
    
    primaryionslist = ask_for_PIList()
    hexVSpis_params = load_hexVSpis_params(drycalibsfile, refMass, primaryionslist)
    mResfinal, mResfinal_PIs = load_and_merge_results(resultfiles, primaryionslist) #PI: (58, 1) #rest of masses: (58, 1195)
    summedPIs = compute_summed_primary_ions(mResfinal_PIs)

    licorDat = load_licor_data(dir_licor_data)
    licor_final = interpolate_licor_to_ptr_time(mResfinal, licorDat)
    calibDF = CSV.read(humcalibfile, DataFrame; delim='\t', header=2) #DataFrame, containing humidity dependent calibration parameters (p1-p5 with errors, double exponential fit) for different compounds including reference compound
    
    dcps_per_ppb, indices = build_calibration_traces(mResfinal, summedPIs, hexVSpis_params, refName, licor_final, calibDF)
    
    plot_calibration_traces(mResfinal, dcps_per_ppb, summedPIs, indices, resultfp, ionization) #for directly calibrated masses
    plot_directly_calibrated_traces(mResfinal, dcps_per_ppb, summedPIs, indices, resultfp, ionization)

    if exportTraces
        export_calibrated_traces(mResfinal, dcps_per_ppb, ionization, HeaderForExportDict, resultfp)
    end
end

function calibrate_traces_main(drycalibsfile, resultfp, resultfiles, ionization, refMass, refName, exportTraces, HeaderForExportDict)
    primaryionslist = ask_for_PIList()
    hexVSpis_params = load_hexVSpis_params(drycalibsfile, refMass, primaryionslist)
    mResfinal, mResfinal_PIs = load_and_merge_results(resultfiles, primaryionslist) #PI: (58, 1) #rest of masses: (58, 1195)
    summedPIs = compute_summed_primary_ions(mResfinal_PIs)

    dcps_per_ppb = build_calibration_traces(mResfinal, summedPIs, hexVSpis_params, refName)

    if exportTraces
        export_calibrated_traces(mResfinal, dcps_per_ppb, ionization, HeaderForExportDict, resultfp)
    end
end


end # module CalibrateTraces


