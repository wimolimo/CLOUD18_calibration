#main TO Do: find out what to do about humparams and cloudhum, where to get it from.
module Calibrate_Traces

# export calibrate_traces_main(CalibrationConfig)

#using ... from ...
using HDF5
using CSV
using DataFrames
using Dates
using Plots
using DelimitedFiles
#import Statistics

#check !
using ToFTracer2
import TOFTracer2.CalibrationFunctions as CalF
import TOFTracer2.ImportFunctions as ImpF
import TOFTracer2.InterpolationFunctions as IntpF
import TOFTracer2.ExportFunctions as ExpF
import TOFTracer2.ResultFileFunctions
import TOFTracer2.PlotFunctions
import TOFTracer2.MasslistFunctions
import TOFTracer2.massLibrary

################################
# define relative filepaths of used data
################################

struct CalibrationConfig
    dir_licor_data::String
    humcalibfile::String
    drycalibsfile::String
    resultfp::String
    resultfiles::Vector{String}
    ionization::String
    primaryionslist #::Vector{Float64}
    refMass #::Float64
    refName::String
    exportTraces::Bool
    HeaderForExportDict::Dict{String,Any}
end

dir_CLOUD18 = joinpath(@__DIR__, "..", "..", "..")
dir_calib_data = joinpath(dir_CLOUD18, "CLOUD18_data", "Calibration")
dir_licor_data = joinpath(dir_CLOUD18, "CLOUD18_data", "Licor")

#dry calibration file
drycalibsfile = joinpath(dir_calib_data, "2025-11-25 08h28m48_1ppb_std_brown.h5") #FILENAME ANPASSEN!!! WIE SIEHT DIESE FILE AUS???

#humidity dependent calibration file
humcalibfile = joinpath(dir_calib_data, "Humidity-dependent_std", "results", "_result.hdf5") #humcalibfp; should be txt!!!! #from humidity dependence calibration script?

#file to be calibrated at once with same mass list
resultfp = joinpath(dir_calib_data, "Test") #change result filepath to data that is analyzed, results of this script are also saved here
resultfiles = ["$(resultfp)/results/_result.hdf5"] #adjust filename, can add multiple files #["$(resultfp)part1/results/_result.hdf5","$(resultfp)part2/results/_result.hdf5"]

ionization = "NH4+" # "NH4+", "H+"...
primaryionslist = [] # leave empty -> default: adding all possible water and ammonium clusters

refMass = massLibrary.HEXANONE_nh4[1] #mass of hexanone + NH4+ from julia package manualMassLibrary.jl
refName = TOFTracer2.MasslistFunctions.sumFormulaStringFromCompositionArray(massLibrary.HEXANONE_nh4[4]; ion = "")

exportTraces = true # if true, check HeaderForExportDict below:
HeaderForExportDict = Dict(
        "title"=>"oVOCs from Runs...",
        "level"=>2,
        "version"=>"01",
        "authorname_mail"=>"Ruth, Clea clea.ruth@student.uibk.ac.at",
        "units"=>"ppt",
        "addcomment"=>"The data have been humidity-depently calibrated with Hexanone as reference (Onr=1), compounds with Onr>1 are calibrated with kinetic limit. All traces have been corrected to the duty-cycle-corrected primary ion trace. Uncertainty roughly factor 3. Not transmission-corrected yet.\n",
        "threshold"=>0,
        "nrrows_addcomment" => 4
        )

#= humparams is not defined yet
humparams=(Float64[],Float64[]," ")
humparams = fitParameters(cloudhumfinal, licorfinal; functiontype="exponential")
cloudhumfinal = IntpF.sortSelectAverageSmoothInterpolate(time2interpolate2, cloudhum.time, cloudhum[!,cloudhumLabel]; returnSTdev=false, selectY=selectY.cloud)
licorfinal = IntpF.sortSelectAverageSmoothInterpolate(time2interpolate2, licorDat.datetime, licorDat.H2O_mmolpermol; returnSTdev=false, selectY=selectY.inlet)
cloudhum = CSV.read(frostpointfile, DataFrame; dateformat=frostpointDatetimeFormat)
    if Year(cloudhum.time[1]) < Year(999)
        cloudhum.time = cloudhum.time .+ Year(2000)
    end =#

# 
#####################################################
# load and prepare metadata of the final calibration
#####################################################

"""
    load_hexVSpis_params(drycalibsfile::String)

Load or create and save hexanone vs. primary ion calibration parameters, using TOFTracer2 with CalibrationFunctions.dryCal_selectPIandRefDataFromIFIG.

# Arguments
- `drycalibsfile::String`: Path to dry calibration file (HDF5 or CSV)

# Returns
- `hexVSpis_params`: Tuple containing parameters, errors, and metadata

# Creates
- CSV file with hexanone vs. primary ion calibration parameters if loaded from HDF5
"""
function load_hexVSpis_params(drycalibsfile::String)
    if HDF5.ishdf5(drycalibsfile)
        hexVSpis_params = CalF.dryCal_selectPIandRefDataFromIFIG(drycalibsfile)
        hexVSpis_params2export =
            vcat(hexVSpis_params[3], ["parameters" "errors"],
                 hcat(hexVSpis_params[1], hexVSpis_params[2]))
        writedlm("$(dirname(drycalibsfile))Hexanone_VS_PIs_params.csv",
                 hexVSpis_params2export)
    else
        a = CSV.read(drycalibsfile, DataFrame, header=[2])
        b = CSV.read(drycalibsfile, DataFrame; footerskip=3, header=false)
        hexVSpis_params = (a.parameters, a.errors, [values(b[1, :])[1], values(b[1, :])[2]])
    end
    return hexVSpis_params
end


#load humidity data inlet
function load_licor_data(dir_licor_data::String)
    return ImpF.createLicorData_fromFiles(dir_licor_data;
        filefilter=r"licor_.*\.txt", #rename file to use with licor_restofthename.txt
        headerrow=2,
        columnNameOfInterest="H₂O_(mmol_mol⁻¹)",
        type_columnOfInterest=Float64)
end


#interpolate licor data to ptr time
function interpolate_licor_to_ptr_time(mResfinal, licorDat)
    licor_final =
        IntpF.sortSelectAverageSmoothInterpolate(
            mResfinal.Times,
            licorDat.time,
            licorDat.value;
            returnSTdev = false,
            selectY = [0.0, Inf] #[-21, 5]
        )
    return licor_final
end


#plot humidity dependent calibration
function plot_humidity_dependent_calibration(humcalibfile, humparams, licorDat, ionization)
    calibDF = CSV.read(humcalibfile, DataFrame; header=2)

    CalF.plot_humdep_fromCalibParameters(
        calibDF = calibDF,
        humparams = humparams,
        cloudhum = licorDat.value,   # uses now Licor directly #TO DO check if ok!!!
        hum4plot = collect(0:0.2:12),
        savefp = dirname(humcalibfile),
        humdepcalibRelationship = "double exponential",
        humidityRelationship = "exponential",
        ionization = ionization
    )

    return calibDF
end



function load_and_merge_results(resultfiles, primaryionslist)
    mResfinal_PIs =
        ResultFileFunctions.loadResults(
            resultfiles[1];
            useAveragesOnly = true,
            massesToLoad = primaryionslist
        )

    mResfinal =
        ResultFileFunctions.loadResults(
            resultfiles[1];
            useAveragesOnly = true
        )

    if length(resultfiles) > 1
        for i in 2:length(resultfiles)
            mResfinal_PIs =
                ResultFileFunctions.joinResultsTime(
                    mResfinal_PIs,
                    ResultFileFunctions.loadResults(
                        resultfiles[i];
                        useAveragesOnly = true,
                        massesToLoad = primaryionslist
                    )
                )

            mResfinal =
                ResultFileFunctions.joinResultsTime(
                    mResfinal,
                    ResultFileFunctions.loadResults(
                        resultfiles[i];
                        useAveragesOnly = true
                    )
                )
        end
    end

    return mResfinal, mResfinal_PIs
end

function compute_summed_primary_ions(mResfinal_PIs)
    summedPIs =
        mResfinal_PIs.Traces .* sqrt.(100 ./ mResfinal_PIs.MasslistMasses)
    summedPIs[summedPIs .<= 0] .= 0
    return summedPIs
end

function build_calibration_traces(
    mResfinal,
    summedPIs,
    licor_final,
    humparams,
    calibDF,
    hexVSpis_params,
    refName)

    #initialize calibration traces
    dcps_per_ppb = zeros(length(mResfinal.Times), length(mResfinal.MasslistMasses))

    fref = findfirst(calibDF[!, "Sumformula"] .== refName)
    refparams = calibDF[fref, [:p1, :p2, :p3, :p4, :p5]]

    #composition based filters
    undeffilter = BitVector(sum(mResfinal.MasslistCompositions; dims=1)[1, :] .== 0)

    oxygen_number = findfirst(==("O"), mResfinal.MasslistElements)

    oneoxygenfilter = BitVector(mResfinal.MasslistCompositions[oxygen_number, :] .== 1)

    twoplusoxygenfilter = BitVector(mResfinal.MasslistCompositions[oxygen_number, :] .>= 2)
    
    #hexanone-primary-ion factor
    f_hex = CalF.applyFunction(summedPIs, hexVSpis_params[1]; functiontype = hexVSpis_params[3][1])

    #dry / kinetic limit calibration for undef and >=2 O
    dcps_per_ppb[:, (undeffilter .| twoplusoxygenfilter)] .=
        f_hex .*
        CalF.applyFunction(
            zeros(length(mResfinal.Times)),
            refparams;
            functiontype = "double exponential"
        )

    #humid / equilibrium calibration for 1 O
    dcps_per_ppb[:, oneoxygenfilter] .=
        f_hex .*
        CalF.applyFunction(
            CalF.applyFunction(
                licor_final,
                humparams[1];
                functiontype = humparams[3][1]
            ),
            refparams;
            functiontype = "double exponential"
        )

    indices = Int[]

 #= for (name, mass) in zip(calibDF[!, "Sumformula"], calibDF[!, "Mass"])
        params =
            calibDF[
                findfirst(calibDF[!, "Sumformula"] .== name),
                [:p1, :p2, :p3, :p4, :p5]
            ]

        index =
            findfirst(
                isapprox.(mResfinal.MasslistMasses, mass, atol=0.0001)
            )
 =#
    for row in eachrow(calibDF)
    params = row[[:p1, :p2, :p3, :p4, :p5]]

    index = findfirst(isapprox.(mResfinal.MasslistMasses, row.Mass; atol=1e-4))

        if index isa Int
            dcps_per_ppb[:, index] =
                f_hex .*
                CalF.applyFunction(
                    CalF.applyFunction(
                        licor_final,
                        humparams[1];
                        functiontype = humparams[3][1]
                    ),
                    params;
                    functiontype = "double exponential"
                )

            push!(indices, index)
        end
    end

    return dcps_per_ppb, indices
end
#what I changed: remove frostpoint dependency, use only licor instead. use humidity dependant calibration for 1 oxygen and kinetic limit calibration for 2+ oxygen.
#TO DO: don't show 0 oxygen. --> later?

function plot_calibration_traces(mResfinal, dcps_per_ppb, summedPIs, indices, resultfp)
    p1 = plot(size=(1000,600),
        xlabel="time [UTC]",
        ylabel="calibration factor [dcps / ppb]")

    for i in indices
        plot!(p1,
            mResfinal.Times,
            dcps_per_ppb[:, i],
            label = "index $i - $(round(mResfinal.MasslistMasses[i],digits=2)), " *
                    MasslistFunctions.sumFormulaStringFromCompositionArray(
                        mResfinal.MasslistCompositions[:,i]
                    )
        )
    end

    plot!(p1, mResfinal.Times, summedPIs, label="summed primary ions", lw=2)

    savefig(p1, "$(resultfp)CalibrationTraces.png")
    savefig(p1, "$(resultfp)CalibrationTraces.pdf")
end

function export_calibrated_traces(
    mResfinal,
    dcps_per_ppb,
    indices,
    ionization,
    HeaderForExportDict,
    resultfp
)
    filterCnr =
        mResfinal.MasslistCompositions[
            findfirst(mResfinal.MasslistElements .== "C"), :
        ] .>= 1

    #masses with 0 oxygen are filtered out
    filterNoCalib = vec(sum(dcps_per_ppb; dims=1) .> 0)
    filterNnr =
        mResfinal.MasslistCompositions[
            findfirst(mResfinal.MasslistElements .== "N"), :
        ] .== 1

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



function calibrate_traces_main(config::CalibrationConfig)

    hexVSpis_params =
        config.hexVSpis_params === nothing ?
            load_hexVSpis_params(config.drycalibsfile) :
            config.hexVSpis_params

    licorDat =
        config.licorDat === nothing ?
            load_licor_data(config.dir_licor_data) :
            config.licorDat

    calibDF =
        plot_humidity_dependent_calibration(
            config.humcalibfile,
            humparams,
            licorDat,
            config.ionization
        )

    primaryionslist =
        isempty(config.primaryionslist) ?
            massLibrary.FullPrimaryionslist_NH4soft :
            config.primaryionslist

    mResfinal, mResfinal_PIs =
        load_and_merge_results(config.resultfiles, primaryionslist)

    summedPIs = compute_summed_primary_ions(mResfinal_PIs)

    licor_final = interpolate_licor_to_ptr_time(mResfinal, licorDat)

    dcps_per_ppb, indices =
        build_calibration_traces(
            mResfinal,
            summedPIs,
            licor_final,
            humparams,
            calibDF,
            hexVSpis_params,
            config.refName
        )

    plot_calibration_traces(
        mResfinal,
        dcps_per_ppb,
        summedPIs,
        indices,
        config.resultfp
    )

    if config.exportTraces
        export_calibrated_traces(
            mResfinal,
            dcps_per_ppb,
            indices,
            config.ionization,
            config.HeaderForExportDict,
            config.resultfp
        )
    end
end



end # module Calibrate_Traces