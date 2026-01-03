module Calibrate_Traces

#export ....

#using ... from ...
using HDF5
using CSV
using DataFrames

using ToFTracer2
import TOFTracer2.CalibrationFunctions as CalF
import TOFTracer2.ImportFunctions as ImpF

################################
# define relative filepaths of used data
################################

dir_CLOUD18 = joinpath(@__DIR__, "..", "..", "..")
dir_data = joinpath(dir_CLOUD18, "CLOUD18_data")
dir_calib_data = joinpath(dir_data, "Calibration")
dir_licor_data = joinpath(dir_data, "Licor")

#dry calibration file
drycalibsfile = joinpath(dir_calib_data, "2025-11-25 08h28m48_1ppb_std_brown.h5") #FILENAME ANPASSEN!!! WIE SIEHT DIESE FILE AUS???

#humidity dependent calibration file
dir_humidcalib = joinpath(dir_calib_data, "Humidity-dependent_std", "results")
humidcalibsfile = joinpath(dir_humidcalib, "_result.hdf5") #STIMMT FILENAME??? WAS FÜR EINE FILE SOLL DAS SEIN; EHER txt ???

ionization = "NH4+" # "NH4+", "H+"...
primaryionslist = [] # leave empty -> default: adding all possible water and ammonium clusters


refMass = massLibrary.HEXANONE_nh4[1] #mass of hexanone + NH4+ from julia package manualMassLibrary.jl
refName = TOFTracer2.MasslistFunctions.sumFormulaStringFromCompositionArray(massLibrary.HEXANONE_nh4[4]; ion = "")

#file to be calibrated at once with same mass list
dir_results_to_calibrate = joinpath(dir_calib_data, "Test")
resultsfiles_to_calibrate = ["$(dir_results_to_calibrate)/results/_result.hdf5"] #adjust filename, can add multiple files

exportTraces = true # if true, check HeaderForExportDict below:
HeaderForExportDict = Dict(
        "title"=>"oVOCs from Runs...",
        "level"=>2,
        "version"=>"01",
        "authorname_mail"=>"Ruth, Clea clea.ruth@student.uibk.ac.at",
        "units"=>"ppt",
        "addcomment"=>"The data have been humidity-depently calibrated with Hexanone as reference (Onr=1),
        compounds with Onr>1 are calibrated with kinetic limit.
        All traces have been corrected to the duty-cycle-corrected primary ion trace.
        Uncertainty roughly factor 3. Not transmission-corrected yet.\n",
        "threshold"=>0,
        "nrrows_addcomment" => 4
        )


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


if !isdefined(Main,:hexVSpis_params)
    hexVSpis_params = load_hexVSpis_params(drycalibsfile)
end

#load humidity data inlet
function load_licor_data(dir_licor_data::String)
    return ImpF.createLicorData_fromFiles(dir_licor_data;
        filefilter=r"licor_.*\.txt", #rename file to use with licor_restofthename.txt
        headerrow=2,
        columnNameOfInterest="H₂O_(mmol_mol⁻¹)",
        type_columnOfInterest=Float64)
end

if !isdefined(Main,:licorDat)
    licorDat = load_licor_data(dir_licor_data)
end

#interpolate licor data to ptr time
function interpolate_licor_to_ptr_time(mResfinal, licorDat)
    licor_final =
        IntpF.sortSelectAverageSmoothInterpolate(
            mResfinal.Times,
            licorDat.time,
            licorDat.value;
            returnSTdev = false,
            selectY = [0.0, Inf]
        )
    return licor_final
end

#plot humidity dependent calibration
function plot_humidity_dependent_calibration(humcalibfp, humparams, licorDat, ionization)
    calibDF = CSV.read(humcalibfp, DataFrame; header=2)

    CalF.plot_humdep_fromCalibParameters(
        calibDF = calibDF,
        humparams = humparams,
        cloudhum = licorDat.value,   # uses now Licor directly
        hum4plot = collect(0:0.2:12),
        savefp = dirname(humcalibfp),
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

    licor_final =
        IntpF.sortSelectAverageSmoothInterpolate(
            mResfinal.Times,
            licorDat.time,
            licorDat.value;
            returnSTdev = false,
            selectY = [-21, 5]
        )

    #initialize calibration traces
    dcps_per_ppb =
        zeros(length(mResfinal.Times), length(mResfinal.MasslistMasses))

    fref = findfirst(calibDF[!, "Sumformula"] .== refName)
    refparams = calibDF[fref, [:p1, :p2, :p3, :p4, :p5]]

    undeffilter =
        vec(sum(mResfinal.MasslistCompositions; dims=1)) .== 0
    greater2oxygenfilter =
        vec(mResfinal.MasslistCompositions[
            mResfinal.MasslistElements .== "O", :] .>= 2)

    #dry / kinetic limit calibration for undef and >=2 O
    dcps_per_ppb[:, (undeffilter .| greater2oxygenfilter)] .=
        CalF.applyFunction(
            summedPIs,
            hexVSpis_params[1];
            functiontype = hexVSpis_params[3][1]
        ) .*
        CalF.applyFunction(
            zeros(length(mResfinal.Times)),
            refparams;
            functiontype = "double exponential"
        )

    #humid / equilibrium calibration for 1 O
    dcps_per_ppb[:,
        vec(mResfinal.MasslistCompositions[
            mResfinal.MasslistElements .== "O", :] .== 1)
    ] .=
        CalF.applyFunction(
            summedPIs,
            hexVSpis_params[1];
            functiontype = hexVSpis_params[3][1]
        ) .*
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

    for (name, mass) in zip(calibDF[!, "Sumformula"], calibDF[!, "Mass"])
        params =
            calibDF[
                findfirst(calibDF[!, "Sumformula"] .== name),
                [:p1, :p2, :p3, :p4, :p5]
            ]

        index =
            findfirst(
                isapprox.(mResfinal.MasslistMasses, mass, atol=0.0001)
            )

        if index isa Int
            dcps_per_ppb[:, index] =
                CalF.applyFunction(
                    summedPIs,
                    hexVSpis_params[1];
                    functiontype = hexVSpis_params[3][1]
                ) .*
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

function main(config::CalibrationConfig)
    hexVSpis_params = load_hexVSpis_params(config.drycalibsfile)

    licorDat = load_licor_data(config.licorFilepath)

    calibDF =
        plot_humidity_dependent_calibration(
            config.humcalibfp,
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