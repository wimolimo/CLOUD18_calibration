#julia
#] activate .
#] precompile
#include("src\\run_CLOUD18_calibration.jl")

include("CLOUD18_calibration.jl")
using CLOUD18_calibration: calibrate_traces_main
using TOFTracer2
using TOFTracer2.MasslistFunctions
using TOFTracer2: massLibrary
using Dates
using CSV
using DataFrames

#################################
# Define parameters for calibration
#################################

dir_CLOUD18 = joinpath(@__DIR__, "..", "..")
dir_calib_data = joinpath(dir_CLOUD18, "CLOUD18_data", "Calibration")
dir_licor_data = joinpath(dir_CLOUD18, "CLOUD18_data", "Licor") #select licor data where times match the measurement to be calibrated: rename file to use with licor_restofthename.txt

#dry calibration file # this ican be either the processed file of the dry calibrations or the CSV file containing the exported hexanone vs primary ion parameters for loading them:
drycalibsfile = joinpath(dir_calib_data, "dry_std", "results", "_result.hdf5") 
#or #drycalibsfile = joinpath(dir_calib_data, "dry_std\\results\\Hexanone_VS_PIs_params.csv")

#files for humidity dependent calibration and resulting file with fit parameters
std_fp = joinpath(dir_calib_data, "Humidity-dependent_std", "results") #folder containing humidity dependent calibration data
std_file = joinpath(std_fp, "_result.hdf5") #file to be used for humidity dependent calibration
bg_file = joinpath(dir_calib_data, "humidity_dependent_BG", "results", "_result.hdf5") #background file for humidity dependent calibration
hum_file = joinpath(dir_licor_data, "2025-11-21.txt") #licor file for humidity dependent calibration 
humcalibfile = joinpath(std_fp, "fitParameters_relative.txt") #humcalibfp

#file to be calibrated at once with same mass list
resultfp = joinpath(dir_CLOUD18, "CLOUD18_data", "Nonanal", "Nonanal_10deg") #change result filepath to data that is analyzed #results of this script are also saved here
resultfiles = ["$(resultfp)/results/_result.hdf5"] #adjust filename, can add multiple files #["$(resultfp)part1/results/_result.hdf5","$(resultfp)part2/results/_result.hdf5"]

#hum dep calibration data as result file for testing if 1ppb hexanone is correctly calibrated to 1ppb at different AHs
#resultfp = std_fp 
#resultfiles = ["$(resultfp)/_result.hdf5"] 

#dry calibs data as result for testing if 1ppb hexanone is correctly calibrated to 1ppb at dry conditions
resultfp = joinpath(dir_calib_data, "dry_std", "results")
resultfiles = [drycalibsfile]

onlyUseAverages = false # false for high-time-resolution

ionization = "NH3H+" #  "H+"...
primaryionslist = [] #chosen by user
refCompound = MasslistFunctions.createCompound(C=6, H=15,O=1,N=1) # Hexanone
refMass = refCompound[1] #mass of refCompound
refName = TOFTracer2.MasslistFunctions.sumFormulaStringFromCompositionArray(refCompound[4]; ion = ionization) #name of refCompound with ionization for export

exportTraces = true # if true, check HeaderForExportDict below:
run = "theRunYouAnalyze_partX"
campaign = "CLOUD18" #for export
HeaderForExportDict = Dict(
        "title"=>"Humidity-dependent calibrated data (for directly calibrated species) of Experiment ... from $(campaign) campaign", # add here specific infos! RunNrs or Topic...
        "level"=>1, # level 0: raw data, level 1: postprocessed but unchecked data, level 2: "ready for use in publications"
        "version"=>"01", # change, if you reprocess after having the data uploaded to cernbox before
        "authorname_mail"=>"LastName, FirstName email@(student.)uibk.ac.at", 
        "units"=>"ppt",
        "addcomment"=>"The data have been humidity-dependently calibrated with Hexanone as reference (Onr=1), \n compounds with Onr>1 are calibrated with the maximum sensitivity (dry Hexanone). \n All traces have been corrected to the duty-cycle-corrected primary ion trace.\n",
        "threshold"=>0, #sigma-Threshold you apply for filtering data. If you do not apply filtering yet, keep it zero.
        "nrrows_addcomment"=>3 # Ensure, this nr corresponds to the nr of '\n' in your addcomment.
        )

# select time window and respective chamber-temperature for inlet loss correction. 
# If only one chambertemperature condition applies to the calibrated data, only one time window can be selected. 
# If multiple chamber temperature conditions apply, multiple time windows can be selected and respective chamber temperatures can be assigned, e.g. for IEPOX runs:
#=
timeranges=[(DateTime(2025,10,29), DateTime(2025,11,6,11,0)),
           (DateTime(2025,11,6,14,0), DateTime(2025,11,8,5,0))] #select time window for inlet loss correction (e.g., for nonanal run with known conditions)
chamberTs = [15,-30] # use a vector even for just one entry
=#
timeranges=[(DateTime(2025,9,1), DateTime(2025,12,2))] # select time window for inlet loss correction to correspond to times of different chamber temperatures within your calibrated dataset (example: basically no filtering (all of CLOUD18 included))
chamberTs = [15] # adjust according to your data. Use a vector, even for only 1 value!

#########################
# Run calibration steps
#########################

# 1. Run humidity dependence calibration
let
        println("Do you want to (re)run the humidity dependence calibration? (y/n)")
        input_humcalib = readline()
        while !(input_humcalib in ["y", "n"])
                println("Invalid input. Please enter 'y' for yes or 'n' for no.")
                input_humcalib = readline()
        end
        if input_humcalib == "y"
        CLOUD18_calibration.HumidityDependenceCalibration.run_humidity_dependence_calibration(std_fp, std_file, bg_file, hum_file, ion=ionization)
        end
end

# 2. Run calibrate_traces
let
        println("Do you want to apply the humidity-dependent calibration (recommended for T>0°C)? (y/n)")
        userinput = readline()
        while !(userinput in ["y", "n"])
                println("Invalid input. Please enter 'y' for yes or 'n' for no.")
                userinput = readline()
        end
        if userinput == "y"
                mResfinal, dcps_per_ppb = CLOUD18_calibration.CalibrateTraces.calibrate_traces_main(dir_licor_data, humcalibfile, drycalibsfile, resultfp, resultfiles, ionization, refMass, refName, exportTraces, HeaderForExportDict;useAverages=onlyUseAverages)
        elseif userinput == "n"
                mResfinal, dcps_per_ppb = CLOUD18_calibration.CalibrateTraces.calibrate_traces_main(drycalibsfile, resultfp, resultfiles, ionization, refMass, refName, exportTraces, HeaderForExportDict;useAverages=onlyUseAverages)
        end
end
# 3. Run Inlet Loss Correction for each chamber temperature condition applicable to the calibrated data
#CLOUD18_calibration.InletLossCorrection.run_inlet_loss_correction(resultfp; timerange=timerange, export_fp=resultfp)
CLOUD18_calibration.InletLossCorrection.run_inlet_loss_correction(resultfp; 
                timeranges=timeranges, 
                export_fp=resultfp, 
                ion = ionization, 
                flow=8.2, # 7 slpm inletflow (logbook CLOUD18) + 1.2 slpm sampleflow
                sampleflow = 1.2,
                inletLength = 0.7,
                chamberTs=chamberTs, 
                HeaderForExportDict = HeaderForExportDict, 
                campaign=campaign,
                run=run)
