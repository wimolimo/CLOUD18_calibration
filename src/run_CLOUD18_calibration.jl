"""
    Entry script to run the CLOUD18 trace calibration, containing all needed parameters and filepaths.
"""
using .CLOUD18_calibration

#################################
# Define parameters for calibration
#################################

"""
    struct CalibrationConfig

A struct to hold configuration parameters for trace calibration.

#Arguments
- `dir_licor_data::String`: Directory path containing Licor data files.
- `hexVSpis_params`: Tuple containing hexanone vs. primary ion parameters.
- `licorDat`: DataFrame containing Licor humidity data.
- `humcalibfile::String`: Path to humidity calibration file (CSV).
- `drycalibsfile::String`: Path to dry calibration file (HDF5 or CSV).
- `resultfp::String`: File path to save the results.
- `resultfiles::Vector{String}`: List of paths to result files that can be calibrated at the same time (HDF5).
- `ionization::String`: Ionization method used (e.g., "NH4+").
- `primaryionslist::Vector{Float64}`: List of primary ion masses to load.
- `refMass::Float64`: Mass of the reference compound.
- `refName::String`: Name of the reference compound.
- `exportTraces::Bool`: Flag to indicate whether to export calibrated traces.
- `HeaderForExportDict::Dict{String,Any}`: Dictionary containing information for export.
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

dir_CLOUD18 = joinpath(@__DIR__, "..", "..")
dir_calib_data = joinpath(dir_CLOUD18, "CLOUD18_data", "Calibration")
dir_licor_data = joinpath(dir_CLOUD18, "CLOUD18_data", "Licor")

#dry calibration file # this ican be either the processed file of the dry calibrations or the CSV file containing the exported hexanone vs primary ion parameters for loading them:
drycalibsfile = joinpath(dir_calib_data, "dry_std", "results", "_result.hdf5")

#humidity dependent calibration file, from humidity_dependence_calibration.jl
# humcalibfile = joinpath(dir_calib_data, "Humidity-dependent_std", "results", "fitParameters_relative.txt") #humcalibfp

#file to be calibrated at once with same mass list
resultfp = joinpath(dir_calib_data, "Test") #change result filepath to data that is analyzed #results of this script are also saved here
resultfiles = ["$(resultfp)/results/_result.hdf5"] #adjust filename, can add multiple files #["$(resultfp)part1/results/_result.hdf5","$(resultfp)part2/results/_result.hdf5"]

ionization = "NH4+" # "NH4+", "H+"...
primaryionslist = [] # leave empty -> default: adding all possible water and ammonium clusters

refMass = massLibrary.HEXANONE_nh4[1] #mass of hexanone + NH4+ from julia package manualMassLibrary.jl
refName = TOFTracer2.MasslistFunctions.sumFormulaStringFromCompositionArray(massLibrary.HEXANONE_nh4[4]; ion = "")

exportTraces = true # if true, check HeaderForExportDict below:
HeaderForExportDict = Dict(
        "title"=>"Exampletitle...",
        "level"=>2,
        "version"=>"01",
        "authorname_mail"=>"Ruth, Clea clea.ruth@student.uibk.ac.at",
        "units"=>"ppt",
        "addcomment"=>"The data have been humidity-depently calibrated with Hexanone as reference (Onr=1), compounds with Onr>1 are calibrated with kinetic limit. All traces have been corrected to the duty-cycle-corrected primary ion trace. Uncertainty roughly factor 3. Not transmission-corrected yet.\n",
        "threshold"=>0,
        "nrrows_addcomment"=>4
        )
      




#########################
# Run calibration steps
#########################

# 1. Run humidity dependence calibration
CLOUD18_calibration.HumidityDependence.humidity_dependence_calibration_main() #ANPASSEN
humcalibfile = humcalibfp #from humidity dependence calibration

# 2. Run calibrate_traces
config = CalibrationConfig(dir_licor_data, humcalibfile, drycalibsfile, resultfp, resultfiles, ionization, primaryionslist, refMass, refName, exportTraces, HeaderForExportDict)
CLOUD18_calibration.CalibrateTraces.calibrate_traces_main(config)

# 3. Run Inlet Loss Correction