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
drycalibsfile = joinpath(dir_calib_data, "CLOUD18_dry_calibrations.h5") #FILENAME ANPASSEN!!! WIE SIEHT DIESE FILE AUS???

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


#####################################################
# load and prepare metadata of the final calibration
#####################################################
### find relation of sum of primary ions to hexanone
if !isdefined(Main,:hexanone_vs_pi_params)
    if HDF5.ishdf5(drycalibsfile)
        hexanone_vs_pi_params = CalF.dryCal_selectPIandRefDataFromIFIG(drycalibsfile)
        hexanone_vs_pi_params_to_export = vcat(hexanone_vs_pi_params[3],["parameters" "errors"],hcat(hexanone_vs_pi_params[1],hexanone_vs_pi_params[2]))
        writedlm("$(dirname(drycalibsfile))hexanone_vs_pi_params_csv.csv",hexanone_vs_pi_params_to_export)
    else #if file is not hdf5, assume csv and read data as DataFrame, and construct tuple
        a = CSV.read(drycalibsfile,DataFrame, header=[2])
        b = CSV.read(drycalibsfile,DataFrame;footerskip=3,header=false)
        hexanone_vs_pi_params = (a.parameters,a.errors,[values(b[1,:])[1],values(b[1,:])[2]])        
    end
end

#load humidity data inlet
if !isdefined(Main,:licorDat)
    licorDat = ImpF.createLicorData_fromFiles(dir_licor_data;
        filefilter=r"licor_.*\.txt", #rename file to use with licor_restofthename.txt
        headerrow=2,
        columnNameOfInterest="H₂O_(mmol_mol⁻¹)",
        type_columnOfInterest=Float64)
end


#-----------------------------------------------------------------------------------
# load calib parameters relative to Hexanone from file and plot humidity-dependence
#-----------------------------------------------------------------------------------

calibDF = CSV.read(humidcalibsfile, DataFrame; header=2)

CalF.plot_humdep_fromCalibParameters(;calibDF=calibDF,
    humparams=humparams, #from CalF.getInletCLOUDHumidityRelation #EXCHANGE
    cloudhum=cloudhum, #REALLY? GET FROM LICOR?
    hum4plot=collect(0:0.2:12), # define range of humidities to plot
    savefp=dirname(dir_humidcalib),
    humdepcalibRelationship="double exponential",
    humidityRelationship="exponential",
    ionization=ionization
    )


###############################################
# load data to create final calibration traces
###############################################

if isempty(primaryionslist)
    primaryionslist = massLibrary.FullPrimaryionslist_NH4soft
end

mResfinal_PIs = ResultFileFunctions.loadResults(resultsfiles_to_calibrate[1]; useAveragesOnly=true, massesToLoad=primaryionslist)
mResfinal = ResultFileFunctions.loadResults(resultsfiles_to_calibrate[1]; useAveragesOnly=true)
if length(resultsfiles_to_calibrate) > 1
    for i in 2:length(resultsfiles_to_calibrate)
        mres_pi_i = ResultFileFunctions.loadResults(resultsfiles_to_calibrate[i]; useAveragesOnly=true, massesToLoad=primaryionslist)
        mResfinal_PIs = ResultFileFunctions.joinResultsTime(mResfinal_PIs, mres_pi_i)
        mres_i = ResultFileFunctions.loadResults(resultsfiles_to_calibrate[i]; useAveragesOnly=true)
        mResfinal = ResultFileFunctions.joinResultsTime(mResfinal, mres_i)
    end
end

summedPIs = mResfinal_PIs.Traces*sqrt.(100 ./ mResfinal_PIs.MasslistMasses)
summedPIs[summedPIs .<= 0] .= 0







end # module Calibrate_Traces