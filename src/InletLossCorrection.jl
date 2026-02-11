"""
    Module InletLossCorrection

Provides functionality to correct trace gas concentrations for diffusional losses in the 
sampling inlet based on the CLOUD18 experimental setup.
"""
module InletLossCorrection
export run_inlet_loss_correction

using HDF5, DataFrames, CSV
using PyPlot
using Dates
using TOFTracer2


"""
    run_inlet_loss_correction(fp; timerange, export_fp)

Performs inlet loss correction on exported CLOUD traces and compositions.

The correction accounts for diffusional wall losses in the sampling line. It assumes a 
total flow of 7 slpm and standard CLOUD sampling parameters (0.7m inlet length). 
The function filters data by the provided `timerange`, applies the transmission 
efficiency factor, and exports the result with a standardized CLOUD metadata header.

# Arguments
- `fp::String`: Directory path containing the `ptr3compositions_CLOUDheader.txt` and `ptr3traces_CLOUDheader.csv` files.
- `timerange::Vector{DateTime}`: A two-element vector `[start, end]` defining the period to correct. 
  Defaults to the full range.
- `export_fp::String`: Destination folder path for the corrected CSV files. Defaults to `fp`.

# Scientific Parameters (Specific for CLOUD18 Inlet)
- `flow`: 7 slpm (total flow).
- `sampleflow`: 1 slpm.
- `inletLength`: 0.7 m.
- `chamberT`: 8 °C.
- `ptrT`: 37 °C.

# Returns
- No return value. The corrected dataset is saved as a new CSV file in the `export_fp` directory.
"""
function run_inlet_loss_correction(fp; timeranges:: Vector{Tuple{DateTime, DateTime}} = [(DateTime(2000,1,1), DateTime(3000,1,1))], export_fp::String = fp, ion = "NH3H+", flow=8, sampleflow = 1.2,inletLength = 0.9, chamberTs=[0], HeaderForExportDict = Dict(), campaign="CLOUD18", run::String="")

    for (time, chamberT) in zip(timeranges, chamberTs)
        println("Running inlet loss correction for chamber temperature: ", chamberT, "°C")      
        fpcompositions = "$(fp)/ptr3compositions_CLOUDheader.txt"
        fptraces = "$(fp)/ptr3traces_CLOUDheader.csv"
        mResult = TOFTracer2.ImportFunctions.importExportedTraces(fptraces, fpcompositions)
        filter = time[1] .< mResult.Times .< time[2]

        transmissions = TOFTracer2.CalibrationFunctions.calculateInletTransmission_CLOUD(mResult.MasslistCompositions; 
        ion = ion, flow=flow, sampleflow = sampleflow,inletLength = inletLength, chamberT=chamberT, roomT=25, ptrT=37)      # change temperature according to cloud data
			
        #  Correction		
        mResult.Traces[filter,:] .= mResult.Traces[filter,:] ./ transpose(transmissions)

        HeaderForExportDict["addcomment"] = string(HeaderForExportDict["addcomment"],"Transmission-corrected for a total flow rate of $(flow)slpm.\n")
        HeaderForExportDict["nrrows_addcomment"] = HeaderForExportDict["nrrows_addcomment"] + 1

        HeaderForExport = TOFTracer2.ExportFunctions.CLOUDheader(mResult.Times[filter];
            title = HeaderForExportDict["title"],
            level=HeaderForExportDict["level"],
            version=HeaderForExportDict["version"],
            authorname_mail=HeaderForExportDict["authorname_mail"],
            units=HeaderForExportDict["units"],
            addcomment=HeaderForExportDict["addcomment"],
            threshold=HeaderForExportDict["threshold"],
            nrrows_addcomment = HeaderForExportDict["nrrows_addcomment"])

        TOFTracer2.ExportFunctions.exportTracesCSV_CLOUD(export_fp,
            mResult.MasslistElements,
            mResult.MasslistMasses,
            mResult.MasslistCompositions,
            mResult.Times[filter],
            mResult.Traces[filter,:];
            transmission=transmissions,
            headers = HeaderForExport,
            ion = ion,
            average=0,
            filenameAddition="_UIBK_oVOCs_$(run)_$(campaign)_inletLossCorr_T$(chamberT)C_V1")

    end
end

end # module InletLossCorrection

# run_inlet_loss_correction(resultfp; timerange=[DateTime(2025,11,24,5,0), DateTime(2025,11,24,10,0)], export_fp=export_fp)