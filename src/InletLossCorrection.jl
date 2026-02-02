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
function run_inlet_loss_correction(fp; timerange::Vector{DateTime} = [DateTime(2000,1,1), DateTime(3000,1,1)], export_fp::String = fp)

    fpcompositions = "$(fp)/ptr3compositions_CLOUDheader.txt"
    fptraces = "$(fp)/ptr3traces_CLOUDheader.csv"

    mResult = TOFTracer2.ImportFunctions.importExportedTraces(fptraces, fpcompositions)

    # Nonanal run (24.11.2025 05:00 - 24.11.2025 10:00): flow = 7 slpm, chamberT = 8°C
    filter1 = timerange[1] .< mResult.Times .< timerange[2]
    transmissions1 = TOFTracer2.CalibrationFunctions.calculateInletTransmission_CLOUD(mResult.MasslistCompositions; 
        ion = "NH4+", flow=7, sampleflow = 1,inletLength = 0.7, chamberT=8, roomT=25, ptrT=37)      # change temperature according to cloud data
											
    #  Correction		
    mResult.Traces[filter1,:] .= mResult.Traces[filter1,:] ./ transpose(transmissions1)

    # exporting
    # run for per filter & transmission

    HeaderForExportDict = Dict(
            "title"=>"oxidized hydrocarbons from Nonanal runs at 8°C",
            "level"=>2,
            "version"=>"01",
            "authorname_mail"=>"Ruth, Clea clea.ruth@uibk.ac.at; Wittler, Timo, wittler.timo@uibk.ac.at",
            "units"=>"ppt",
            "addcomment"=>"The data have been humidity-depently calibrated with Hexanone as reference (Onr=[1,2]),
            compounds with Onr>2 are calibrated with kinetic limit.
            All traces have been corrected to the duty-cycle-corrected primary ion trace.
            Uncertainty roughly factor 3. Transmission-corrected for a total flow of 8 slpm.\n",
            "threshold"=>0,
            "nrrows_addcomment" => 4
            )

    HeaderForExport = TOFTracer2.ExportFunctions.CLOUDheader(mResult.Times[filter1];
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
            mResult.Times[filter1],
            mResult.Traces[filter1,:];
            transmission=transmissions1,
            headers = HeaderForExport,
            ion = "NH4+",
            average=0)
end

end # module InletLossCorrection

# run_inlet_loss_correction(resultfp; timerange=[DateTime(2025,11,24,5,0), DateTime(2025,11,24,10,0)], export_fp=export_fp)