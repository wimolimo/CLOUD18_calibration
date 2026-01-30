using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
using Dates
using TOFTracer2


fp = "/media/wiebke/Extreme SSD/CLOUD16/PTR3/Surfactants/data/result_CHO/"
fpcompositions_CHO = "$(fp)ptr3compositions_CHOproducts_NH4+.txt"
fptraces_CHO = "$(fp)ptr3traces_CHOproducts_NH4+.csv"



mResult_CHON = TOFTracer2.ImportFunctions.importExportedTraces(fptraces_CHO, fpcompositions_CHO)
mResult_CHO = TOFTracer2.ImportFunctions.importExportedTraces(fptraces_CHON, fpcompositions_CHON)
mResult = ResultFileFunctions.joinResultsMasses(mResult_CHON,mResult_CHO)

# start until 26.10.2023 02:33: 0.77slpm+3.3 = 4 slpm, chamberT = -15°C
filter1 = mResult.Times .< DateTime(2023,10,26,2,33)
transmissions1 = TOFTracer2.CalibrationFunctions.calculateInletTransmission_CLOUD(mResult.MasslistMasses, mResult.MasslistCompositions; 
    ion = "NH4+", flow=4, sampleflow = 1.2,inletLength = 0.7, chamberT=-15, roomT=25, ptrT=37)
											
# 26.10.2023 02:33 until 30.10.23 19:00:  6.77slpm+3.3 = 10 slpm, chamberT = -15°C							
filter2 = DateTime(2023,10,26,2,33) .< mResult.Times .< DateTime(2023,10,30,19)
transmissions2 = TOFTracer2.CalibrationFunctions.calculateInletTransmission_CLOUD(mResult.MasslistMasses, mResult.MasslistCompositions; 
    ion = "NH4+", flow=10, sampleflow = 1.2,inletLength = 0.7, chamberT=-15, roomT=25, ptrT=37)

# T-transition
filter2a = DateTime(2023,10,30,19) .< mResult.Times .< DateTime(2023,10,30,22)

# 30.10.23 22:00 until end:  6.77slpm+3.3 = 10 slpm, chamberT = +10°C											
filter3 = DateTime(2023,10,30,22) .< mResult.Times
transmissions3 = TOFTracer2.CalibrationFunctions.calculateInletTransmission_CLOUD(mResult.MasslistMasses, mResult.MasslistCompositions; 
    ion = "NH4+", flow=10, sampleflow = 1.2,inletLength = 0.7, chamberT=10, roomT=25, ptrT=37)

# Correction		
mResult.Traces[filter1,:] .= mResult.Traces[filter1,:] ./ transpose(transmissions1)
mResult.Traces[filter2,:] .= mResult.Traces[filter2,:] ./ transpose(transmissions2)
mResult.Traces[filter2a,:] .= NaN
mResult.Traces[filter3,:] .= mResult.Traces[filter3,:] ./ transpose(transmissions3)


# exporting
# run for per filter & transmission

HeaderForExportDict = Dict(
        "title"=>"oxidized hydrocarbons from Nonanal runs at -15°C",
        "level"=>2,
        "version"=>"01",
        "authorname_mail"=>"Scholz, Wiebke wiebke.scholz@uibk.ac.at",
        "units"=>"ppt",
        "addcomment"=>"The data have been humidity-depently calibrated with Hexanone as reference (Onr=[1,2]),
        compounds with Onr>2 are calibrated with kinetic limit.
        All traces have been corrected to the duty-cycle-corrected primary ion trace.
        Uncertainty roughly factor 3. Transmission-corrected for a total flow of 4 slpm.\n",
        "threshold"=>0,
        "nrrows_addcomment" => 5
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

TOFTracer2.ExportFunctions.exportTracesCSV_CLOUD("/media/wiebke/Extreme SSD/CLOUD16/PTR3/Surfactants/data/",
        mResult.MasslistElements,
        mResult.MasslistMasses,
        mResult.MasslistCompositions,
        mResult.Times[filter1],
        mResult.Traces[filter1,:];
        transmission=transmissions1,
        headers = HeaderForExport,
        ion = "NH4+",
        average=0)
        