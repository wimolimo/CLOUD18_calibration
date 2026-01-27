using HDF5
using PyCall
using PyPlot
using Dates
using CSV
using DataFrames
using Statistics
import LsqFit
using TOFTracer2
import TOFTracer2.InterpolationFunctions as IntpF
import TOFTracer2.CalibrationFunctions as CalF
import TOFTracer2.ExportFunctions as ExpF
import TOFTracer2.ImportFunctions as ImpF

random_file = joinpath("C://Users//c7441399//Documents//Atemluft", "2026-01-22-beginn-der-aufzeichnungen.txt")
flight_fp = joinpath("C:\\Users\\c7441399\\Documents\\Atemluft\\flights")
partector_file = joinpath("C:\\Users\\c7441399\\Documents\\Atemluft", "Paterctor_TEST.csv")

plotStart = DateTime(2000, 1, 1, 0, 0, 0)
plotEnd = DateTime(3000, 1, 1, 0, 0, 0)

#######################
#Atemluft file load
##########################
    # parse ISO timestamps with optional "+HH:MM" or "Z" offset.
    # If apply_offset==true the timezone offset is APPLIED to the base time
    # (e.g. "2026-01-20T00:00:00+01:00" -> 2026-01-20T01:00:00).
    # If apply_offset==false the offset is ignored and only the base timestamp is returned.
    function parse_iso_with_offset(s::AbstractString, apply_offset::Bool=true)
        # handle trailing Z (UTC marker) -> treat as no offset
        if endswith(s, "Z")
            s2 = replace(s, "Z" => "")
            return DateTime(s2, dateformat"yyyy-mm-ddTHH:MM:SS")
        end

        m = match(r"^(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})([+-]\d{2}:\d{2})$", s)

        if m === nothing
            # fallback: strip any timezone suffix and parse
            s2 = replace(s, r"[+-]\d{2}:\d{2}$" => "")
            return DateTime(s2, dateformat"yyyy-mm-ddTHH:MM:SS")
        end
        
        base = DateTime(m.captures[1], dateformat"yyyy-mm-ddTHH:MM:SS")

        if !apply_offset
            return base
        end

        off = m.captures[2]               # e.g. "+01:00"
        sign = off[1]
        hh = parse(Int, off[2:3])
        mm = parse(Int, off[5:6])
        minutes = hh*60 + mm

        return sign == '+' ? base + Minute(minutes) : base - Minute(minutes)
    end


"""
    load_plotLicorData(humfile;ax="None", header=1)
    
loads and plots the given licor file. 

- header gives the line, in which the header is located (typically ==1 or ==2)
- if ax (PyCall.PyObject) is given, it will plot the data in that axis, else, if will create a new figure
"""
function load_plotAtemluftData(humfile, flight_filepath, partector_file; ax="None", header_co2=1, apply_partector_offset::Bool=true)

    humdat=DataFrame(CSV.File(humfile, header = header_co2))
    humtime = humdat[!,"System_Date_(Y-M-D)"] .+ humdat[!,"System_Time_(h:m:s)"]
    humdat[!,"DateTime"] = humtime

    # list CSV files in flight directory (full paths)
    flight_files = readdir(flight_filepath; join=true)
    flight_files = filter(f -> isfile(f) && endswith(lowercase(f), ".csv"), flight_files)
    sort!(flight_files)

    particle_number = CSV.read(partector_file, DataFrame; header=4, delim=',')
    particle_col = string.(particle_number[!,"dateTime"])
    

    # broadcast the boolean as a Ref so each call receives the same bool
    particle_time = parse_iso_with_offset.(particle_col, Ref(apply_partector_offset))
    particle_counts = particle_number[!,"particle_number_concentration"]

    fig = figure()
    ax2 = subplot()
    h2o_mmol = humdat[!,"CO₂_(µmol_mol⁻¹)"]

    ax2.plot(humtime,h2o_mmol, label = "CO₂", linewidth=1)
    ax3 = ax2.twinx()

    for file in flight_files

        # read CSV (ensure header row and delimiter), then take column named "actual_on"
        flightdata = CSV.read(file, DataFrame; header=1, delim=',')
        # column contains ISO-like timestamps "2026-01-23T19:20:04Z"
        flightcol = string.(flightdata[!,"actual_on"])
        # remove trailing 'Z' (UTC marker) and parse to DateTime
        flightcol = replace.(flightcol, "Z" => "")
        flighttime = DateTime.(flightcol, dateformat"yyyy-mm-ddTHH:MM:SS")
        ax2.axvline.(flighttime, color="red", linestyle="--", alpha=0.5, linewidth=1)

    end

    ax3.plot(particle_time, particle_counts, label="Partector particle number concentration", color="green", alpha=0.7, linewidth=1)
    ax3.set_ylabel("Partector particle number concentration [#/cm³]")
    ax2.set_ylabel("CO₂ [mmol mol⁻¹]")
    ax2.legend(loc=1)
    ax2.set_yscale("linear")
    return (humdat, fig, ax2)
 end

(humDat_random, fig_random, ax_random) = load_plotAtemluftData(random_file, flight_fp, partector_file; header_co2=2, apply_partector_offset=false)



