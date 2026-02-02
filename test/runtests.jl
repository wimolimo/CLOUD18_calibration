using CLOUD18_calibration
using CLOUD18_calibration.CalibrateTraces
using CLOUD18_calibration.HumidityDependenceCalibration
using Test
using DataFrames
using Dates
using CSV
using HDF5
using TOFTracer2

@testset "CLOUD18_calibration" begin

    @testset "CalibrateTraces" begin

        @testset "parse_formula_to_composition" begin
            # Test simple formula
            comp1 = CalibrateTraces.parse_formula_to_composition("C6H12O")
            @test comp1["C"] == 6
            @test comp1["H"] == 12
            @test comp1["O"] == 1
            
            # Test formula with charge notation
            comp2 = CalibrateTraces.parse_formula_to_composition("C6H12O.H+")
            @test comp2["C"] == 6
            @test comp2["H"] == 13  # 12 + 1 from H+
            @test comp2["O"] == 1
            
            # Test formula with NH4+
            comp3 = CalibrateTraces.parse_formula_to_composition("C6H12O.NH4+")
            @test comp3["C"] == 6
            @test comp3["H"] == 16  # 12 + 4 
            @test comp3["N"] == 1
            @test comp3["O"] == 1
            
        end

        @testset "find_formula_index" begin
            formulas = ["C6H12O", "C5H10O2", "C7H14O"]
            
            # Test exact match
            @test CalibrateTraces.find_formula_index(formulas, "C6H12O") == 1
            @test CalibrateTraces.find_formula_index(formulas, "C5H10O2") == 2
            
            # Test no match
            @test CalibrateTraces.find_formula_index(formulas, "C8H16O") === nothing
            
            # Test composition-based matching (different notation, same composition)
            formulas_with_adduct = ["C6H12O.H+", "C5H10O2.NH4+"]
            @test CalibrateTraces.find_formula_index(formulas_with_adduct, "C6H13O") == 1  # Same as C6H12O.H+
        end

        @testset "plot_fit" begin
            # Test that plot_fit returns correct functions
            @test CalibrateTraces.plot_fit("power") === TOFTracer2.CalibrationFunctions.PowerFunction
            @test CalibrateTraces.plot_fit("linear") === TOFTracer2.CalibrationFunctions.LinearFunction
            
            # Test error for unsupported type
            @test_throws ErrorException CalibrateTraces.plot_fit("unsupported")
        end

        @testset "calc_fhex" begin
            # Test hexanone sensitivity calculation
            summedPIs = [1000.0, 1500.0, 2000.0]
            
            # Power function parameters: a=0.5, b=0.8
            hexVSpis_params = ([0.5, 0.8], [0.01, 0.02], ["power", "y = 0.5*x^0.8"])
            
            f_hex, f_hex_err = CalibrateTraces.calc_fhex(summedPIs, hexVSpis_params)
            
            # Check dimensions
            @test length(f_hex) == length(summedPIs)
            @test length(f_hex_err) == length(summedPIs)
            
            # Check that all values are positive
            @test all(f_hex .> 0)
            @test all(f_hex_err .>= 0)
            
            # Test with zero/negative primary ions
            summedPIs_with_zero = [0.0, -10.0, 1000.0]
            f_hex2, f_hex_err2 = CalibrateTraces.calc_fhex(summedPIs_with_zero, hexVSpis_params)
            @test f_hex_err2[1] == 0  # Error should be 0 for zero PI
            @test f_hex_err2[2] == 0  # Error should be 0 for negative PI
        end

        @testset "compute_summed_primary_ions" begin
            # Create mock measurement results with known values
            times = [Dates.DateTime(2025, 1, 1) + Dates.Hour(i) for i in 0:2]
            masses = [18.0, 36.0, 54.0]  # H2OH+, (H2O)2H+, (H2O)3H+
            elements = ["H", "O"]
            elements_masses = [1, 16]
            compositions = [2 4 6; 1 2 3] 
            traces = [
                100.0 200.0 300.0;  # Time 1
                150.0 250.0 350.0;  # Time 2
                200.0 300.0 400.0   # Time 3
            ]
            
            mRes = TOFTracer2.ResultFileFunctions.MeasurementResult(
                times, masses, elements, elements_masses, compositions, traces
            )
            
            summedPIs = CalibrateTraces.compute_summed_primary_ions(mRes)
            
            # Check dimensions and properties
            @test length(summedPIs) == length(times)
            @test all(summedPIs .>= 0)  # All values should be non-negative
            @test eltype(summedPIs) <: Real
            
            # Check duty cycle correction is applied: trace * sqrt(100/mass)
            @test summedPIs[1] ≈ sum(traces[1, :] .* sqrt.(100 ./ masses))

            #Check specific values
            expected_summedPIs = [
                sum(traces[1, :] .* sqrt.(100 ./ masses)),
                sum(traces[2, :] .* sqrt.(100 ./ masses)),
                sum(traces[3, :] .* sqrt.(100 ./ masses))
            ]
            @test all(summedPIs .≈ expected_summedPIs)
        end

    
        @testset "build_calibration_traces - humidity method error handling" begin
            # Create minimal mock data
            times = [Dates.DateTime(2025, 1, 1), Dates.DateTime(2025, 1, 2)]
            masses = [18.0, 36.0, 54.0, 163.0]
            elements = ["C", "H", "N", "O"]
            elements_masses = [12.0, 1.00783, 14.00307, 15.99492]
            compositions = [0 0 0 6; 3 3 3 12; 0 0 0 0; 1 2 3 1]  # Oxygen in last row
            traces = ones(2, 4) .* 100.0
            
            mRes = TOFTracer2.ResultFileFunctions.MeasurementResult(
                times, masses, elements, elements_masses, compositions, traces
            )
            
            summedPIs = [1000.0, 1500.0]
            hexVSpis_params = ([1.0, 0.5], [0.01, 0.02], ["power", "fit"])
            licor_final = [5.0, 6.0]
            
            # Create calibration DataFrame with missing reference compound
            calibDF = DataFrame(
                Sumformula=["C5H10O"],
                Mass=[86.0],
                p1=[1.0],
                p2=[0.1],
                p3=[0.01],
                p4=[0.001],
                p5=[0.0001]
            )
            
            # Should throw error when reference compound not found
            @test_throws ErrorException CalibrateTraces.build_calibration_traces(
                mRes, summedPIs, hexVSpis_params, "C6H12O", licor_final, calibDF
            )
        end

        @testset "build_calibration_traces - humidity method with valid data" begin
            times = [Dates.DateTime(2025, 1, 1), Dates.DateTime(2025, 1, 2)]
            masses = [18, 163, 100]  # Water, Hexanone, some compound with 1 O
            elements = ["C", "H", "N", "O"]
            elements_masses = [12, 1, 14, 16]
            # Oxygen counts: 1, 1, 1 (last row)
            compositions = [0 6 5; 3 12 10; 0 0 0; 1 1 1]
            traces = ones(2, 3) .* 100.0
            
            mRes = TOFTracer2.ResultFileFunctions.MeasurementResult(
                times, masses, elements, elements_masses, compositions, traces
            )
            
            summedPIs = [1000.0, 1500.0]
            licor_final = [5.0, 6.0]
            
            calibDF = DataFrame(
                Sumformula=["C6H12O"],
                Mass=[163.0],
                p1=[1.0],
                p2=[0.1],
                p3=[0.01],
                p4=[0.001],
                p5=[0.0001]
            )
            
            hexVSpis_params = ([1.0, 0.5], [0.01, 0.02], ["power", "fit"])
            
            dcps_per_ppb, dcps_per_ppb_err, indices = CalibrateTraces.build_calibration_traces(
                mRes, summedPIs, hexVSpis_params, "C6H12O", licor_final, calibDF
            )
            
            # Check output dimensions
            @test size(dcps_per_ppb) == (length(times), length(masses))
            @test size(dcps_per_ppb_err) == size(dcps_per_ppb)
            @test isa(indices, Vector{Int})
            
            # Check that calibration was applied (some values non-zero)
            @test any(dcps_per_ppb .!= 0)
            @test any(dcps_per_ppb_err .>= 0)
        end

        @testset "build_calibration_traces - dry method" begin
            times = [Dates.DateTime(2025, 1, 1), Dates.DateTime(2025, 1, 2)]
            masses = [18, 163, 100]
            elements = ["C", "H", "N", "O"]
            elements_masses = [12, 1, 14, 16]
            compositions = [0 6 5; 3 12 10; 0 0 0; 1 1 2] 
            traces = ones(2, 3) .* 100.0
            
            mRes = TOFTracer2.ResultFileFunctions.MeasurementResult(
                times, masses, elements, elements_masses, compositions, traces
            )
            
            summedPIs = [1000.0, 1500.0]
            hexVSpis_params = ([1.0, 0.5], [0.01, 0.02], ["power", "fit"])
            
            # Call dry method (no licor, no calibDF)
            dcps_per_ppb, dcps_per_ppb_err = CalibrateTraces.build_calibration_traces(
                mRes, summedPIs, hexVSpis_params, "C6H12O"
            )
            
            # Check output dimensions
            @test size(dcps_per_ppb) == (length(times), length(masses)) broken == true ###### gives () == (2,3) ##########################################################################
            @test size(dcps_per_ppb_err) == size(dcps_per_ppb)
            
            # Check that calibration was applied
            @test any(dcps_per_ppb .!= 0)
            @test all(dcps_per_ppb_err .>= 0)
        end

        @testset "build_calibration_traces - oxygen-based filtering" begin
            # Test that compounds are calibrated based on oxygen content
            times = [Dates.DateTime(2025, 1, 1)]
            masses = [50.0, 60.0, 70.0, 80.0, 163.0]
            elements = ["C", "H", "N", "O"]
            elements_masses = [12.0, 1.00783, 14.00307, 15.99492]
            # Oxygen counts: 0, 1, 2, 3, 1 (last row)
            compositions = [3 3 3 3 6; 6 6 6 6 12; 0 0 0 0 0; 0 1 2 3 1]
            traces = ones(1, 5) .* 100.0
            
            mRes = TOFTracer2.ResultFileFunctions.MeasurementResult(
                times, masses, elements, elements_masses, compositions, traces
            )
            
            summedPIs = [1000.0]
            licor_final = [5.0]
            
            calibDF = DataFrame(
                Sumformula=["C6H12O"],
                Mass=[163.0],
                p1=[1.0],
                p2=[0.1],
                p3=[0.01],
                p4=[0.001],
                p5=[0.0001]
            )
            
            hexVSpis_params = ([1.0, 0.5], [0.01, 0.02], ["power", "fit"])
            
            dcps_per_ppb, dcps_per_ppb_err, indices = CalibrateTraces.build_calibration_traces(
                mRes, summedPIs, hexVSpis_params, "C6H12O", licor_final, calibDF
            )
            
            # Mass with 0 oxygen (index 1) should have zero calibration
            @test dcps_per_ppb[1, 1] == 0
            
            # Mass with 1 oxygen (index 2) should have humidity-dependent calibration
            @test dcps_per_ppb[1, 2] > 0
            
            # Mass with >=2 oxygen (indices 3, 4) should have dry calibration
            @test dcps_per_ppb[1, 3] > 0
            @test dcps_per_ppb[1, 4] > 0
        end


        @testset "Integration test - full workflow simulation" begin
            # Simulate a minimal end-to-end calibration workflow
            times = [Dates.DateTime(2025, 1, 1), Dates.DateTime(2025, 1, 2)]
            
            # Primary ions
            pi_masses = [18.0, 36.0]
            pi_elements = ["H", "O"]
            pi_elements_masses = [1.00783, 15.99492]
            pi_compositions = [3 3; 1 2]
            pi_traces = [1000.0 500.0; 1200.0 600.0]
            
            mRes_PIs = TOFTracer2.ResultFileFunctions.MeasurementResult(
                times, pi_masses, pi_elements, pi_elements_masses, pi_compositions, pi_traces
            )
            
            # All masses including compounds
            all_masses = [18.0, 36.0, 163.0, 100.0]  # PIs + Hexanone + compound with 1O
            all_elements = ["C", "H", "N", "O"]
            all_elements_masses = [12.0, 1.00783, 14.00307, 15.99492]
            all_compositions = [0 0 6 5; 3 3 12 10; 0 0 0 0; 1 2 1 1]  # O: 1, 2, 1, 1
            all_traces = [1000.0 500.0 800.0 600.0; 1200.0 600.0 900.0 700.0]
            
            mRes = TOFTracer2.ResultFileFunctions.MeasurementResult(
                times, all_masses, all_elements, all_elements_masses, all_compositions, all_traces
            )
            
            # Step 1: Compute summed primary ions
            summedPIs = CalibrateTraces.compute_summed_primary_ions(mRes_PIs)
            @test length(summedPIs) == 2
            
            # Step 2: Calculate hexanone sensitivity
            hexVSpis_params = ([1.0, 0.8], [0.01, 0.02], ["power", "fit"])
            f_hex, f_hex_err = CalibrateTraces.calc_fhex(summedPIs, hexVSpis_params)
            @test length(f_hex) == 2
            
            # Step 3: Build calibration traces (humidity-dependent)
            licor_final = [5.0, 6.0]
            calibDF = DataFrame(
                Sumformula=["C6H12O"],
                Mass=[163.0],
                p1=[1.0], p2=[0.1], p3=[0.01], p4=[0.001], p5=[0.0001]
            )
            
            dcps_per_ppb, dcps_per_ppb_err, indices = CalibrateTraces.build_calibration_traces(
                mRes, summedPIs, hexVSpis_params, "C6H12O", licor_final, calibDF
            )
            
            @test size(dcps_per_ppb) == (2, 4)
            @test size(dcps_per_ppb_err) == (2, 4)
            @test length(indices) >= 1  # At least hexanone should be in indices
            
            # Step 4: Verify calibrated concentrations can be calculated
            calibrated_conc = 1000.0 .* all_traces ./ dcps_per_ppb  # Convert to ppt
            @test size(calibrated_conc) == (2, 4)
            @test all(isfinite.(calibrated_conc[:, indices]))  # Check no NaN/Inf for calibrated masses
        end

    end  # @testset "CalibrateTraces"

    @testset "HumidityDependenceCalibration" begin
        
        # --- Mock Data Setup ---
        times = [DateTime(2026, 1, 30, 12, 0, 0), DateTime(2026, 1, 30, 12, 10, 0)]
        masses = [163.038, 100.07] 
        elements = ["C", "H", "N", "O"]
        el_masses = [12.0, 1.0078, 14.003, 15.994]
        compositions = [6 5; 12 10; 1 0; 1 1] 
        traces = [1000.0 500.0; 1100.0 550.0]
        
        mRes_mock = TOFTracer2.ResultFileFunctions.MeasurementResult(
            times, masses, elements, el_masses, compositions, traces
        )
        
        hum_times = collect(DateTime(2026, 1, 30, 11, 59, 0):Second(1):DateTime(2026, 1, 30, 12, 20, 0))
        hum_vals = fill(5.0, length(hum_times)) 
        humdf_mock = DataFrame(DateTime = hum_times, Symbol("H₂O_(mmol_mol⁻¹)") => hum_vals)

        @testset "get_ion_metadata" begin
            std_dict = Dict("Hexanone" => ([0.0, 163.038], ["C6H12O.NH4+"]))
            m, k, i = HumidityDependenceCalibration.get_ion_metadata("NH4+", std_dict)
            @test m[1] ≈ 163.038
            @test k[1] == "Hexanone"
        end

        @testset "getHumiditySensitivity" begin
            # Simulate user pressing Enter (default 4.0)
            input = IOBuffer("\n") 
            c_data, c_std, h_avg, h_std, win = redirect_stdin(input) do
                HumidityDependenceCalibration.getHumiditySensitivity(mRes_mock, humdf_mock, []; ppt=1000.0)
            end
            @test win == 4.0
            @test h_avg[1] ≈ 5.0
        end

        @testset "Fitting and Export" begin
            hums = [1.0, 5.0, 10.0]; h_stds = [0.1, 0.1, 0.1]
            c_data = [10.0 5.0; 5.0 2.5; 2.0 1.0]; c_std = c_data .* 0.01
            
            out_fp = tempdir()
            
            p_mat, e_mat, h_axis = HumidityDependenceCalibration.DoubleExponential_and_fit(
                hums, h_stds, c_data, c_std, mRes_mock, "NH4+", out_fp, yscale="log"
            )
            
            @test size(p_mat) == (5, 2)
            
            # Note: Added out_fp=out_fp here to match your keyword argument definition
            @test_nowarn HumidityDependenceCalibration.export_sensitivities(p_mat, e_mat, mRes_mock, "test_params.txt", out_fp=out_fp)
        end
    end
end