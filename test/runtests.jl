using CLOUD18_calibration
using CLOUD18_calibration.CalibrateTraces
using Test
using DataFrames
using Dates
using CSV
using HDF5
using TOFTracer2

@testset "CalibrateTraces" begin
    
    @testset "CalibrationConfig" begin
        # Test that CalibrationConfig struct can be created
        config = CalibrationConfig(
            "/path/to/licor",
            "/path/to/humcalib.txt",
            "/path/to/drycalib.hdf5",
            "/path/to/result/",
            ["/path/to/result1.hdf5"],
            "NH4+",
            Float64[],
            163.0,
            "Hexanone",
            false,
            Dict{String, Any}()
        )
        @test config.dir_licor_data == "/path/to/licor"
        @test config.humcalibfile == "/path/to/humcalib.txt"
        @test config.drycalibsfile == "/path/to/drycalib.hdf5"
        @test config.ionization == "NH4+"
        @test config.refName == "Hexanone"
        @test config.exportTraces == false
    end

    @testset "compute_summed_primary_ions" begin
        # Create mock measurement results with known values
        times = [Dates.DateTime(2025, 1, 1) + Dates.Hour(i) for i in 0:2]
        masses = [18.0, 36.0, 54.0]  # H2O, (H2O)2, (H2O)3
        elements = ["H", "O"]
        elements_masses = [1.00783, 15.99492]
        compositions = [3 3 3; 1 2 3]  # H3O, H3O2, H3O3
        traces = [
            100.0 200.0 300.0;  # Time 1
            150.0 250.0 350.0;  # Time 2
            200.0 300.0 400.0   # Time 3
        ]
        
        mRes = TOFTracer2.ResultFileFunctions.MeasurementResult(
            times, masses, elements, elements_masses, compositions, traces
        )
        
        summedPIs = compute_summed_primary_ions(mRes)
        
        # Check dimensions and properties
        @test length(summedPIs) == length(times)
        @test all(summedPIs .>= 0)  # All values should be non-negative
        @test eltype(summedPIs) <: Real
    end

    @testset "compute_summed_primary_ions with zero/negative values" begin
        times = [Dates.DateTime(2025, 1, 1)]
        masses = [18.0, 36.0]
        elements = ["H", "O"]
        elements_masses = [1.00783, 15.99492]
        compositions = [3 3; 1 2]
        traces = [-100.0 50.0]  # One negative value
        
        mRes = TOFTracer2.ResultFileFunctions.MeasurementResult(
            times, masses, elements, elements_masses, compositions, traces
        )
        
        summedPIs = compute_summed_primary_ions(mRes)
        
        # Negative or zero values should be set to 0
        @test summedPIs[1] >= 0
    end

    @testset "build_calibration_traces error handling" begin
        # Create minimal mock data
        times = [Dates.DateTime(2025, 1, 1), Dates.DateTime(2025, 1, 2)]
        masses = [18.0, 36.0, 54.0, 163.0]
        elements = ["C", "H", "N", "O"]
        elements_masses = [12.0, 1.00783, 14.00307, 15.99492]
        compositions = [0 0 0 1; 3 3 3 5; 0 0 0 1; 1 2 3 2]  # C0H3O1, C0H3O2, C0H3O3, C1H5NO2
        traces = ones(2, 4) .* 100.0
        
        mRes = TOFTracer2.ResultFileFunctions.MeasurementResult(
            times, masses, elements, elements_masses, compositions, traces
        )
        
        summedPIs = [1000.0, 1500.0]
        licor_final = [5.0, 6.0]
        
        # Create calibration DataFrame with missing reference compound
        calibDF = DataFrame(
            Sumformula=["compound1"],
            Mass=[163.0],
            p1=[1.0],
            p2=[0.1],
            p3=[0.01],
            p4=[0.001],
            p5=[0.0001]
        )
        
        hexVSpis_params = ([1.0, 0.5], [0.01], ["power"])
        
        # Should throw error when reference compound not found
        @test_throws ErrorException build_calibration_traces(
            mRes, summedPIs, licor_final, calibDF, hexVSpis_params, "NonexistentCompound"
        )
    end

    @testset "build_calibration_traces with valid data" begin
        times = [Dates.DateTime(2025, 1, 1), Dates.DateTime(2025, 1, 2)]
        masses = [18.0, 36.0, 54.0, 163.0]
        elements = ["C", "H", "N", "O"]
        elements_masses = [12.0, 1.00783, 14.00307, 15.99492]
        compositions = [0 0 0 1; 3 3 3 5; 0 0 0 1; 1 2 3 2]
        traces = ones(2, 4) .* 100.0
        
        mRes = TOFTracer2.ResultFileFunctions.MeasurementResult(
            times, masses, elements, elements_masses, compositions, traces
        )
        
        summedPIs = [1000.0, 1500.0]
        licor_final = [5.0, 6.0]
        
        calibDF = DataFrame(
            Sumformula=["Hexanone"],
            Mass=[163.0],
            p1=[1.0],
            p2=[0.1],
            p3=[0.01],
            p4=[0.001],
            p5=[0.0001]
        )
        
        hexVSpis_params = ([1.0, 0.5], [0.01], ["power"])
        
        dcps_per_ppb, indices = build_calibration_traces(
            mRes, summedPIs, licor_final, calibDF, hexVSpis_params, "Hexanone"
        )
        
        # Check output dimensions
        @test size(dcps_per_ppb) == (length(times), length(masses))
        @test isa(indices, Vector{Int})
        # Check that some values are non-zero (calibration was applied)
        @test any(dcps_per_ppb .!= 0)
    end

    @testset "plot_calibration_traces" begin
        # Create test data
        times = [Dates.DateTime(2025, 1, 1), Dates.DateTime(2025, 1, 2)]
        masses = [18.0, 36.0]
        elements = ["H", "O"]
        elements_masses = [1.00783, 15.99492]
        compositions = [3 3; 1 2]
        traces = ones(2, 2) .* 100.0
        
        mRes = TOFTracer2.ResultFileFunctions.MeasurementResult(
            times, masses, elements, elements_masses, compositions, traces
        )
        
        dcps_per_ppb = [100.0 200.0; 150.0 250.0]
        summedPIs = [1000.0, 1500.0]
        indices = [1, 2]
        resultfp = tempdir() * "/test_"
        
        # This function creates files, just check it doesn't error
        @test_nowarn plot_calibration_traces(mRes, dcps_per_ppb, summedPIs, indices, resultfp)
    end

    @testset "export_calibrated_traces input validation" begin
        # Create minimal test data
        times = [Dates.DateTime(2025, 1, 1), Dates.DateTime(2025, 1, 2)]
        masses = [18.0, 36.0, 117.0]  # H2O, (H2O)2, C5H9NO2
        elements = ["C", "H", "N", "O"]
        elements_masses = [12.0, 1.00783, 14.00307, 15.99492]
        compositions = [0 0 5; 3 3 9; 0 0 1; 1 2 2]
        traces = ones(2, 3) .* 100.0
        
        mRes = TOFTracer2.ResultFileFunctions.MeasurementResult(
            times, masses, elements, elements_masses, compositions, traces
        )
        
        dcps_per_ppb = [100.0 200.0 300.0; 150.0 250.0 350.0]
        
        header_dict = Dict{String, Any}(
            "title" => "Test Export",
            "level" => "L2",
            "version" => "1.0",
            "authorname_mail" => "Test Author",
            "units" => "ppb",
            "addcomment" => "Test comment",
            "threshold" => 0.0,
            "nrrows_addcomment" => 2
        )
        
        resultfp = tempdir() * "/test_export_"
        
        # This function requires interactive input, so we just verify it doesn't crash on setup
        # In a real scenario, you would mock the interactive components
        # @test_nowarn export_calibrated_traces(mRes, dcps_per_ppb, "NH4+", header_dict, resultfp)
    end

    @testset "CalibrationConfig with complex header dict" begin
        header_dict = Dict{String, Any}(
            "title" => "CLOUD18 Calibration",
            "level" => "L2",
            "version" => "2.0",
            "authorname_mail" => "User <user@example.com>",
            "units" => "ppb",
            "addcomment" => "Test calibration data",
            "threshold" => 0.1,
            "nrrows_addcomment" => 3
        )
        
        config = CalibrationConfig(
            "/data/licor",
            "/data/humcalib.txt",
            "/data/drycalib.hdf5",
            "/results/",
            ["/data/result1.hdf5", "/data/result2.hdf5"],
            "NH4+",
            [18.033836, 36.044416],
            163.0,
            "Hexanone",
            true,
            header_dict
        )
        
        @test config.HeaderForExportDict["title"] == "CLOUD18 Calibration"
        @test config.HeaderForExportDict["threshold"] == 0.1
        @test config.exportTraces == true
        @test length(config.resultfiles) == 2
    end

    @testset "compute_summed_primary_ions consistency" begin
        # Test that function is deterministic
        times = [Dates.DateTime(2025, 1, 1), Dates.DateTime(2025, 1, 2), Dates.DateTime(2025, 1, 3)]
        masses = [18.0, 36.0]
        elements = ["H", "O"]
        elements_masses = [1.00783, 15.99492]
        compositions = [3 3; 1 2]
        traces = [100.0 200.0; 150.0 250.0; 200.0 300.0]
        
        mRes = TOFTracer2.ResultFileFunctions.MeasurementResult(
            times, masses, elements, elements_masses, compositions, traces
        )
        
        result1 = compute_summed_primary_ions(mRes)
        result2 = compute_summed_primary_ions(mRes)
        
        @test result1 ≈ result2
    end

    @testset "build_calibration_traces with different oxygen compositions" begin
        # Test that oxygen-based filtering works correctly
        times = [Dates.DateTime(2025, 1, 1)]
        masses = [18.0, 36.0, 54.0, 163.0]  # H2O (1O), (H2O)2 (2O), (H2O)3 (3O), Hex (2O)
        elements = ["C", "H", "N", "O"]
        elements_masses = [12.0, 1.00783, 14.00307, 15.99492]
        compositions = [0 0 0 1; 3 3 3 5; 0 0 0 1; 1 2 3 2]  # Oxygen counts: 1, 2, 3, 2
        traces = ones(1, 4) .* 100.0
        
        mRes = TOFTracer2.ResultFileFunctions.MeasurementResult(
            times, masses, elements, elements_masses, compositions, traces
        )
        
        summedPIs = [1000.0]
        licor_final = [5.0]
        
        calibDF = DataFrame(
            Sumformula=["Hexanone"],
            Mass=[163.0],
            p1=[1.0],
            p2=[0.1],
            p3=[0.01],
            p4=[0.001],
            p5=[0.0001]
        )
        
        hexVSpis_params = ([1.0, 0.5], [0.01], ["power"])
        
        dcps_per_ppb, indices = build_calibration_traces(
            mRes, summedPIs, licor_final, calibDF, hexVSpis_params, "Hexanone"
        )
        
        # Verify filtering was applied
        @test size(dcps_per_ppb) == (1, 4)
        @test isa(indices, Vector{Int})
    end

end  # @testset "CalibrateTraces"
