## Tests about lidar data processing and feature classification

# Define test data and expected results
hfile = normpath(@__DIR__, "..", "data", "CPro_Lidar_Altitudes_m.dat")
h5file = normpath(@__DIR__, "data", "caliop", "CPro", "2012", "CPro_23.h5")
hprofile = CSV.read(hfile, DataFrame, copycols = false, types=Float32)
lidarprofile = TrackMatcher.get_lidarheights((Inf, -Inf), Float32)

fcf = [0x0001, 0x0019, 0x041b, 0x061b, 0x0a1b, 0x021a, 0x081a, 0x0a1a, 0x0c1a, 0x0e1a]
feature = [clear, clear, dust, polluted, polluted_dust, low_opaque, ac, as, ci, cb]

## Test sets
@testset "lidar altitude profile" begin
    @testset "example test" begin
        @test lidarprofile.coarse == hprofile.CPro
        @test length(lidarprofile.coarse) == 399
        @test length(lidarprofile.fine) == 545
        @test lidarprofile.i.f.h30m == 253
        @test lidarprofile.i.top == 1
        @test lidarprofile.i.bottom == 399
        @test lidarprofile.fine[end] ≈ -4.41219f2
        @test typeof(lidarprofile.coarse) == typeof(lidarprofile.fine) == Vector{Float32}
    end
    @testset "error handling" begin
        h5open(h5file, "r") do fid
            @test_throws("selected values must overlap with",
                TrackMatcher.get_lidarheights((-1000, -Inf), Float32))
            @test_throws("selected values must overlap with",
                TrackMatcher.get_lidarheights((Inf, 30_000), Float32))
            @test_throws "altitude bounds are inverted" TrackMatcher.get_lidarheights((0, 20_000), Float32)
            # Test special cases currently not used in TrackMatcher
            @test_throws "invalid input: 2D array provided for fine resolution heights; " *
                "set coarse to true to resolve" TrackMatcher.get_lidarcolumn(Float32,
                read(fid, "Extinction_Coefficient_532"), lidarprofile, coarse=false)
            @test_nowarn TrackMatcher.get_lidarcolumn(Float32, read(fid, "Extinction_Coefficient_532"),
                lidarprofile)
            @test_nowarn TrackMatcher.get_lidarcolumn(UInt16, read(fid, "Atmospheric_Volume_Description"),
                lidarprofile, coarse=true)
            @test_nowarn TrackMatcher.get_lidarcolumn(UInt16, read(fid, "Atmospheric_Volume_Description"),
                lidarprofile, coarse=true, missingvalues = 9999)
        end
    end
    @testset "fine grid special cases" begin
        @testset "lidar range completely above fine grid" begin
            lp_fine = @test_nowarn TrackMatcher.get_lidarheights((7000, 130), Float32)
            @test lp_fine.i.f.h30m == 1
            @test length(lp_fine.fine) == 231
            @test length(lp_fine.coarse) == 117
            @test lp_fine.fine[1] ≈ 7.0133823f3
            @test lp_fine.fine[end] ≈ 1.276055f2
            @test lp_fine.coarse[1] ≈ 7.043321f3
            @test lp_fine.coarse[end] ≈ 9.7667f1
            @test lp_fine.i.top == 274
            @test lp_fine.i.bottom == 390
            @test lp_fine.i.f.top == 2
            @test lp_fine.i.f.bottom == 2
        end
        @testset "lidar range completely within fine grid" begin
            lp_coarse = @test_nowarn TrackMatcher.get_lidarheights((25_000, 10_000), Float32)
            @test lp_coarse.i.f.h30m == 200
            @test length(lp_coarse.fine) == length(lp_coarse.coarse) == 199
            @test lp_coarse.coarse == lp_coarse.fine
            @test lp_coarse.fine[1] ≈ 2.5125969f4
            @test lp_coarse.fine[end] ≈ 9.977261f3
            @test lp_coarse.i.top == 27
            @test lp_coarse.i.bottom == 225
            @test lp_coarse.i.f.top == 1
            @test lp_coarse.i.f.bottom == 0
        end
    end
end

@testset "feature classification" begin
    for i in eachindex(fcf)
        @test TrackMatcher.feature_classification(TrackMatcher.classification(fcf[i])...) ==
            feature[i]
    end
end
