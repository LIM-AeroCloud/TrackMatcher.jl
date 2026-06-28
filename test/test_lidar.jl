## Tests about lidar data processing and feature classification

# Define test data and expected results
hfile = normpath(@__DIR__, "..", "data", "CPro_Lidar_Altitudes_m.dat")
h5file = normpath(@__DIR__, "data", "caliop", "cpro", "CPro_23.h5")
hprofile = CSV.read(hfile, DataFrame, copycols = false, types=Float32)
lidarprofile = TrackMatcher.get_lidarheights((Inf, -Inf), Float32)

fcf = [0x0001, 0x0019, 0x041b, 0x061b, 0x0a1b, 0x021a, 0x081a, 0x0a1a, 0x0c1a, 0x0e1a]
feature = [clear, clear, dust, polluted, polluted_dust, low_opaque, ac, as, ci, cb]

## Test sets
@testset "lidar altitude profile" begin
    @test lidarprofile.coarse == hprofile.CPro
    @test length(lidarprofile.coarse) == 399
    @test length(lidarprofile.fine) == 545
    @test lidarprofile.i.f.h30m == 253
    @test lidarprofile.i.top == 1
    @test lidarprofile.i.bottom == 399
    @test lidarprofile.fine[end] ≈ -4.41219f2
    @test typeof(lidarprofile.coarse) == typeof(lidarprofile.fine) == Vector{Float32}
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

@testset "feature classification" begin
    for i in eachindex(fcf)
        @test TrackMatcher.feature_classification(TrackMatcher.classification(fcf[i])...) ==
            feature[i]
    end
end
