## Tests about lidar data processing and feature classification

# Define test data and expected results
hfile = normpath(@__DIR__, "..", "data", "CPro_Lidar_Altitudes_m.dat")
h5file = normpath(@__DIR__, "data", "caliop", "cpro", "CPro_23.h5")
hprofile = CSV.read(hfile, DataFrame, copycols = false, types=Float32)
lidarprofile = TrackMatcher.get_lidarheights((Inf, -Inf), Float32)
lp_wrongtop = TrackMatcher.get_lidarheights((-1000, -Inf), Float32)
lp_wrongbottom = TrackMatcher.get_lidarheights((Inf, 30_000), Float32)
lp_inverted = TrackMatcher.get_lidarheights((0, 20_000), Float32)

fcf = [0x0001, 0x0019, 0x041b, 0x061b, 0x0a1b, 0x021a, 0x081a, 0x0a1a, 0x0c1a, 0x0e1a]
feature = [clear, clear, dust, polluted, polluted_dust, low_opaque, ac, as, ci, cb]

## Test sets
@testset "lidar altitude profile" begin
    @test lidarprofile.coarse == hprofile.CPro
    @test length(lidarprofile.coarse) == 399
    @test length(lidarprofile.fine) == 546
    @test lidarprofile.i30 == 253
    @test lidarprofile.itop == 1
    @test lidarprofile.ibottom == 399
    @test lidarprofile.maxtop ≈ 29796.324f0
    @test lidarprofile.maxbottom ≈ -441.219f0
    @test lidarprofile.fine[end] ≈ -471.15698f0
    @test typeof(lidarprofile.coarse) == typeof(lidarprofile.fine) == Vector{Float32}
    @test isempty(lp_wrongtop.coarse) && isempty(lp_wrongtop.fine)
    @test isempty(lp_wrongbottom.coarse) && isempty(lp_wrongbottom.fine)
    @test isempty(lp_inverted.coarse) && isempty(lp_inverted.fine)
    h5open(h5file, "r") do fid
        @test_throws "selected values must overlap with" TrackMatcher.get_lidarcolumn(
            UInt16, read(fid, "Atmospheric_Volume_Description"), lp_wrongtop, coarse=false)
        @test_throws ArgumentError TrackMatcher.get_lidarcolumn(Float32,
            read(fid, "Extinction_Coefficient_532"), lp_wrongtop, missingvalues = -9999)
        @test_throws "selected values must overlap with" TrackMatcher.get_lidarcolumn(UInt16,
            read(fid, "Atmospheric_Volume_Description"), lp_wrongbottom, coarse=false)
        @test_throws ArgumentError TrackMatcher.get_lidarcolumn(Float32,
            read(fid, "Extinction_Coefficient_532"), lp_wrongbottom, missingvalues = -9999)
        @test_throws "the altitude range may be inverted" TrackMatcher.get_lidarcolumn(
            UInt16, read(fid, "Atmospheric_Volume_Description"), lp_inverted, coarse=false)
        @test_throws ArgumentError TrackMatcher.get_lidarcolumn(Float32,
            read(fid, "Extinction_Coefficient_532"), lp_inverted, missingvalues = -9999)
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
