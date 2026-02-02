@testset "conversions" begin
    # Define input data
    df_empty = DataFrame(A = String[], B = Int[], C = AbstractFloat[], D = Union{Missing, Float64}[],
        E = Vector{Float16}[], F = Vector{<:Union{Missing,Float64}}[], G = Float32[])
    df_full = DataFrame(
        A = ["a", "b", "c"],
        B = [1, 2, 3],
        C = [1.1, 2.2, 3.3],
        D = [missing, 2.5, 3.5],
        E = [[1.0f0, 2.0f0], [3.0f0], [4.0f0, 5.0f0, 6.0f0]],
        F = [[missing, 2.0], [3.0], [4.0, missing, 6.0]],
        G = Float32[1.5, 2.5, 3.5],
        H = [missing, missing, missing]
    )
    t_pre2010 = 60612.03728585648
    t_poost2010 = 190102.00791336806
    @testset "floats" begin
        @test typeof.(eachcol(TrackMatcher.convert_floats!(df_empty, Float32))) ==
            [Vector{String}, Vector{Int64}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}},
            Vector{Vector{Union{Missing, Float32}}}, Vector{Vector{Union{Missing, Float32}}}, Vector{Union{Missing, Float32}}]
        @test typeof.(eachcol(TrackMatcher.convert_floats!(df_full, Float16))) ==
            [Vector{String}, Vector{Int64}, Vector{Union{Missing, Float16}}, Vector{Union{Missing, Float16}},
            Vector{Vector{Union{Missing, Float16}}}, Vector{Vector{Union{Missing, Float16}}}, Vector{Union{Missing, Float16}}, Vector{Union{Missing, Float16}}]
    end
    @testset "UTC" begin
        @test TrackMatcher.convert_utc(t_pre2010) == DateTime(2006,6,12,0,53,41,498) broken=true
        @test TrackMatcher.convert_utc(t_poost2010) == DateTime(2019,1,2,0,11,23,715)
    end
    @testset "earth radius" begin
        @test TrackMatcher.earthradius(Float32(0.)) ≈ 6.378137e6 atol=0.5
        @test TrackMatcher.earthradius(Float16(30)) ≈ 6.372824e6 atol=0.5
        @test TrackMatcher.earthradius(Float32(30)) ≈ 6.372824e6 atol=0.5
        @test TrackMatcher.earthradius(Float64(30)) ≈ 6.372824e6 atol=0.5
        @test typeof(TrackMatcher.earthradius(Float16(30))) == Float64
        @test typeof(TrackMatcher.earthradius(Float32(30))) == Float64
        @test typeof(TrackMatcher.earthradius(Float64(30))) == Float64
        @test TrackMatcher.earthradius(Float32(45)) ≈ 6.367489e6 atol=0.5
        @test TrackMatcher.earthradius(Float32(60)) ≈ 6.362132e6 atol=0.5
        @test TrackMatcher.earthradius(Float32(90)) ≈ 6.356752e6 atol=0.5
    end
    @testset "units" begin
        @test TrackMatcher.ft2m(Float32(1000.)) ≈ Float32(304.8)
        @test TrackMatcher.ft2m(Float64(1000.)) ≈ Float64(304.8)
        @test typeof(TrackMatcher.ft2m(Float32(1000.))) == Float32
        @test typeof(TrackMatcher.ft2m(Float64(1000.))) == Float64
        @test TrackMatcher.ft2m(missing) === missing
        @test TrackMatcher.ft2m(Float32(40_000.)) ≈ 12192.0
        @test TrackMatcher.ftpmin2mps(Float32(100.)) ≈ 0.508
        @test TrackMatcher.ftpmin2mps(Float64(100.)) ≈ 0.508
        @test typeof(TrackMatcher.ftpmin2mps(Float32(100.))) == Float32
        @test typeof(TrackMatcher.ftpmin2mps(Float64(100.))) == Float64
        @test TrackMatcher.ftpmin2mps(missing) === missing
        @test TrackMatcher.ftpmin2mps(Float32(6000.)) ≈ 30.48
        @test TrackMatcher.knot2mps(Float32(10.)) ≈ 5.14444
        @test TrackMatcher.knot2mps(Float64(10.)) ≈ 5.14444 atol=5e-6
        @test typeof(TrackMatcher.knot2mps(Float32(10.))) == Float32
        @test typeof(TrackMatcher.knot2mps(Float64(10.))) == Float64
        @test TrackMatcher.knot2mps(missing) === missing
        @test TrackMatcher.knot2mps(Float32(180.)) ≈ 92.59992
    end
end
