## Define expected results

cpro_files = [joinpath("Level1", "CPro1.h5"), joinpath("Level1", "Level2", "CPro2.h5")]
clay_files = [joinpath("Level1", "CLay.h5")]
mixed_files = [joinpath("Level1", "CLay.h5"), joinpath("Level1", "CPro1.h5"),
    joinpath("Level1", "Level2", "CPro2.h5")]
empty_files = String[]
mini_cpro, mini_clay = ["CPro_01.h5", "CPro_02.h5"], ["CLay_01.h5", "CLay_02.h5"]
utc = [
    DateTime(2012, 2, 6, 0, 13, 5, 668),
    DateTime(2012, 2, 6, 0, 13, 6, 412),
    DateTime(2012, 2, 6, 0, 13, 7, 156),
    DateTime(2012, 2, 6, 0, 13, 7, 900),
    DateTime(2012, 2, 6, 0, 13, 8, 644),
    DateTime(2012, 2, 6, 0, 13, 9, 388),
    DateTime(2012, 2, 6, 0, 13, 10, 132),
    DateTime(2012, 2, 6, 0, 13, 10, 876),
    DateTime(2012, 2, 6, 0, 13, 11, 620),
    DateTime(2012, 2, 6, 0, 13, 12, 364),
    DateTime(2012, 2, 6, 0, 13, 13, 108),
    DateTime(2012, 2, 6, 0, 13, 13, 852),
    DateTime(2012, 2, 6, 0, 13, 14, 596),
    DateTime(2012, 2, 6, 0, 13, 15, 340),
    DateTime(2012, 2, 6, 0, 13, 16, 084),
    DateTime(2012, 2, 6, 0, 13, 16, 828),
    DateTime(2012, 2, 6, 0, 13, 17, 572),
    DateTime(2012, 2, 6, 0, 13, 18, 316),
    DateTime(2012, 2, 6, 0, 13, 19, 060),
    DateTime(2012, 2, 6, 0, 13, 19, 804),
    DateTime(2012, 2, 6, 0, 13, 20, 548),
    DateTime(2012, 2, 6, 0, 13, 21, 292),
    DateTime(2012, 2, 6, 0, 13, 22, 036),
    DateTime(2012, 2, 6, 0, 13, 22, 780),
    DateTime(2012, 2, 6, 0, 13, 23, 524),
    DateTime(2012, 2, 6, 0, 13, 24, 267),
    DateTime(2012, 2, 6, 0, 13, 25, 011),
    DateTime(2012, 2, 6, 0, 13, 25, 755),
    DateTime(2012, 2, 6, 0, 13, 26, 499),
    DateTime(2012, 2, 6, 0, 13, 27, 243),
    DateTime(2012, 2, 6, 0, 13, 27, 987)
]
lat = Float32[6.2779856, 6.233021, 6.188323, 6.1434226, 6.0986094, 6.053796, 6.0087543, 5.9637613,
    5.918778, 5.8737893, 5.8287735, 5.7838264, 5.7388263, 5.693794, 5.648815, 5.604101, 5.559198,
    5.514163, 5.4691825, 5.424208, 5.3792043, 5.3341794, 5.289179, 5.2441573, 5.1991415, 5.1542006,
    5.109221, 5.064237, 5.019097, 4.973801, 4.928818]
lon = Float32[23.687, 23.67742, 23.667685, 23.657978, 23.648378, 23.638838, 23.629349, 23.619642,
    23.609898, 23.600185, 23.59066, 23.580997, 23.571411, 23.561863, 23.552141, 23.54242, 23.53282,
    23.523283, 23.513565, 23.504055, 23.494406, 23.484657, 23.47527, 23.465548, 23.4559, 23.446463,
    23.436754, 23.427029, 23.417404, 23.40788, 23.398176]


## Helper functions to evaluate tests

"""
    test_sat_datafiles(files, type, satfiles, expected_type=type) -> Bool

Test function `TrackMatcher.sat_datafiles` to return the `expected_type` from the given `type`
and clean `files` to match `satfiles`.
Returns `true` if both tests pass, `false` otherwise.
"""
function test_sat_datafiles(files, type, satfiles, expected_type=type)::Bool
    success = false
    input_files = copy(files)
    sattype = TrackMatcher.sat_datafiles!(input_files, type)
    success |= sattype == expected_type
    success |= input_files == satfiles
    return success
end

missing_aware_approx(x, y; atol=0.0, rtol=sqrt(eps(Float64))) =
    ismissing(x) ? ismissing(y) : (!ismissing(y) && isapprox(x, y; atol, rtol))

nested_missing_aware_approx(a, b; atol=0.0, rtol=sqrt(eps(Float64))) =
    length(a) == length(b) && all(zip(a, b)) do (u, v)
        length(u) == length(v) && all(missing_aware_approx.(u, v; atol, rtol))
    end


## Testsets

@testset "read satellite data" begin
    @test TrackMatcher.scandir(joinpath("data", "caliop", "correct"), ".h5") == mixed_files
    @test TrackMatcher.scandir(joinpath("data", "caliop", "correct"), [".h5"]) == mixed_files
    @test TrackMatcher.scandir(joinpath("data", "caliop", "correct"), [".h5", ".hdf"]) == mixed_files
    @test isempty(TrackMatcher.scandir(joinpath("data", "caliop", "no_satdata"), ".h5"))
    @test @test_logs(
        (:info, "satellite data type auto-detected as CPro"),
        (:warn, "the given folder contains 1 non-CPro files which will be ignored"),
        test_sat_datafiles(mixed_files, :undef, cpro_files, :CPro)
    )
    @test @test_logs(
        (:warn, "the given folder contains 1 non-CPro files which will be ignored"),
        test_sat_datafiles(mixed_files, :CPro, cpro_files)
    )
    @test @test_logs(
        (:warn, "the given folder contains 2 non-CLay files which will be ignored"),
        test_sat_datafiles(mixed_files, :CLay, clay_files)
    )
    @test test_sat_datafiles(cpro_files, :CPro, cpro_files)
    @test @test_logs(
        (:warn, "the given folder contains 2 non-CLay files which will be ignored"),
        (:warn, "no CLay data files found, TrackMatcher requires CALIOP data"),
        test_sat_datafiles(cpro_files, :CLay, empty_files)
    )
    @test @test_logs(
        (:warn, "no satellite data files found, TrackMatcher requires CALIOP data"),
        test_sat_datafiles(empty_files, :undef, empty_files, :undef)
    )
end

@testset "SatSet" begin
    @testset "CPro" begin
        mktempdir() do root
            src = joinpath(@__DIR__, "data", "caliop", "cpro")
            cp.(joinpath.(src, mini_cpro), joinpath.(root, mini_cpro))
            sat = SatSet(root, type=:CPro)
            @test length(sat.granules) == 2
            @test sat.granules.lat isa Vector{Vector{Float32}}
            @test sat.granules.lon isa Vector{Vector{Float32}}
            @test sat.metadata.roots[0x0001] == realpath(root)
            @test sat.metadata.date.start == DateTime(2012, 2, 5, 5, 40, 3, 659)
            @test sat.metadata.date.stop == DateTime(2012, 2, 5, 7, 18, 50, 827)
            @test all(sat.metadata.granules.latmin .≈ [-69.07832, -81.821655])
            @test all(sat.metadata.granules.latmax .≈ [81.82099, 78.0213])
            @test all(sat.metadata.granules.elonmin .≈ [0.013624359, 53.707783])
            @test all(sat.metadata.granules.elonmax .≈ [77.15055, 179.803])
            @test all(sat.metadata.granules.wlonmin .≈ [-92.6476, -179.90944])
            @test all(sat.metadata.granules.wlonmax .≈ [-0.22770761, -92.70058])
            @test sat.metadata.type == :CPro
        end
        mktempdir() do root
            touch(joinpath(root, "CPro_01.h5"))
            sat = @test_logs(
                (:warn, "read error; data skipped"),
                (:warn, "no satellite data files successfully loaded"),
                SatSet(root, type=:CPro)
            )
            @test sat.metadata.date == (start=DateTime(9999), stop=DateTime(9999))
        end
    end
    @testset "CLay" begin
        mktempdir() do root
            src = joinpath(@__DIR__, "data", "caliop", "clay")
            cp.(joinpath.(src, mini_clay), joinpath.(root, mini_clay))
            sat = SatSet(root, type=:CLay)
            @test length(sat.granules) == 2
            @test sat.granules.lat isa Vector{Vector{Float32}}
            @test sat.granules.lon isa Vector{Vector{Float32}}
            @test sat.metadata.roots[0x0001] == realpath(root)
            @test sat.metadata.date.start == DateTime(2012, 2, 5, 5, 40, 3, 659)
            @test sat.metadata.date.stop == DateTime(2012, 2, 5, 7, 18, 50, 827)
            @test all(sat.metadata.granules.latmin .≈ [-69.07832, -81.821655])
            @test all(sat.metadata.granules.latmax .≈ [81.82099, 78.0213])
            @test all(sat.metadata.granules.elonmin .≈ [0.013624359, 53.707783])
            @test all(sat.metadata.granules.elonmax .≈ [77.15055, 179.803])
            @test all(sat.metadata.granules.wlonmin .≈ [-92.6476, -179.90944])
            @test all(sat.metadata.granules.wlonmax .≈ [-0.22770761, -92.70058])
            @test sat.metadata.type == :CLay
        end
        mktempdir() do root
            touch(joinpath(root, "CLay_01.h5"))
            sat = @test_logs(
                (:warn, "read error; data skipped"),
                (:warn, "no satellite data files successfully loaded"),
                SatSet(root, type=:CLay)
            )
            @test sat.metadata.date == (start=DateTime(9999), stop=DateTime(9999))
        end
    end
    @testset "constructors" begin
        sat = SatSet()
        @test sat.granules == StructArray{SatData{Float32}}(undef, 0)
        sat = SatData()
        @test sat.time == DateTime[]
        @test sat.lat == Float32[]
        @test sat.lon == Float32[]
        sd = SatData(joinpath(@__DIR__, "data", "caliop", "cpro", "CPro_01.h5"))
        st = SatTrack(joinpath(@__DIR__, "data", "caliop", "cpro", "CPro_01.h5"))
        st32 = SatTrack{Float32}(joinpath(@__DIR__, "data", "caliop", "cpro", "CPro_01.h5"))
        @test sd.time == st.time == st32.time
        @test sd.lat == st.lat == st32.lat
        @test sd.lon == st.lon == st32.lon
        sat = SecondaryMetadata()
        @test sat.type == :undef
        @test isempty(sat.roots)
        @test isempty(sat.loadtime.periods)
        @test isnothing(sat.attachments)
        @test sat.granules == DataFrame(file = String[], root = UInt16[], tstart = DateTime[], tstop = DateTime[],
            latmin = Float32[], latmax = Float32[], elonmin = Float32[], elonmax = Float32[],
            wlonmin = Float32[], wlonmax = Float32[])
        mktempdir() do root
            src = joinpath(@__DIR__, "data", "caliop", "cpro")
            cp.(joinpath.(src, mini_cpro), joinpath.(root, mini_cpro))
            sat = SatSet(root, type=:CPro)
            s64 = SatSet{Float64}(sat)
            sat64 = SatSet{Float64}(root, type=:CPro)
            sec64 = SecondarySet{Float64}(root, type=:CPro)
            sec32 = SecondarySet(root, type=:CPro)
            @test s64.granules.lat == sat64.granules.lat
            @test s64.granules.lon == sat64.granules.lon
            @test s64.granules.time == sat64.granules.time
            @test sec64.granules.lat == sat64.granules.lat
            @test sec64.granules.lon == sat64.granules.lon
            @test sec64.granules.time == sat64.granules.time
            @test sat64.granules.lat isa Vector{Vector{Float64}}
            @test sat64.granules.lon isa Vector{Vector{Float64}}
            @test s64.granules.lat isa Vector{Vector{Float64}}
            @test s64.granules.lon isa Vector{Vector{Float64}}
            @test sec32.granules.lat isa Vector{Vector{Float32}}
            @test sec32.granules.lon isa Vector{Vector{Float32}}
            @test sat64.metadata.date.start == DateTime(2012, 2, 5, 5, 40, 3, 659)
            @test sat64.metadata.date.stop == DateTime(2012, 2, 5, 7, 18, 50, 827)
            @test all(isapprox.(sat64.metadata.granules.latmin, [-69.07832, -81.821655]; atol=1e-5))
            @test all(isapprox.(sat64.metadata.granules.latmax, [81.82099, 78.0213]; atol=1e-5))
            @test all(isapprox.(sat64.metadata.granules.elonmin, [0.013624359, 53.707783]; atol=1e-5))
            @test all(isapprox.(sat64.metadata.granules.elonmax, [77.15055, 179.803]; atol=1e-5))
            @test all(isapprox.(sat64.metadata.granules.wlonmin, [-92.6476, -179.90944]; atol=1e-5))
            @test all(isapprox.(sat64.metadata.granules.wlonmax, [-0.22770761, -92.70058]; atol=1e-5))
            @test sat64.metadata.type == :CPro
        end
    end
end

@testset "Observations" begin
    timeindex = [2035:2065]
    @testset "CPro" begin
        lidarprofile = TrackMatcher.get_lidarheights((15_000, -Inf), Float32)
            cpro = CPro([joinpath(@__DIR__, "data", "caliop", "cpro", "CPro_23.h5")], timeindex, lidarprofile)
        cpro_empty = CPro([joinpath(@__DIR__, "data", "caliop", "cpro", "CPro_23.h5")],
            timeindex, lidarprofile, false)
        cpro64 = CPro{Float64}(cpro)
        @test cpro.time == utc
        @test cpro.lat == lat
        @test cpro.lon == lon
        @test cpro.atmos_state isa Vector{Vector{Enum{UInt16}}} && length(cpro.atmos_state) == 31
        @test cpro.EC532 isa Vector{Vector{<:Union{Missing,Float32}}} && length(cpro.EC532) == 31
        @test cpro.h_tropo isa Vector{Float32} && length(cpro.h_tropo) == 31
        @test cpro.temp isa Vector{<:Vector{<:Union{Missing,Float32}}} && length(cpro.temp) == 31
        @test cpro.pressure isa Vector{<:Vector{<:Union{Missing,Float32}}} && length(cpro.pressure) == 31
        @test cpro.rH isa Vector{<:Vector{<:Union{Missing,Float32}}} && length(cpro.rH) == 31
        @test cpro.IWC isa Vector{<:Vector{<:Union{Missing,Float32}}} && length(cpro.IWC) == 31
        @test cpro.deltap isa Vector{<:Vector{<:Union{Missing,Float32}}} && length(cpro.deltap) == 31
        @test cpro.CADscore isa Vector{<:Vector{<:Union{Missing,Int8}}} && length(cpro.CADscore) == 31
        @test cpro.night isa BitVector && length(cpro.night) == 31
        @test isempty(cpro_empty.time) && isempty(cpro_empty.lat) && isempty(cpro_empty.lon) && isempty(cpro_empty.atmos_state) &&
            isempty(cpro_empty.EC532) && isempty(cpro_empty.h_tropo) && isempty(cpro_empty.temp) &&
            isempty(cpro_empty.pressure) && isempty(cpro_empty.rH) && isempty(cpro_empty.IWC) &&
            isempty(cpro_empty.deltap) && isempty(cpro_empty.CADscore) && isempty(cpro_empty.night)
        @test cpro64.lat ≈ cpro.lat && cpro64.lat isa Vector{Float64}
        @test cpro64.lon ≈ cpro.lon && cpro64.lon isa Vector{Float64}
        @test cpro64.EC532 isa Vector{<:Vector{<:Union{Missing,Float64}}} &&
            nested_missing_aware_approx(cpro64.EC532, cpro.EC532)
        @test cpro64.h_tropo ≈ cpro.h_tropo && cpro64.h_tropo isa Vector{Float64}
        @test cpro64.temp isa Vector{<:Vector{<:Union{Missing,Float64}}} &&
            nested_missing_aware_approx(cpro64.temp, cpro.temp)
        @test cpro64.pressure isa Vector{<:Vector{<:Union{Missing,Float64}}} &&
            nested_missing_aware_approx(cpro64.pressure, cpro.pressure)
        @test cpro64.rH isa Vector{<:Vector{<:Union{Missing,Float64}}} &&
            nested_missing_aware_approx(cpro64.rH, cpro.rH)
        @test cpro64.IWC isa Vector{<:Vector{<:Union{Missing,Float64}}} &&
            nested_missing_aware_approx(cpro64.IWC, cpro.IWC)
        @test cpro64.deltap isa Vector{<:Vector{<:Union{Missing,Float64}}} &&
            nested_missing_aware_approx(cpro64.deltap, cpro.deltap)
        sat = @test_logs(
            (:error, "ReadError: profile observations could not be read from file, skipping"),
            CPro([joinpath(@__DIR__, "data", "caliop", "incorrect", "CPro.h5")], timeindex, lidarprofile)
        )
        @test isempty(sat.time) && isempty(sat.lat) && isempty(sat.lon) && isempty(sat.atmos_state) &&
            isempty(sat.EC532) && isempty(sat.h_tropo) && isempty(sat.temp) &&
            isempty(sat.pressure) && isempty(sat.rH) && isempty(sat.IWC) &&
            isempty(sat.deltap) && isempty(sat.CADscore) && isempty(sat.night)
    end
    @testset "CLay" begin
        clay = CLay([joinpath(@__DIR__, "data", "caliop", "clay", "CLay_23.h5")], timeindex, (15_000, -Inf))
        clay_empty = CLay([joinpath(@__DIR__, "data", "caliop", "clay", "CLay_23.h5")],
            timeindex, (15_000, -Inf), 5000, false)
        clay64 = CLay{Float64}(clay)
        @test clay.time == utc
        @test clay.lat == lat
        @test clay.lon == lon
        @test clay.layer_top isa Vector{Vector{Float32}} && length(clay.layer_top) == 31
        @test clay.layer_base isa Vector{Vector{Float32}} && length(clay.layer_base) == 31
        @test clay.atmos_state isa Vector{Vector{Enum{UInt16}}} && length(clay.atmos_state) == 31
        @test clay.OD isa Vector{Vector{Float32}} && length(clay.OD) == 31
        @test clay.IWP isa Vector{<:Vector{<:Union{Missing,Float32}}} && length(clay.IWP) == 31
        @test clay.Ttop isa Vector{Vector{Float32}} && length(clay.Ttop) == 31
        @test clay.h_tropo isa Vector{Float32} && length(clay.h_tropo) == 31
        @test clay.night isa BitVector && length(clay.night) == 31
        @test clay.averaging isa Vector{Int8} && length(clay.averaging) == 31
        @test isempty(clay_empty.time) && isempty(clay_empty.lat) && isempty(clay_empty.lon) && isempty(clay_empty.layer_top) &&
            isempty(clay_empty.layer_base) && isempty(clay_empty.atmos_state) && isempty(clay_empty.OD) &&
            isempty(clay_empty.IWP) && isempty(clay_empty.Ttop) && isempty(clay_empty.h_tropo) &&
            isempty(clay_empty.night) && isempty(clay_empty.averaging)
        @test clay64.lat ≈ clay.lat && clay64.lat isa Vector{Float64}
        @test clay64.lon ≈ clay.lon && clay64.lon isa Vector{Float64}
        @test clay64.layer_top ≈ clay.layer_top && clay64.layer_top isa Vector{Vector{Float64}}
        @test clay64.layer_base ≈ clay.layer_base && clay64.layer_base isa Vector{Vector{Float64}}
        @test clay64.OD ≈ clay.OD && clay64.OD isa Vector{Vector{Float64}}
        @test clay64.IWP isa Vector{<:Vector{<:Union{Missing,Float64}}} &&
            nested_missing_aware_approx(clay64.IWP, clay.IWP)
        @test clay64.Ttop ≈ clay.Ttop && clay64.Ttop isa Vector{Vector{Float64}}
        @test clay64.h_tropo ≈ clay.h_tropo && clay64.h_tropo isa Vector{Float64}
        sat = @test_logs(
            (:error, "ReadError: layer observations could not be read from file, skipping"),
            CLay([joinpath(@__DIR__, "data", "caliop", "incorrect", "CLay.h5")], timeindex, (15_000, -Inf))
        )
        @test isempty(sat.time) && isempty(sat.lat) && isempty(sat.lon) && isempty(sat.layer_top) &&
            isempty(sat.layer_base) && isempty(sat.atmos_state) && isempty(sat.OD) &&
            isempty(sat.IWP) && isempty(sat.Ttop) && isempty(sat.h_tropo) &&
            isempty(sat.night) && isempty(sat.averaging)
    end
end
