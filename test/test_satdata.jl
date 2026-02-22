mini_cpro, mini_clay = ["CPro_01.h5", "CPro_02.h5"], ["CLay_01.h5", "CLay_02.h5"]
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
