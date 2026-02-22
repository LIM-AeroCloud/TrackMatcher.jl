mini_cpro, mini_clay = ["CPro_01.h5", "CPro_02.h5"], ["CLay_01.h5", "CLay_02.h5"]
@testset "SatSet" begin
    @testset "CPro" begin
        mktempdir() do root
            src = joinpath("test", "data", "caliop", "cpro")
            cp.(joinpath.(src, mini_cpro), joinpath.(root, mini_cpro))
            sat =SatSet(root)
            @test length(sat.granules) == 2
            @test sat.granules.lat isa Vector{Vector{Float32}}
            @test sat.granules.lon isa Vector{Vector{Float32}}
            @test sat.metadata.roots[0x0001] == realpath(root)
            @test sat.metadata.date.start == DateTime(2012, 2, 5, 5, 40, 3, 659)
            @test sat.metadata.date.stop == DateTime(2012, 2, 5, 7, 18, 50, 827)
            @test all(sat.metadata.granules.latmin .≈ [Float32(-69.0783), Float32(-81.8217)])
            @test all(sat.metadata.granules.latmax .≈ [Float32(81.821), Float32(78.0213)])
            @test all(sat.metadata.granules.elonmin .≈ [Float32(0.0136244), Float32(53.7078)])
            @test all(sat.metadata.granules.elonmax .≈ [Float32(77.1506), Float32(179.803)])
            @test all(sat.metadata.granules.wlonmin .≈ [Float32(-92.6476), Float32(-179.909)])
            @test all(sat.metadata.granules.wlonmax .≈ [Float32(-0.227708), Float32(-92.7006)])
            @test sat.metadata.type == :CPro
        end
    end
    @testset "CLay" begin
        mktempdir() do root
            src = joinpath("test", "data", "caliop", "clay")
            cp.(joinpath.(src, mini_clay), joinpath.(root, mini_clay))
            sat =SatSet(root)
            @test length(sat.granules) == 2
            @test sat.granules.lat isa Vector{Vector{Float32}}
            @test sat.granules.lon isa Vector{Vector{Float32}}
            @test sat.metadata.roots[0x0001] == realpath(root)
            @test sat.metadata.date.start == DateTime(2012, 2, 5, 5, 40, 3, 659)
            @test sat.metadata.date.stop == DateTime(2012, 2, 5, 7, 18, 50, 827)
            @test all(sat.metadata.granules.latmin .≈ [Float32(-69.0783), Float32(-81.8217)])
            @test all(sat.metadata.granules.latmax .≈ [Float32(81.821), Float32(78.0213)])
            @test all(sat.metadata.granules.elonmin .≈ [Float32(0.0136244), Float32(53.7078)])
            @test all(sat.metadata.granules.elonmax .≈ [Float32(77.1506), Float32(179.803)])
            @test all(sat.metadata.granules.wlonmin .≈ [Float32(-92.6476), Float32(-179.909)])
            @test all(sat.metadata.granules.wlonmax .≈ [Float32(-0.227708), Float32(-92.7006)])
            @test sat.metadata.type == :CLay
        end
    end
end
