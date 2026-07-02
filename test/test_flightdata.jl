@testset "flight data" begin
    flight = FlightSet(volpe=joinpath(@__DIR__, "..", "test", "data", "volpe"))
    @test length(flight.volpe) == 13
    @test isempty(flight.flightaware)
    @test isempty(flight.webdata)
    @test all([flight.volpe.alt...;] .≥ 5000)
    @test all(length.(getproperty.(Ref(flight.volpe[2]), propertynames(flight.volpe))[1:end-1]) .== 77)
    @test flight.volpe.time isa Vector{<:Vector{DateTime}}
    @test flight.volpe.lat isa Vector{<:Vector{Float32}}
    @test flight.volpe.lon isa Vector{<:Vector{Float32}}
    @test flight.volpe.alt isa Vector{<:Vector{<:Union{Missing,Float32}}}
end
