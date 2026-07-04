@testset "flight data" begin
    @testset "VOLPE" begin
        flight = FlightSet(volpe=joinpath(@__DIR__, "..", "test", "data", "volpe"))
        flight64 = FlightSet{Float64}(volpe=joinpath(@__DIR__, "..", "test", "data", "volpe"))
        flight16 = FlightSet{Float16}(flight64)
        flight_empty = FlightSet()
        @testset "data integrity" begin
            @test length(flight.volpe) == 13
            @test isempty(flight.flightaware)
            @test isempty(flight.webdata)
            @test minimum([flight.volpe.alt...;]) ≥ 5000
            @test all(length.(getproperty.(Ref(flight.volpe[2]), propertynames(flight.volpe))[1:end-1]) .== 77)
            @test flight.volpe.time isa Vector{<:Vector{DateTime}}
            @test flight.volpe.lat isa Vector{<:Vector{Float32}}
            @test flight.volpe.lon isa Vector{<:Vector{Float32}}
            @test flight.volpe.alt isa Vector{<:Vector{<:Union{Missing,Float32}}}
            @test isempty(flight_empty.volpe) && isempty(flight_empty.flightaware) &&
                isempty(flight_empty.webdata)
        end
        @testset "data precision" begin
            @test flight64.volpe.time isa Vector{<:Vector{DateTime}}
            @test flight64.volpe.lat isa Vector{<:Vector{Float64}}
            @test flight16.volpe.lat isa Vector{<:Vector{Float16}}
            @test flight64.volpe.lon isa Vector{<:Vector{Float64}}
            @test flight16.volpe.lon isa Vector{<:Vector{Float16}}
            @test flight64.volpe.alt isa Vector{<:Vector{<:Union{Missing,Float64}}}
            @test flight16.volpe.alt isa Vector{<:Vector{<:Union{Missing,Float16}}}
            @test all(flight.volpe.lat[1] .≈ flight64.volpe.lat[1])
            @test all(flight.volpe.lon[1] .≈ flight64.volpe.lon[1])
            @test all(flight.volpe.alt[1] .≈ flight64.volpe.alt[1])
        end
        @testset "constructors" begin
            primary = PrimarySet(volpe=joinpath(@__DIR__, "..", "test", "data", "volpe"))
            primary64 = PrimarySet{Float64}(volpe=joinpath(@__DIR__, "..", "test", "data", "volpe"))
            flightprimary = PrimarySet{Float64}(flight)
            alt_atol = 2e-3
            approx_vec(v1, v2; atol=1e-5) = all(((ismissing(x) || ismissing(y)) ?
                x === y : isapprox(x, y; atol=atol, rtol=0)) for (x, y) in zip(v1, v2))
            fields = []
            for field in fieldnames(FlightData)
                push!(fields, getproperty(flight.volpe[1], field))
            end
            track = FlightTrack(fields...)
            track64 = FlightTrack{Float64}(track)
            track_converted = FlightTrack{Float64}(track)
            @test primary isa FlightSet{Float32}
            @test primary64 isa FlightSet{Float64}
            @test flightprimary isa FlightSet{Float64}
            @test primary.volpe.time == flight.volpe.time
            @test primary.volpe.lat == flight.volpe.lat
            @test primary.volpe.lon == flight.volpe.lon
            @test primary.volpe.alt == flight.volpe.alt
            @test primary64.volpe.time == flight64.volpe.time
            @test primary64.volpe.lat == flight64.volpe.lat
            @test primary64.volpe.lon == flight64.volpe.lon
            @test primary64.volpe.alt == flight64.volpe.alt
            @test flightprimary.volpe.time == flight64.volpe.time
            @test all(approx_vec(a, b) for (a, b) in zip(flightprimary.volpe.lat, flight64.volpe.lat))
            @test all(approx_vec(a, b) for (a, b) in zip(flightprimary.volpe.lon, flight64.volpe.lon))
            @test all(approx_vec(a, b; atol=alt_atol) for (a, b) in zip(flightprimary.volpe.alt, flight64.volpe.alt))
            @test track isa FlightData{Float32}
            @test track64 isa FlightData{Float64}
            @test track_converted isa FlightData{Float64}
            @test track.time == flight.volpe.time[1]
            @test track.lat == flight.volpe.lat[1]
            @test track.lon == flight.volpe.lon[1]
            @test track.alt == flight.volpe.alt[1]
            @test track64.time == flight64.volpe.time[1]
            @test approx_vec(track64.lat, flight64.volpe.lat[1])
            @test approx_vec(track64.lon, flight64.volpe.lon[1])
            @test approx_vec(track64.alt, flight64.volpe.alt[1]; atol=alt_atol)
            @test track_converted.time == flight64.volpe.time[1]
            @test approx_vec(track_converted.lat, flight64.volpe.lat[1])
            @test approx_vec(track_converted.lon, flight64.volpe.lon[1])
            @test approx_vec(track_converted.alt, flight64.volpe.alt[1]; atol=alt_atol)
        end
    end
end
