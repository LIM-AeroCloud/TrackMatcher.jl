@testset "flight data" begin
    @testset "VOLPE" begin
        flight = FlightSet(volpe=joinpath(@__DIR__, "data", "volpe"))
        flight64 = FlightSet{Float64}(volpe=joinpath(@__DIR__, "data", "volpe"))
        flight16 = FlightSet{Float16}(flight64)
        t0 = now()
        flight_empty = FlightSet()
        flight_nothing = FlightSet(volpe=joinpath(@__DIR__, "data", "caliop", "clay"))
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
            @test t0 ≤ flight_empty.metadata.date.start == flight_empty.metadata.date.stop ≤
                DateTime(flight_empty.metadata.created)
            @test isempty(flight_nothing.volpe)
            @test t0 ≤ flight_nothing.metadata.date.start == flight_nothing.metadata.date.stop ≤
                DateTime(flight_nothing.metadata.created)
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
            # Define datasets
            primary = PrimarySet(volpe=joinpath(@__DIR__, "data", "volpe"))
            primary64 = PrimarySet{Float64}(volpe=joinpath(@__DIR__, "data", "volpe"))
            flightprimary = PrimarySet{Float64}(flight)
            alt_atol = 2e-3
            approx_vec(v1, v2; atol=1e-5) = all(((ismissing(x) || ismissing(y)) ?
                x === y : isapprox(x, y; atol=atol, rtol=0)) for (x, y) in zip(v1, v2))
            fields, fields64 = [], []
            for field in fieldnames(FlightData)
                push!(fields, getproperty(flight.volpe[1], field))
                push!(fields64, getproperty(flight64.volpe[1], field))
            end
            track = FlightTrack(fields...)
            track64 = FlightTrack{Float64}(fields64...)
            track_converted = FlightTrack{Float64}(track)
            track_empty = FlightTrack()
            flighttrack = FlightData(fields...)
            primmeta = PrimaryMetadata()
            flightmeta = FlightMetadata()

            # Test constructors and datasets
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
            @test track_empty isa FlightData{Float32}
            @test isempty(track_empty.time) && isempty(track_empty.lat) &&
                isempty(track_empty.lon) && isempty(track_empty.alt)
            @test flighttrack isa FlightData{Float32}
            @test flighttrack.lat == flight.volpe.lat[1]
            @test flighttrack.lon == flight.volpe.lon[1]
            @test flighttrack.alt == flight.volpe.alt[1]
            @test flighttrack.lat isa Vector{Float32} && flighttrack.lon isa Vector{Float32} &&
                flighttrack.alt isa Vector{<:Union{Missing,Float32}}
            @test primmeta isa PrimaryMetadata{Float32}
            @test flightmeta isa FlightMetadata{Float32}
            @test isempty(flightmeta.id)
        end
    end
    @testset "FlightAware" begin
        flight_old = FlightSet(flightaware=joinpath(@__DIR__, "data", "archive", "old"))
        flight_new = FlightSet(flightaware=joinpath(@__DIR__, "data", "archive", "new"))
        @test length(flight_old.flightaware) == 9
        @test isempty(flight_old.volpe)
        @test isempty(flight_old.webdata)
        @test minimum(skipmissing([flight_old.flightaware.alt...;])) ≥ 5000
        @test length(flight_new.flightaware) == 3
        @test isempty(flight_new.volpe)
        @test isempty(flight_new.webdata)
        @test minimum(skipmissing([flight_new.flightaware.alt...;])) ≥ 5000
        @test flight_old.flightaware.lat isa Vector{<:Vector{Float32}} &&
            flight_old.flightaware.lon isa Vector{<:Vector{Float32}} &&
            flight_old.flightaware.alt isa Vector{<:Vector{<:Union{Missing,Float32}}}
        @test flight_new.flightaware.lat isa Vector{<:Vector{Float32}} &&
            flight_new.flightaware.lon isa Vector{<:Vector{Float32}} &&
            flight_new.flightaware.alt isa Vector{<:Vector{<:Union{Missing,Float32}}}
    end
    @testset "web data" begin
        # Load web data
        web, webok = nothing, nothing
        @testset "error handling and logs" begin
            # Test that loading files with incorrect names produces warnings and skips data
            web = @test_logs(
                (:warn, r"Unknown time zone."),
                (:warn, r"Unknown time zone format"),
                (:error, r"Invalid file name format. Data skipped."),
                (:error, r"Error reading file. Try to specify column delimiter. Data skipped."),
                (:info, r"FlightSet loaded"),
                FlightSet(webdata=joinpath(@__DIR__, "data", "webdata"))
            )
            webok = @test_logs(
                (:info, r"FlightSet loaded"),
                FlightSet(webdata=joinpath(@__DIR__, "data", "webdata", "ok"), delim='\t')
            )
        end
        # Test data integrity and metadata
        @testset "data integrity and metadata" begin
            @test length(web.webdata) == 4
            @test isempty(web.volpe) && isempty(web.flightaware)
            @test length(webok.webdata) == 3
            @test minimum(skipmissing([webok.webdata.alt...;])) ≥ 5000
            @test webok.webdata.lat isa Vector{<:Vector{Float32}} &&
                webok.webdata.lon isa Vector{<:Vector{Float32}} &&
                webok.webdata.alt isa Vector{<:Vector{<:Union{Missing,Float32}}}
            @test webok.webdata.lat isa Vector{<:Vector{Float32}} &&
                webok.webdata.lon isa Vector{<:Vector{Float32}} &&
                webok.webdata.alt isa Vector{<:Vector{<:Union{Missing,Float32}}}
            @test web.webdata.metadata[1].flight_num == "ABC123"
            @test web.webdata.metadata[1].route.orig == "ABCD" && web.webdata.metadata[1].route.dest == "EFGH"

            @test webok.webdata.metadata[1].date.start == DateTime(2010, 6, 6, 7, 40, 40) &&
                webok.webdata.metadata[1].date.stop == DateTime(2010, 6, 6, 12, 36, 29)
            @test web.webdata.metadata[1].date.start == Dates.DateTime(2016, 09, 24, 0, 28, 41) &&
                web.webdata.metadata[1].date.stop == Dates.DateTime(2016, 09, 24, 5, 27, 07)
        end
    end
    flight = FlightSet(
        volpe=joinpath(@__DIR__, "data", "volpe"),
        flightaware=joinpath(@__DIR__, "data", "archive", "new"),
        webdata=joinpath(@__DIR__, "data", "webdata", "ok"), delim='\t'
    )
    @test length(flight.volpe) == 13 && length(flight.flightaware) == 3 && length(flight.webdata) == 3
end
