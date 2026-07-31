# ¡ Needs data from test_flightdata.jl, test_clouddata.jl, and test_satdata.jl to run

## Setup and helper functions
# Debug test data
# flight = FlightSet(volpe=joinpath(@__DIR__, "data", "volpe"))
# flight64 = FlightSet{Float64}(volpe=joinpath(@__DIR__, "data", "volpe"))
# flight_empty = FlightSet()
# cloud = CloudSet(joinpath(@__DIR__, "data", "cloud"))

# cpro_src = joinpath(@__DIR__, "data", "caliop", "CPro")
# clay_src = joinpath(@__DIR__, "data", "caliop", "CLay")
# cpro = SatSet(cpro_src, type=:CPro)
# clay = SatSet(clay_src, type=:CLay)

function xresults(intersections::XData{T}; obs::Vector{Bool}=[true, true, true], approx::Union{Bool,Int}=true) where T
    data = DataFrame(
        id=["V-2-1", "V-2-2"],
        lat=T[5.589219, 11.3293705],
        lon=T[14.462531, 12.077299],
        alt=T[11582.462, 11581.832],
        tdiff=Dates.CompoundPeriod[
            Dates.CompoundPeriod(Minute(-29), Second(-56)),
            Dates.CompoundPeriod(Minute(23), Second(34))
        ],
        tprim=[DateTime(2012, 2, 6, 0, 43, 13), DateTime(2012, 2, 6, 1, 27, 1)],
        tsec=[DateTime(2012, 2, 6, 0, 13, 17), DateTime(2012, 2, 6, 1, 50, 35)],
        atmos_state=[clear, ci]
    )
    accuracy = DataFrame(
        id=["V-2-1", "V-2-2"],
        intersection=T[1.0103808f6, 1.337841f6],
        primdist=T[89937.97f0, 39393.39f0],
        secdist=T[1559.8745, 3761.1736],
        primtime=[Dates.CompoundPeriod(Minute(-3), Second(-53)), Dates.CompoundPeriod(Second(21))],
        sectime=[Dates.CompoundPeriod(Millisecond(172)), Dates.CompoundPeriod(Millisecond(284))]
    )
    cpro = obs[2] ? CPro{T}([joinpath(@__DIR__, "data", "caliop", "CPro", "CPro_4.h5")],
        [2035:2065], TrackMatcher.get_lidarheights((15000, -Inf)), true) : CPro{T}()
    clay = obs[3] ? CLay{T}([joinpath(@__DIR__, "data", "caliop", "CLay", "CLay_4.h5")],
        [2035:2065]) : CLay{T}()
    cpro2 = obs[2] ? CPro{T}([joinpath(@__DIR__, "data", "caliop", "CPro", "CPro_6.h5")],
        [1906:1936], TrackMatcher.get_lidarheights((15000, -Inf)), true) : CPro{T}()
    clay2 = obs[3] ? CLay{T}([joinpath(@__DIR__, "data", "caliop", "CLay", "CLay_6.h5")],
        [1906:1936]) : CLay{T}()

    track = obs[1] ? FlightData{T}(
        (getproperty(flight.volpe[1], prop)[39:39] for prop in propertynames(flight.volpe[1])[1:(end-1)])...,
        flight.volpe.metadata[1]
    ) : FlightData{T}()
    track2 =  obs[1] ? FlightData{T}(
        (getproperty(flight.volpe[1], prop)[43:43] for prop in propertynames(flight.volpe[1])[1:(end-1)])...,
        flight.volpe.metadata[1]
    ) : FlightData{T}()
    observations = DataFrame(
        id = ["V-2-1", "V-2-2"],
        primary=[track, track2],
        CPro=[cpro, cpro2],
        CLay=[clay, clay2]
    )
    expected = XData{T}(data, observations, accuracy, intersections.metadata)
    results = if approx isa Int
        isapprox(intersections, expected, atol=10.0^-approx)
    elseif approx === true
        isapprox(intersections, expected)
    else
        intersections == expected
    end
    return results
end

## Test sets

@testset "intersections" begin
    # Run intercept finding routines
    xf_pro_lay = Intersection(flight, cpro, true)
    xf_lay = XData(flight, clay)
    # Test results
    @testset "data integrity" begin
        @test xresults(xf_pro_lay, approx=false)
        @test xresults(xf_lay, obs=[true, false, true], approx=false)
    end
end

#=
T = Float32



c = xf.observations.CPro[1]

@testset "intercept finding" begin
    @testset "VOLPE" begin
        xfcpro = XData(flight, cpro)
        xccpro = Intersection(cloud, cpro)
        xf_promoted = XData(flight64, cpro)
        xf64 = Intersection{Float64}(flight, cpro)
        xfclay = Intersection(flight, clay)
        xcclay = XData(cloud, clay)
        xc64 = XData{Float64}(cloud, clay, true)
        xfcproclay = Intersection(flight, cpro, true)
        xcclaycpro = XData(cloud, clay, true)
        xfempty = XData(flight_empty, cpro, true)
        @testset "data integrity" begin
            @test length(flight.volpe) == 2
            @test isempty(flight.flightaware)
            @test isempty(flight.webdata)
            @test minimum([flight.volpe.alt...;]) ≥ 5000
            @test all(length.(getproperty.(Ref(flight.volpe[1]), propertynames(flight.volpe))[1:(end-1)]) .== 77)
            @test flight.volpe.time isa Vector{<:Vector{DateTime}}
            @test flight.volpe.lat isa Vector{<:Vector{Float32}}
            @test flight.volpe.lon isa Vector{<:Vector{Float32}}
            @test flight.volpe.alt isa Vector{<:Vector{<:Union{Missing,Float32}}}
            @test flight.volpe.metadata[1].id == 2 # ℹ flight 1 is filtered out because of altmin=5000
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
            primary = PrimarySet(volpe=joinpath(@__DIR__, "data", "volpe", "hit"))
            primary64 = PrimarySet{Float64}(volpe=joinpath(@__DIR__, "data", "volpe", "hit"))
            flightprimary = PrimarySet{Float64}(flight)
            alt_atol = 2e-3
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
                (:error, r"Unable to parse date. Data skipped."),
                (:error, r"Unknown file format"),
                (:warn, r"Unknown time zone."),
                (:warn, r"Unknown time zone format"),
                (:error, r"Invalid file name format. Data skipped."),
                (:error, r"Error reading file. Try to specify column delimiter. Data skipped."),
                (:info, r"FlightSet loaded"),
                FlightSet(webdata=joinpath(@__DIR__, "data", "webdata"))
            )
            webok = @test_logs(
                (:info, r"FlightSet loaded"),
                FlightSet(webdata=joinpath(@__DIR__, "data", "webdata", "ok"), delim=('\t'))
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
            @test web.webdata.metadata[1].date.start ==
                  DateTime(ZonedDateTime(2016, 09, 24, 2, 28, 41, localzone()), UTC) &&
                  web.webdata.metadata[1].date.stop ==
                  DateTime(ZonedDateTime(2016, 09, 24, 7, 27, 07, localzone()), UTC)
        end
    end
    flight_all = FlightSet(
        volpe=joinpath(@__DIR__, "data", "volpe", "hit"),
        flightaware=joinpath(@__DIR__, "data", "archive", "new"),
        webdata=joinpath(@__DIR__, "data", "webdata", "ok"), delim=('\t')
    )
    @test length(flight_all.volpe) == 2 && length(flight_all.flightaware) == 3 && length(flight_all.webdata) == 3
end
=#
