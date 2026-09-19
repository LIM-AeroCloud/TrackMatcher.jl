# ¡ Needs general test data from init.jl
## Setup helper functions

function xdata_matches(
    intersections::XData{T},
    expected::XData{T},
    obs::Vector{Bool};
    data_atol::Real=1e-3,
    accuracy_atol::Real=2e-1
) where T
    intersections.data.id == expected.data.id || return false
    isapprox(intersections.data.lat, expected.data.lat; atol=data_atol) || return false
    isapprox(intersections.data.lon, expected.data.lon; atol=data_atol) || return false
    isapprox(intersections.data.alt, expected.data.alt; atol=data_atol) || return false
    intersections.data.tdiff == expected.data.tdiff || return false
    intersections.data.tprim == expected.data.tprim || return false
    intersections.data.tsec == expected.data.tsec || return false
    intersections.data.atmos_state == expected.data.atmos_state || return false

    intersections.accuracy.id == expected.accuracy.id || return false
    isapprox(intersections.accuracy.intersection, expected.accuracy.intersection;
        atol=accuracy_atol) || return false
    isapprox(intersections.accuracy.primdist, expected.accuracy.primdist;
        atol=accuracy_atol) || return false
    isapprox(intersections.accuracy.secdist, expected.accuracy.secdist;
        atol=accuracy_atol) || return false
    intersections.accuracy.primtime == expected.accuracy.primtime || return false
    intersections.accuracy.sectime == expected.accuracy.sectime || return false

    intersections.observations.id == expected.observations.id || return false
    names(intersections.observations) == names(expected.observations) || return false
    obs[1] || all(isempty, intersections.observations.primary) || return false
    obs[2] || all(isempty, intersections.observations.CPro) || return false
    obs[3] || all(isempty, intersections.observations.CLay) || return false
    return true
end


function xresults(
    intersections::XData{T},
    primary::PrimarySet{T},
    obsdata::Tuple{CPro,CLay}=(cpro, clay);
    obs::Vector{Bool}=[true, true, true],
    approx::Union{Bool,Int}=true
) where T
    data = primary isa FlightSet ? DataFrame(
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
    ) : DataFrame(
        id=["C-1-1"],
        lat=T[5.24175],
        lon=T[23.464722],
        alt=T[NaN],
        tdiff=Dates.CompoundPeriod[Dates.CompoundPeriod(Minute(-25), Second(-6))],
        tprim=[DateTime(2012, 2, 6, 0, 38, 29)],
        tsec=[DateTime(2012, 2, 6, 0, 13, 23)],
        atmos_state=[clear]
    )
    accuracy = primary isa FlightSet ? DataFrame(
        id=["V-2-1", "V-2-2"],
        intersection=T[1.0103808f6, 1.337841f6],
        primdist=T[89937.97f0, 39393.39f0],
        secdist=T[1559.8745, 3761.1736],
        primtime=[Dates.CompoundPeriod(Minute(-3), Second(-53)), Dates.CompoundPeriod(Second(21))],
        sectime=[Dates.CompoundPeriod(Millisecond(172)), Dates.CompoundPeriod(Millisecond(284))]
    ) : DataFrame(
        id=["C-1-1"],
        intersection=T[34.18339],
        primdist=T[NaN],
        secdist=T[NaN],
        primtime=[Dates.CompoundPeriod()],
        sectime=[Dates.CompoundPeriod()]
    )
    observations = primary isa FlightSet ? DataFrame(
        id=["V-2-1", "V-2-2"],
        # Non-empty dummy data for observations
        primary=flight.volpe,
        CPro=[obsdata[1], obsdata[1]],
        Clay=[obsdata[2], obsdata[2]]
    ) : DataFrame(
        id=["C-1-1"],
        primary=flight[1:1],
        CPro=[obsdata[1]],
        Clay=[obsdata[2]]
    )
    expected = XData{T}(data, observations, accuracy, intersections.metadata)
    results = if approx isa Int
        xdata_matches(intersections, expected, obs;
            data_atol=10.0^-approx, accuracy_atol=10.0^-approx)
    elseif approx === true
        xdata_matches(intersections, expected, obs;
            data_atol=1e-3, accuracy_atol=1e-1)
    else
        xdata_matches(intersections, expected, obs;
            data_atol=1e-3, accuracy_atol=1e-1)
    end
    return results
end

## Test sets

@testset "intersections" begin
    @testset "longitude-axis intersection coordinates" begin
        primary = (track = x -> 2x, min = 0.0, max = 2.0)
        secondary = (track = x -> 3 .- x, min = 0.0, max = 2.0)

        primary_coords, secondary_coords = TrackMatcher.findXcoords(
            primary, secondary, 0.1, true, Float64)

        @test length(primary_coords) == 1
        @test length(secondary_coords) == 1
        @test primary_coords[1][1] ≈ 2.0
        @test primary_coords[1][2] ≈ 1.0
        @test secondary_coords[1][1] ≈ 2.0
        @test secondary_coords[1][2] ≈ 1.0
    end

    # Run intercept finding routines
    xf_cpro = Intersection(flight, sat_cpro) # ℹ reused in constructor testset
    xf_clay = XData(flight, sat_clay, true)
    xf64 = Intersection{Float64}(xf_cpro)
    xf64_promoted = Intersection(flight64, sat_cpro)
    xf64_forced = XData{Float64}(flight, sat_cpro)
    xc_pro = XData(cloud, sat_cpro)
    xc_lay = XData(cloud, sat_clay)
    xf_empty = XData(flight_empty, sat_cpro)
    # Test results
    @testset "data integrity" begin
        @test xresults(xf_cpro, flight, obs=[true, true, false], approx=false)
        @test xresults(xf_clay, flight, approx=false)
        @test xf64 isa XData{Float64}
        @test xf64 ≈ xf_cpro
        @test xf64_promoted isa XData{Float64}
        @test xdata_matches(xf64_promoted, xf64, [true, true, false])
        @test xf64_forced isa XData{Float64}
        @test xdata_matches(xf64_forced, xf64, [true, true, false])
        @test isempty(xf_empty)
        @test xc_lay isa XData{Float32}
        @test xc_pro isa XData{Float32}
    end
    @testset "promotion constructors preserve original data" begin
        original_data = deepcopy(xf_cpro.data)
        original_observations = deepcopy(xf_cpro.observations)
        original_accuracy = deepcopy(xf_cpro.accuracy)
        original_metadata = deepcopy(xf_cpro.metadata)

        promoted = XData{Float64}(xf_cpro)

        @test xf_cpro.data == original_data
        @test xf_cpro.observations == original_observations
        @test xf_cpro.accuracy == original_accuracy
        @test xf_cpro.metadata == original_metadata
        @test promoted isa XData{Float64}

        track = flight.volpe[1]
        promoted_track = FlightData{Float64}(track)
        @test promoted_track.time !== track.time
        @test promoted_track.lat !== track.lat
        @test promoted_track.metadata !== track.metadata

        sat = sat_cpro.granules[1]
        promoted_sat = SatData{Float64}(sat)
        @test promoted_sat.time !== sat.time
        @test promoted_sat.lat !== sat.lat

        cpro_obs = xf_cpro.observations.CPro[1]
        promoted_cpro = CPro{Float64}(cpro_obs)
        @test promoted_cpro.time !== cpro_obs.time
        @test promoted_cpro.lat !== cpro_obs.lat

        clay_obs = CLay{Float32}()
        promoted_clay = CLay{Float64}(clay_obs)
        @test promoted_clay.time !== clay_obs.time
        @test promoted_clay.lat !== clay_obs.lat
    end
    @testset "constructors" begin
        # Instantiate with convenience constructors
        mdata = MeasuredData([
            "volpe" => joinpath(@__DIR__, "data", "volpe"),
            "cloudtracks" => joinpath(@__DIR__, "data", "cloud"),
            "sat" => cpro_src
        ])
        mset = MeasuredSet([
            "volpe" => joinpath(@__DIR__, "data", "volpe"),
            "cloudtracks" => joinpath(@__DIR__, "data", "cloud"),
            "sat" => cpro_src
        ])
        mdata64 = MeasuredData{Float64}(mdata)
        data = Data([
            "volpe" => joinpath(@__DIR__, "data", "volpe"),
            "cloudtracks" => joinpath(@__DIR__, "data", "cloud"),
            "sat" => cpro_src
        ])
        dset = DataSet([
            "volpe" => joinpath(@__DIR__, "data", "volpe"),
            "cloudtracks" => joinpath(@__DIR__, "data", "cloud"),
            "sat" => cpro_src
        ])
        data64 = Data{Float64}(data)
        # Test convenience constructors
        @test mdata == mset
        @test mdata.flight.volpe isa StructArray{FlightData{Float32}} && length(mdata.flight.volpe) == 2
        @test mdata.flight.flightaware isa StructArray{FlightData{Float32}} && isempty(mdata.flight.flightaware)
        @test mdata.flight.webdata isa StructArray{FlightData{Float32}} && isempty(mdata.flight.webdata)
        @test mdata.cloud.tracks isa StructArray{CloudData{Float32}} && length(mdata.cloud.tracks) == 3
        @test mdata.sat.granules isa StructArray{SatData{Float32}} && length(mdata.sat.granules) == 7

        @test mdata64 isa MeasuredData{Float64}
        @test mdata64.flight.volpe isa StructArray{FlightData{Float64}} && length(mdata64.flight.volpe) == 2
        @test mdata64.flight.flightaware isa StructArray{FlightData{Float64}} && isempty(mdata64.flight.flightaware)
        @test mdata64.flight.webdata isa StructArray{FlightData{Float64}} && isempty(mdata64.flight.webdata)
        @test mdata64.cloud.tracks isa StructArray{CloudData{Float64}} && length(mdata64.cloud.tracks) == 3
        @test mdata64.sat.granules isa StructArray{SatData{Float64}} && length(mdata64.sat.granules) == 7

        @test data == dset
        @test data.trackdata.flight.volpe isa StructArray{FlightData{Float32}} &&
            length(data.trackdata.flight.volpe) == 2
        @test data.trackdata.flight.flightaware isa StructArray{FlightData{Float32}} &&
            isempty(data.trackdata.flight.flightaware)
        @test data.trackdata.flight.webdata isa StructArray{FlightData{Float32}} &&
            isempty(data.trackdata.flight.webdata)
        @test data.trackdata.cloud.tracks isa StructArray{CloudData{Float32}} &&
            length(data.trackdata.cloud.tracks) == 3
        @test data.trackdata.sat.granules isa StructArray{SatData{Float32}} &&
            length(data.trackdata.sat.granules) == 7
        @test data.intersection.flight isa XData{Float32} && size(data.intersection.flight.data) == (2, 8) &&
            size(data.intersection.flight.observations) == (2, 4) && size(data.intersection.flight.accuracy) == (2, 6)
        @test data.intersection.cloud isa XData{Float32} && size(data.intersection.cloud.data) == (1, 8) &&
            size(data.intersection.cloud.observations) == (1, 4) && size(data.intersection.cloud.accuracy) == (1, 6)

        @test data64 isa Data{Float64}
        @test data64.trackdata.flight.volpe isa StructArray{FlightData{Float64}} &&
            length(data64.trackdata.flight.volpe) == 2
        @test data64.trackdata.flight.flightaware isa StructArray{FlightData{Float64}} &&
            isempty(data64.trackdata.flight.flightaware)
        @test data64.trackdata.flight.webdata isa StructArray{FlightData{Float64}} &&
            isempty(data64.trackdata.flight.webdata)
        @test data64.trackdata.cloud.tracks isa StructArray{CloudData{Float64}} &&
            length(data64.trackdata.cloud.tracks) == 3
        @test data64.trackdata.sat.granules isa StructArray{SatData{Float64}} &&
            length(data64.trackdata.sat.granules) == 7
        @test data64.intersection.flight isa XData{Float64} && data64.intersection.cloud isa XData{Float64}

        @test XMetadata(getfield.(Ref(xf_cpro.metadata), fieldnames(XMetadata))...) isa XMetadata{Float32}
    end
    @testset "exception handling" begin
        # Force an exception in the interpolation path
        TrackMatcher.interpolate_trackdata(::FlightTrack) =
            throw(ErrorException("forced failure"))
        @test_logs min_level = Debug match_mode = :all (
            :debug, "primary track ID: 2"
        ) (
            :warn, r"Track data and/or time could not be interpolated"
        ) (
            :info, r"Intersection data \(0 matches\) loaded"
        ) begin
            x = XData(flight, sat_cpro, true)
            @test x isa XData
            @test isempty(x)
        end
        overlap, range = @test_logs (
            :warn, r"no sufficient satellite data"
        ) TrackMatcher.findoverlap(flight.volpe[1], SatSet(), 30, 0.1)
        @test overlap isa Vector{DataFrame} && isempty(overlap)
        @test range ==0:-1
    end
end
