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
    atmos_state = intersections.metadata.sattype == :CLay ?
        [invalid, ci, invalid, invalid, invalid] : [clear, ci, clear, clear, clear]
    data = primary isa FlightSet ? DataFrame(
        id=["V-42-1", "V-42-2", "V-42-3", "V-42-4", "V-43-1"],
        lat=T[20.850319, 29.305733, 40.152092, 49.438873, 49.438488],
        lon=T[26.92516, 29.0, 32.091255, 35.443886, 35.44373],
        alt=T[11004.194, 11004.804, 11507.114, 11502.542, 11503.152],
        tdiff=Dates.CompoundPeriod[
            Dates.CompoundPeriod(Minute(2), Second(14)),
            Dates.CompoundPeriod(Minute(1), Second(36)),
            Dates.CompoundPeriod(Second(45)),
            Dates.CompoundPeriod(),
            Dates.CompoundPeriod()
        ],
        tprim=[DateTime(2012, 2, 6, 0, 6, 50), DateTime(2012, 2, 6, 0, 5, 8),
            DateTime(2012, 2, 6, 0, 2, 58), DateTime(2012, 2, 6, 0, 1, 7),
            DateTime(2012, 2, 6, 0, 1, 7)],
        tsec=[DateTime(2012, 2, 6, 0, 9, 4), DateTime(2012, 2, 6, 0, 6, 44),
            DateTime(2012, 2, 6, 0, 3, 43), DateTime(2012, 2, 6, 0, 1, 7),
            DateTime(2012, 2, 6, 0, 1, 7)],
        atmos_state=atmos_state
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
        id=["V-42-1", "V-42-2", "V-42-3", "V-42-4", "V-43-1"],
        intersection=T[0.21223556, 0.0, 0.0, 0.4238313, 0.3452892],
        primdist=T[96857.97, 67541.35, 17558.139, 51154.125, 0.0],
        secdist=T[2590.4082, 872.12195, 1785.7188, 479.24576, 440.39847],
        primtime=[Dates.CompoundPeriod(Second(-10)), Dates.CompoundPeriod(Second(8)),
            Dates.CompoundPeriod(Second(-2)), Dates.CompoundPeriod(Second(7)),
            Dates.CompoundPeriod()],
        sectime=[Dates.CompoundPeriod(Millisecond(125)), Dates.CompoundPeriod(Millisecond(-7)),
            Dates.CompoundPeriod(Millisecond(-220)), Dates.CompoundPeriod(Millisecond(16)),
            Dates.CompoundPeriod(Millisecond(16))]
    ) : DataFrame(
        id=["C-1-1"],
        intersection=T[34.18339],
        primdist=T[NaN],
        secdist=T[NaN],
        primtime=[Dates.CompoundPeriod()],
        sectime=[Dates.CompoundPeriod()]
    )
    observations = primary isa FlightSet ? DataFrame(
        id=["V-42-1", "V-42-2", "V-42-3", "V-42-4", "V-43-1"],
        # Non-empty dummy data for observations
        primary=[flight.volpe[1], flight.volpe[1], flight.volpe[1], flight.volpe[1], flight.volpe[2]],
        CPro=[obsdata[1] for _ in 1:5],
        Clay=[obsdata[2] for _ in 1:5]
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
    @testset "PCHIP interpolation" begin
        # Interpolation of scalar points
        coorddist, coorddist_derivative = TrackMatcher.pchip_difference_callbacks(
            [0.0, 1.0, 2.0], [0.0, 1.0, 4.0])
        x = 0.5

        @test coorddist(x) ∈ coorddist(x .. x)
        @test coorddist_derivative(x) ∈ coorddist_derivative(x .. x)

        # Interpolation of interval enclosure
        primary = (track=x -> zero.(x), min=0.0, max=2.0)
        secondary = (track=x -> 1 .- 2 .* (x .- 1) .^ 2, min=0.0, max=2.0)

        primary_coords, secondary_coords = TrackMatcher.findXcoords(
            primary, secondary, 1.0, true, Float64)

        @test length(primary_coords) == 2
        @test length(secondary_coords) == 2
        @test all(0.0 .< getindex.(primary_coords, 2) .< 2.0)
    end
    @testset "data integrity" begin
        overlap, isat = TrackMatcher.findoverlap(flight.volpe[2], sat_cpro, 30)
        lon_tracks = TrackMatcher.interpolate_trackdata(flight.volpe[2])
        lon_sat_tracks = TrackMatcher.interpolate_satdata(overlap, isat, true)
        lon_primary, lon_secondary = TrackMatcher.findXcoords(
            lon_tracks[1], lon_sat_tracks[1], 0.01, true, Float32)
        @test !isempty(lon_primary)
        @test !isempty(lon_secondary)
        @test xresults(xf_cpro, flight, obs=[true, true, false], approx=false)
        @test xresults(xf_clay, flight, approx=false)
        @test all(id -> !startswith(id, "V-44"), xf_cpro.data.id)
        @test xf64 isa XData{Float64}
        @test xf64 ≈ xf_cpro
        @test xf64_promoted isa XData{Float64}
        @test isapprox(xf64_promoted.data.lat, xf64.data.lat; atol=1e-3)
        @test isapprox(xf64_promoted.data.lon, xf64.data.lon; atol=1e-3)
        @test xf64_promoted.data.tprim == xf64.data.tprim
        @test xf64_promoted.data.tsec == xf64.data.tsec
        @test xf64_promoted.data.atmos_state == xf64.data.atmos_state
        @test xf64_forced isa XData{Float64}
        @test xf64_forced.data.id == xf64.data.id
        @test isapprox(xf64_forced.data.lat, xf64.data.lat; atol=1e-3)
        @test isapprox(xf64_forced.data.lon, xf64.data.lon; atol=1e-3)
        @test xf64_forced.data.tprim == xf64.data.tprim
        @test xf64_forced.data.tsec == xf64.data.tsec
        @test xf64_forced.data.atmos_state == xf64.data.atmos_state
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
        @test mdata.flight.volpe isa StructArray{FlightData{Float32}} && length(mdata.flight.volpe) == 3
        @test mdata.flight.flightaware isa StructArray{FlightData{Float32}} && isempty(mdata.flight.flightaware)
        @test mdata.flight.webdata isa StructArray{FlightData{Float32}} && isempty(mdata.flight.webdata)
        @test mdata.cloud.tracks isa StructArray{CloudData{Float32}} && length(mdata.cloud.tracks) == 3
        @test mdata.sat.granules isa StructArray{SatData{Float32}} && length(mdata.sat.granules) == 7

        @test mdata64 isa MeasuredData{Float64}
        @test mdata64.flight.volpe isa StructArray{FlightData{Float64}} && length(mdata64.flight.volpe) == 3
        @test mdata64.flight.flightaware isa StructArray{FlightData{Float64}} && isempty(mdata64.flight.flightaware)
        @test mdata64.flight.webdata isa StructArray{FlightData{Float64}} && isempty(mdata64.flight.webdata)
        @test mdata64.cloud.tracks isa StructArray{CloudData{Float64}} && length(mdata64.cloud.tracks) == 3
        @test mdata64.sat.granules isa StructArray{SatData{Float64}} && length(mdata64.sat.granules) == 7

        @test data == dset
        @test data.trackdata.flight.volpe isa StructArray{FlightData{Float32}} &&
            length(data.trackdata.flight.volpe) == 3
        @test data.trackdata.flight.flightaware isa StructArray{FlightData{Float32}} &&
            isempty(data.trackdata.flight.flightaware)
        @test data.trackdata.flight.webdata isa StructArray{FlightData{Float32}} &&
            isempty(data.trackdata.flight.webdata)
        @test data.trackdata.cloud.tracks isa StructArray{CloudData{Float32}} &&
            length(data.trackdata.cloud.tracks) == 3
        @test data.trackdata.sat.granules isa StructArray{SatData{Float32}} &&
            length(data.trackdata.sat.granules) == 7
        @test data.intersection.flight isa XData{Float32} && size(data.intersection.flight.data) == (5, 8) &&
            size(data.intersection.flight.observations) == (5, 4) && size(data.intersection.flight.accuracy) == (5, 6)
        @test data.intersection.cloud isa XData{Float32} && size(data.intersection.cloud.data) == (2, 8) &&
            size(data.intersection.cloud.observations) == (2, 4) && size(data.intersection.cloud.accuracy) == (2, 6)

        @test data64 isa Data{Float64}
        @test data64.trackdata.flight.volpe isa StructArray{FlightData{Float64}} &&
            length(data64.trackdata.flight.volpe) == 3
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
        @test_logs min_level = Debug match_mode = :any (
            :debug, "primary track ID: 42"
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
