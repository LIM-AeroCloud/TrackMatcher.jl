# ¡ Needs general test data from init.jl and helper functions from setup.jl

## Test sets

@testset "intersections" begin

    # Run intercept finding routines
    xf_cpro = Intersection(flight, sat_cpro, expdist=80_000) # ℹ reused in constructor testset
    xf64 = Intersection{Float64}(xf_cpro)
    xf64_promoted = Intersection(flight64, sat_cpro, true, expdist=80_000)
    xf64_forced = XData{Float64}(flight, sat_cpro, expdist=80_000)
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
        @test all(abs.(getindex.(primary_coords, 1) .- getindex.(secondary_coords, 1)) .≤ 1e-6)
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
        @test all(id -> !startswith(id, "V-44"), xf_cpro.data.id)
        @test xf64 isa XData{Float64}
        @test xf64 ≈ xf_cpro
        @test xf64_promoted isa XData{Float64}
        @test xdata_matches(xf64_promoted, xf64, [true, true, true], accuracy_atol=1)
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
        test_logger = Test.TestLogger()
        xf_clay = with_logger(test_logger) do
            XData(flight, sat_clay, true, expdist=80_000)
        end
        @test test_logger.logs[1].message == "maximum distance of intersection to next track point exceeded; data excluded"
        @test test_logger.logs[1].kwargs[:trackID] == 42
        @test test_logger.logs[1].level == Logging.Info
        @test test_logger.logs[2].message == "no sufficient satellite data for time index 2012-02-06T02:30:00...2012-02-06T03:39:00"
        @test test_logger.logs[2].level == Logging.Warn
        @test startswith(test_logger.logs[3].message, "Intersection data (4 matches) loaded")
        @test xresults(xf_clay, flight, approx=false)
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
        mock = which(TrackMatcher.interpolate_trackdata, Tuple{FlightTrack})
        try
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
        finally
            Base.delete_method(mock)
        end
        overlap, range = @test_logs (
            :warn, r"no sufficient satellite data"
        ) TrackMatcher.findoverlap(flight.volpe[1], SatSet(), 30, 0.1)
        @test overlap isa Vector{DataFrame} && isempty(overlap)
        @test range ==0:-1
    end
    @testset "exception handling for secondary sat data" begin
        # Force an exception loading the secondary CPro data (primary type CLay)
        TrackMatcher.CPro{Float32}(files::Vector{String}, timeindex::Vector{UnitRange{Int64}},
            lidarprofile::NamedTuple, saveobs::Bool) = throw(ErrorException("forced failure"))
        mock_cpro = which(TrackMatcher.CPro{Float32},
            Tuple{Vector{String},Vector{UnitRange{Int64}},NamedTuple,Bool})
        try
            logger = Test.TestLogger()
            xclay = Logging.with_logger(logger) do
                XData(flight, sat_clay, true, expdist=80_000)
            end
            warns = filter(r -> r.message == "could not load additional profile data", logger.logs)
            @test !isempty(warns)
            @test all(r -> r.kwargs[:trackID] in (Int32(42), Int32(43)), warns)
            @test all(c -> isempty(c.time), xclay.observations.CPro)
        finally
            Base.delete_method(mock_cpro)
        end

        # Force an exception loading the secondary CLay data (primary type CPro)
        TrackMatcher.CLay{Float32}(files::Vector{String}, timeindex::Vector{UnitRange{Int64}},
            lidarrange::Tuple{Real,Real}, altmin::Real) = throw(ErrorException("forced failure"))
        mock_clay = which(TrackMatcher.CLay{Float32},
            Tuple{Vector{String},Vector{UnitRange{Int64}},Tuple{Real,Real},Real})
        try
            logger = Test.TestLogger()
            xcpro = Logging.with_logger(logger) do
                XData(flight, sat_cpro, true, expdist=80_000)
            end
            warns = filter(r -> r.message == "could not load additional layer data", logger.logs)
            @test !isempty(warns)
            @test all(r -> r.kwargs[:trackID] in (Int32(42), Int32(43)), warns)
            @test all(c -> isempty(c.time), xcpro.observations.CLay)
        finally
            Base.delete_method(mock_clay)
        end
    end
end
