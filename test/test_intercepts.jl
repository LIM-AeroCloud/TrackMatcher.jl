# ¡ Needs data from test_flightdata.jl, test_clouddata.jl, and test_satdata.jl to run

## Setup and helper functions
# Debug test data
flight = FlightSet(volpe=joinpath(@__DIR__, "data", "volpe"))
flight64 = FlightSet{Float64}(volpe=joinpath(@__DIR__, "data", "volpe"))
flight_empty = FlightSet()
cloud = CloudSet(joinpath(@__DIR__, "data", "cloud"))

cpro_src = joinpath(@__DIR__, "data", "caliop", "CPro")
clay_src = joinpath(@__DIR__, "data", "caliop", "CLay")
cpro = SatSet(cpro_src, type=:CPro)
clay = SatSet(clay_src, type=:CLay)

# TODO adjust to cloud data
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

# function promoted_xdata_matches(a::XData{T}, b::XData{T}; data_atol::Real=1e-3, accuracy_atol::Real=0.2) where T
#     a.data.id == b.data.id || return false
#     isapprox(a.data.lat, b.data.lat; atol=data_atol, rtol=0.0) || return false
#     isapprox(a.data.lon, b.data.lon; atol=data_atol, rtol=0.0) || return false
#     isapprox(a.data.alt, b.data.alt; atol=data_atol, rtol=0.0) || return false
#     a.data.tdiff == b.data.tdiff || return false
#     a.data.tprim == b.data.tprim || return false
#     a.data.tsec == b.data.tsec || return false
#     a.data.atmos_state == b.data.atmos_state || return false

#     a.accuracy.id == b.accuracy.id || return false
#     isapprox(a.accuracy.intersection, b.accuracy.intersection; atol=accuracy_atol, rtol=0.0) || return false
#     isapprox(a.accuracy.primdist, b.accuracy.primdist; atol=accuracy_atol, rtol=0.0) || return false
#     isapprox(a.accuracy.secdist, b.accuracy.secdist; atol=accuracy_atol, rtol=0.0) || return false
#     a.accuracy.primtime == b.accuracy.primtime || return false
#     a.accuracy.sectime == b.accuracy.sectime || return false

#     a.observations.id == b.observations.id || return false
#     names(a.observations) == names(b.observations) || return false
#     return true
# end

function xresults(intersections::XData{T}, primary::PrimarySet{T}; obs::Vector{Bool}=[true, true, true], approx::Union{Bool,Int}=true) where T
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
    expected = XData{T}(data, observations, accuracy, intersections.metadata) # bug: observations
    results = if approx isa Int
        xdata_matches(intersections, expected, obs;
            data_atol=10.0^-approx, accuracy_atol=10.0^-approx, rtol=0.0)
    elseif approx === true
        xdata_matches(intersections, expected, obs;
            data_atol=1e-3, accuracy_atol=1e-1, rtol=0.0)
    else
        xdata_matches(intersections, expected, obs;
            data_atol=1e-3, accuracy_atol=1e-1, rtol=0.0)
    end
    return results
end

## Test sets

@testset "intersections" begin
    # Run intercept finding routines
    xf_cpro = Intersection(flight, cpro)
    xf_clay = XData(flight, clay, true)
    xf64 = Intersection{Float64}(xf_cpro)
    xf64_promoted = Intersection(flight64, cpro)
    xf64_forced = XData{Float64}(flight, cpro)
    xc_pro = XData(cloud, cpro)
    xc_lay = XData(cloud, clay)
    xf_empty = XData(flight_empty, cpro)
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

        sat = cpro.granules[1]
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
    @testset "exception handling" begin
        # Force an exception in the interpolation path
        TrackMatcher.interpolate_trackdata(::TrackMatcher.PrimaryTrack) =
            throw(ErrorException("forced failure"))
        @test_logs (
            :warn, r"Track data and/or time could not be interpolated"
        ) (
            :info, r"Intersection data \(0 matches\) loaded"
        ) begin
            x = XData(flight, cpro, true)
            @test x isa XData
        end
    end
end
