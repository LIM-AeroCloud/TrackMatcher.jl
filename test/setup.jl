## Helper functions for satellite data

function metadata(meta::SecondaryMetadata)
    expected = (
        start=DateTime(2012, 2, 5, 21, 16, 31, 602),
        stop=DateTime(2012, 2, 6, 3, 5, 33, 437),
        latmin=Float32[-81.82188, -68.64167],
        latmax=Float32[78.177505, 81.821526],
        elonmin=Float32[0.09188683, 20.686209],
        elonmax=Float32[44.556667, 179.92279],
        wlonmin=Float32[-169.5283, -179.86685],
        wlonmax=Float32[-0.17311835, -170.95845]
    )
    return meta.date.start == expected.start && meta.date.stop == expected.stop &&
        all(isapprox.(meta.granules.latmin[1:2], expected.latmin; atol=1e-6)) &&
        all(isapprox.(meta.granules.latmax[1:2], expected.latmax; atol=1e-6)) &&
        all(isapprox.(meta.granules.elonmin[1:2], expected.elonmin; atol=1e-6)) &&
        all(isapprox.(meta.granules.elonmax[1:2], expected.elonmax; atol=1e-6)) &&
        all(isapprox.(meta.granules.wlonmin[1:2], expected.wlonmin; atol=1e-6)) &&
        all(isapprox.(meta.granules.wlonmax[1:2], expected.wlonmax; atol=1e-6))
end


"""
    test_sat_datafiles(files, type, satfiles, expected_type=type) -> Bool

Test function `TrackMatcher.sat_datafiles` to return the `expected_type` from the given `type`
and clean `files` to match `satfiles`.
Returns `true` if both tests pass, `false` otherwise.
"""
function test_sat_datafiles(files, type, satfiles, expected_type=type)::Bool
    success = true
    input_files = copy(files)
    sattype = TrackMatcher.sat_datafiles!(input_files, type)
    success &= sattype == expected_type
    success &= input_files == satfiles
    return success
end


## Helper functions for intercept finding

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
                  [ci, invalid, invalid, invalid] : [ci, clear, clear, clear]
    data = primary isa FlightSet ? DataFrame(
        id=["V-42-1", "V-42-2", "V-42-3", "V-43-1"],
        lat=T[29.305733, 40.152092, 49.438873, 49.438488],
        lon=T[29.0, 32.091255, 35.443886, 35.44373],
        alt=T[11004.804, 11507.114, 11502.542, 11503.152],
        tdiff=Dates.CompoundPeriod[
            Dates.CompoundPeriod(Minute(1), Second(36)),
            Dates.CompoundPeriod(Second(45)),
            Dates.CompoundPeriod(),
            Dates.CompoundPeriod()
        ],
        tprim=[DateTime(2012, 2, 6, 0, 5, 8), DateTime(2012, 2, 6, 0, 2, 58),
            DateTime(2012, 2, 6, 0, 1, 7), DateTime(2012, 2, 6, 0, 1, 7)],
        tsec=[DateTime(2012, 2, 6, 0, 6, 44), DateTime(2012, 2, 6, 0, 3, 43),
            DateTime(2012, 2, 6, 0, 1, 7), DateTime(2012, 2, 6, 0, 1, 7)],
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
        id=["V-42-1", "V-42-2", "V-42-3", "V-43-1"],
        intersection=T[0.0, 0.0, 0.4238313, 0.3452892],
        primdist=T[67541.35, 17558.139, 51154.125, 0.0],
        secdist=T[872.12195, 1785.7188, 479.24576, 440.39847],
        primtime=[Dates.CompoundPeriod(Second(8)), Dates.CompoundPeriod(Second(-2)),
            Dates.CompoundPeriod(Second(7)), Dates.CompoundPeriod()],
        sectime=[Dates.CompoundPeriod(Millisecond(-7)), Dates.CompoundPeriod(Millisecond(-220)),
            Dates.CompoundPeriod(Millisecond(16)), Dates.CompoundPeriod(Millisecond(16))]
    ) : DataFrame(
        id=["C-1-1"],
        intersection=T[34.18339],
        primdist=T[NaN],
        secdist=T[NaN],
        primtime=[Dates.CompoundPeriod()],
        sectime=[Dates.CompoundPeriod()]
    )
    observations = primary isa FlightSet ? DataFrame(
        id=["V-42-1", "V-42-2", "V-42-3", "V-43-1"],
        # Non-empty dummy data for observations
        primary=[flight.volpe[1], flight.volpe[1], flight.volpe[1], flight.volpe[2]],
        CPro=[obsdata[1] for _ in 1:4],
        Clay=[obsdata[2] for _ in 1:4]
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
