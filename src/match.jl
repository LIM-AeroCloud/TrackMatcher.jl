### functions related to track-matching/finding intersections

"""
    find_intersections(
        track::PrimaryTrack,
        primtracks::Vector,
        altmin::Real,
        sat::SatSet,
        sectracks::Vector,
        dataset::AbstractString,
        trackID::Union{Missing,Int,AbstractString},
        maxtimediff::Int,
        stepwidth::Real,
        Xradius::Real,
        lidarprofile::NamedTuple,
        lidarrange::Tuple{Real,Real},
        primspan::Int,
        secspan::Int,
        expdist::Real,
        savedir::Union{String,Bool},
        savesecondsattype::Bool,
        Float::DataType=Float32
    ) -> Tuple{DataFrame,DataFrame,DataFrame}

Using the interpolated flight or cloud `primtracks` and satellite `sectracks`,
add new spatial and temporal coordinates of all intersections of the current primary
track with the satellite secondary track to `Xdata`, if the overpass of the aircraft/cloud
and the satellite at the intersection is within `maxtimediff` minutes and above `altmin`
 in meters.
Additionally, save the measured flight/cloud `track` and `sat` track data near the intersection
(±`primspan`/±`secspan` datapoints of the intersection) in a `observations` DataFrame
and information about the `accuracy` in another DataFrame.

When `savesecondsattype` is set to `true`, the additional satellite data type not
used to derive the intersections from the `SatData` is stored as well in `Intersection`.
Satellite column data is stored over the `lidarrange` as defined by the `lidarprofile`.

The algorithm finds intersections, by finding roots of a distance function
`d(x) = primtrack(x) - sectrack(x)` with the `IntervalRootFinding` package.
Therefore, flight/cloud `primtracks` and satellite `sectracks` are interpolated
using the defined `stepwidth`. For duplicate intersection finds within an `Xradius`,
only the one with the highest accuracy (lowest `accuracy.intersection`) is saved
unless the distance to the nearest measured track point exceeds the maximum
`expdist` threshold of either track.

All floating point numbers are saved with single precision unless otherwise
specified by `Float`. The `dataset` and `trackID` identifiers are used for
the identification of associated data in the different dataframes of `Intersection`
and as flags in error messages.
"""
function find_intersections(
    track::PrimaryTrack,
    primtracks::Vector,
    altmin::Real,
    sat::SatSet,
    sectracks::Vector,
    dataset::AbstractString,
    trackID::Union{Missing,Int32,AbstractString},
    maxtimediff::Int,
    stepwidth::Real,
    Xradius::Real,
    lidarprofile::NamedTuple,
    lidarrange::Tuple{Real,Real},
    primspan::Int,
    secspan::Int,
    expdist::Real,
    savedir::Union{String,Bool},
    savesecondsattype::Bool,
    Float::DataType=Float32
)
    # Initialise DataFrames for current flight
    Xdata = DataFrame(id=String[], lat=Float[], lon=Float[],
        alt=Union{Missing,Float}[], tdiff=Dates.CompoundPeriod[], tprim = DateTime[],
        tsec = DateTime[], atmos_state = Enum{UInt16}[])
    observations = DataFrame(id=String[], primary=PrimaryTrack[], CPro=CPro[], CLay=CLay[])
    accuracy = DataFrame(id=String[], intersection=Float[], primdist=Float[],
        secdist=Float[], primtime=Dates.CompoundPeriod[], sectime=Dates.CompoundPeriod[])
    counter = 1 # for intersections within the same flight used in the id

     # Loop over sat and flight tracks
    for st in sectracks, (n, pt) in enumerate(primtracks)
        # Continue only for sufficient overlap between flight/sat data
        pt.min < st.max && pt.max > st.min || continue
        # Find intersection coordinates
        Xp, Xs = findXcoords(pt, st, stepwidth, track.metadata.use_lon, Float)
        for i = 1:length(Xp)
        # Get precision of Intersection
        dx = dist.haversine(Xp[i], Xs[i], earthradius(Xp[i][1])) |> Float
        # Determine time difference between aircraf/satellite at intersection
        # Use only the current flight segment for the time interpolation
        # from the flex data in the flight metadata, which coincides with the flighttracks
        # For each flight segment between flex points, one interpolated flight-track exists
        tmf = interpolate_time(track, Xp[i], n)
        tms, obsindex = interpolate_time(sat, st.timeindex, Xs[i])
        dt = Dates.canonicalize(Dates.CompoundPeriod(tms-tmf))
        # Skip intersections that exceed allowed time difference
        abs(tmf - tms) < Dates.Minute(maxtimediff) || continue
        # Add intersection and nearby measurements
        id = string(dataset,-,trackID,-,counter)
        counter = track isa FlightTrack ?
            add_intersections!(Xdata, observations, accuracy, track, sat, Xp[i], Xs[i],
            obsindex, counter, id, dx, dt, tmf, tms, primspan, secspan, altmin, trackID,
            Xradius, expdist, lidarprofile, lidarrange, savedir, savesecondsattype) :
            add_intersections!(Xdata, observations, accuracy, sat, Xp[i], Xs[i],
            obsindex, counter, id, dx, dt, tmf, tms, secspan, altmin, trackID,
            Xradius, expdist, lidarprofile, lidarrange, savedir, savesecondsattype)
        end #loop over intersections of current flight
    end #loop over flight and sat tracks
    # Return intersection data of current flight
    return Xdata, observations, accuracy
end #function find_intersections


"""
    function findoverlap(
        primtrack::PrimaryTrack,
        sectrack::SatSet,
        maxtimediff::Int,
        atol::Real=0.1
    ) -> Vector{Vector{DataFrame}}

Find all track segments in `sectrack` that are within the time frame ± `maxtimediff` of the
`primtrack` and within the spatial bounding box of `primtrack` and an absolute tolerance `atol`.
"""
function findoverlap(
    primtrack::PrimaryTrack,
    sectrack::SatSet,
    maxtimediff::Int,
    atol::Real=0.1
)
    # Select granules within flight time frame ± tolerance
    t1 = findlast(primtrack.metadata.date.start .≥ sectrack.metadata.granules.tstart)
    t2 = findfirst(primtrack.metadata.date.stop .≤ sectrack.metadata.granules.tstop)
        if isnothing(t1) || isnothing(t2)
        @warn string("no sufficient satellite data for time index ",
            "$(primtrack.time[1] - Dates.Minute(maxtimediff))...",
            "$(primtrack.time[end] + Dates.Minute(maxtimediff))")
        return DataFrame[], 0:-1
        end
    dt = t1 < t2 ? (t1:t2) : (t2:t1)
    # Filter granules without an overlapping area
    inarea = [!(granule.elonmin - atol > primtrack.metadata.area.elonmax ||
        granule.elonmax + atol < primtrack.metadata.area.elonmin ||
        granule.wlonmin - atol > primtrack.metadata.area.wlonmax ||
        granule.wlonmax + atol < primtrack.metadata.area.wlonmin)
        for granule in df.eachrow(sectrack.metadata.granules[t1:t2,:])]
    segments = [DataFrame(time = granule.time, lat = granule.lat, lon = granule.lon)
        for granule in sectrack.granules[dt][inarea]]
    # filter granules for segments within a bounding box of the flight track
    filter(df -> size(df, 1) > 1, filter.(withinlimits(primtrack.metadata.area, atol), segments)), dt
end #function findoverlap


"""
    interpolate_trackdata(track::PrimaryTrack)

Using the `track` data, construct a PCHIP polynomial and return it together with
the x data range (`min`/`max` values).
"""
function interpolate_trackdata(track::PrimaryTrack)

    # Define x and y data based on use_lon
    x, y = track.metadata.use_lon ?
        (track.lon, track.lat) : (track.lat, track.lon)
    # Interpolate tracks and times for all segments
    idata = []
    for f in track.metadata.flex
        # Interpolate track data with PCHIP
        pf = pchip(x[f.range],y[f.range])
        # Save the interpolating polynomial to a vector
        first(x[f.range]) < last(x[f.range]) ?
        push!(idata, (track = interpolate(pf), min = first(x[f.range]), max = last(x[f.range]))) :
        push!(idata, (track = interpolate(pf), min = last(x[f.range]), max = first(x[f.range])))
    end

    # Return the interplated data
    return idata
end #function interpolate_trackdata


"""
    interpolate_satdata(
        trackdata::Vector{DataFrame},
        isat::UnitRange{Int},
        use_lon::Bool
    ) -> Vector{Any}

Interpolate `trackdata` with a PCHIP polynomial and define the range of the x data
(`min`/`max` values). X data is defined by the prevailing flight direction from the
`use_lon` flag. Save the time index `isat` of the satellite track for later retrieval
of the observation data.
"""
function interpolate_satdata(
    trackdata::Vector{DataFrame},
    isat::UnitRange{Int},
    use_lon::Bool
)
    # Define x and y data based on use_lon
    x, y = use_lon ? (:lon, :lat) : (:lat, :lon)

    # Interpolate satellite tracks and flight times for all segments of interest
    idata = []
    # Loop over satellite data
    for segment in trackdata
        # Split track in segments with monotonic x data
        seg = findall((segment[1:end-2, x] .< segment[2:end-1, x] .> segment[3:end, x]) .|
            (segment[1:end-2, x] .> segment[2:end-1, x] .< segment[3:end, x])) .+ 1
        pushfirst!(seg, 1); push!(seg, size(segment, 1))
        # Interpolate track segments with PCHIP and save all interpolated segments together with min/max
        for i = 2:length(seg)
        currseg = segment[seg[i-1]:seg[i], :]
        ps = pchip(currseg[:, x], currseg[:, y])
        x0, x1 = currseg[1, x] < currseg[end, x] ?
            (currseg[1, x], currseg[end, x]) : (currseg[end, x], currseg[1, x])
            push!(idata, (track = interpolate(ps), min = x0, max = x1, timeindex = isat))
        end
    end #loop over sat ranges

    # Return a vector with interpolation functions for each dataset
    return idata
end #function interpolate_satdata


"""
    function findXcoords(
        track::NamedTuple,
        sat::NamedTuple,
        stepwidth::Real,
        use_lon::Bool,
        Float::DataType=Float32
    ) -> Xp::, Xs::Tuple{Float,Float}[]

From the interpolated primary `track` and `sat` track (PCHIP polynomial), construct
a function `track - sat` using the PCHIP method on interpolated data points constructed
from both PCHIP polynomials in the common track range with the defined `stepwidth`.
X data is defined by the flag `use_lon` (true for longitude as x data) depending on
the prevailing primary track direction.

Find all intersections between both tracks by finding all roots in the defined
function using the `IntervalRootFinding` package.

Floating point numbers are of the precision set by `Float`, by default `Float32`.
"""
function findXcoords(
    track::NamedTuple,
    sat::NamedTuple,
    stepwidth::Real,
    use_lon::Bool,
    Float::DataType=Float32
)
    # Define common interpolated x and y data
    xstart, xend = max(track.min, sat.min), min(track.max, sat.max)
    xdata = xstart < xend ? collect(Float, xstart:stepwidth:xend) : collect(Float, xend:stepwidth:xstart)
    ydata = track.track(xdata) .- sat.track(xdata)

    length(xdata) > 1 || return Tuple{Float,Float}[], Tuple{Float,Float}[]

    # Use double precision for intersection finding
    xdata, ydata = Float64.(xdata), Float64.(ydata)
    # Define function to find minimum distance between both tracks
    pc = pchip(xdata, ydata)
    function coorddist(x)
        if x isa intar.Interval
            # Evaluate at interval bounds to get a conservative interval result
            xlo, xhi = intar.inf(x), intar.sup(x)
            ylo, yhi = interpolate(pc, xlo), interpolate(pc, xhi)
            return min(ylo, yhi) .. max(ylo, yhi)
        else
            return interpolate(pc, x)
        end
    end
    # Find minimum distance by solving primary track - sat track = 0
    # Disable automatic differentiation to avoid ForwardDiff with intervals
    rts = root.roots(coorddist, xdata[1] .. xdata[end], derivative=Returns(1 .. 1))
    X = Float.([intar.mid(r.region) for r in rts])

    # Return Vector with coordinate pairs
    Xp, Xs = Tuple{Float,Float}[], Tuple{Float,Float}[]
    for x in X
        isnan(x) && continue
        if use_lon
            push!(Xp, (track.track(x), x))
            push!(Xs, (sat.track(x), x))
        else
            push!(Xp, (x, track.track(x)))
            push!(Xs, (x, sat.track(x)))
        end
    end

    return Xp, Xs
end #function findXcoords
