"""
    function intersection(
      trackdata::T where T<:Union{Vector{FlightTrack},Vector{CloudTrack}},
      dbmetadata::DBMetadata,
      sat::SatData,
      savesecondsattype::Bool=false,
      maxtimediff::Int=30,
      primspan::Int=0,
      secspan::Int=15,
      lidarrange::Tuple{Real,Real}=(15_000,-Inf),
      stepwidth::Real=0.01,
      Xradius::Real=20_000,
      expdist::Real=Inf,
      Float::DataType=Float32,
      remarks=nothing
    )

Instantiate `Intersection` from either flight or cloud `trackdata` and `sat` tracks
as defined by the kwargs of `Intersection`. The `dbmetadata` is mainly needed for
defining IDs and the source of the track data.
"""
function intersection(
  trackdata::T where T<:Union{<:Vector{<:FlightData},<:Vector{<:CloudTrack}},
  dbmetadata::DBMetadata,
  sat::SatData,
  savesecondsattype::Bool=false,
  maxtimediff::Int=30,
  primspan::Int=0,
  secspan::Int=15,
  lidarrange::Tuple{Real,Real}=(15_000,-Inf),
  stepwidth::Real=0.01,
  Xradius::Real=20_000,
  expdist::Real=Inf,
  Float::DataType=Float32,
  remarks=nothing
)
  # Initialise DataFrames with Intersection data and monitor start time
  tstart = Dates.now()
  Xdata = DataFrame(id=String[], lat=Float[], lon=Float[],
    alt=Float[], tdiff=Dates.CompoundPeriod[], tflight = DateTime[],
    tsat = DateTime[], feature = Union{Missing,Symbol}[])
  tracked = DataFrame(id=String[], flight=FlightTrack[], CPro=CPro[], CLay=CLay[])
  accuracy = DataFrame(id=String[], intersection=Float[], flightcoord=Float[],
    satcoord=Float[], flighttime=Dates.CompoundPeriod[], sattime=Dates.CompoundPeriod[])
  # Get lidar altitude levels
  lidarprofile = get_lidarheights(lidarrange, Float)
  # New MATLAB session
  ms = mat.MSession()
  # Loop over data from different datasets and interpolate track data and time, throw error on failure
  prog = pm.Progress(length(trackdata), "find intersections...")
  for track in trackdata
    # Get database source and track ID for labelling data/issues
    dataset = track isa FlightTrack ? track.metadata.source : "CloudTrack"
    ID = track isa FlightTrack ? track.metadata.dbID : track.metadata.ID
    try
      # Find sat tracks in the vicinity of flight tracks, where intersections are possible
      overlap = findoverlap(track, sat, maxtimediff, ID)
      if isempty(overlap)
        pm.next!(prog, showvalues = [(:hits, length(Xdata.id)),
          (:featured, length(Xdata.id[.!ismissing.(Xdata.feature) .&
          (Xdata.feature .≠ :no_signal) .& (Xdata.feature .≠ :clear)]))])
        continue
      end
      # Interpolate trajectories with PCHIP method
      primtracks = interpolate_trackdata(track)
      sectracks = interpolate_satdata(sat, overlap, track.metadata.useLON)
      # Calculate intersections and store data and metadata in DataFrames
      currdata, currtrack, curraccuracy = find_intersections(ms, track,
        primtracks, dbmetadata.altmin, sat, sectracks, dataset, ID, maxtimediff,
        stepwidth, Xradius, lidarprofile, lidarrange, primspan, secspan,
        expdist, savesecondsattype, Float)
      append!(Xdata, currdata); append!(tracked, currtrack)
      append!(accuracy, curraccuracy)
    catch err
      @debug begin
        @show ID
        rethrow(err)
      end
      # Issue warning on failure of interpolating track or time data
      @warn("Track data and/or time could not be interpolated. Data ignored.",
        dataset, ID)
    end
    # Monitor progress for progress bar
    pm.next!(prog, showvalues = [(:hits, length(Xdata.id)),
      (:featured, length(Xdata.id[.!ismissing.(Xdata.feature) .&
      (Xdata.feature .≠ :no_signal) .& (Xdata.feature .≠ :clear)]))])
  end #loop over flights
  pm.finish!(prog)
  # Close MATLAB session after looping over all data
  mat.close(ms)
  # Calculate load time
  tend = Dates.now()
  tc = tz.ZonedDateTime(tend, tz.localzone())
  loadtime = Dates.canonicalize(Dates.CompoundPeriod(tend - tstart))
  # Return Intersections after completion
  @info string("Intersection data ($(length(Xdata[!,1])) matches) loaded in ",
    "$(join(loadtime.periods[1:min(2,length(loadtime.periods))], ", ")) to",
    "\n▪ data\n▪ tracked\n▪ accuracy\n▪ metadata")
  Intersection(Xdata, tracked, accuracy, XMetadata(maxtimediff, stepwidth, Xradius,
      expdist, lidarrange, lidarprofile, sat.metadata.type, sat.metadata.date,
    dbmetadata.altmin, dbmetadata.date, tc, loadtime, remarks))
end #function intersection


"""
    function find_intersections(
      ms::mat.MSession,
      track::Union{FlightTrack,CloudTrack},
      primtracks::Vector,
      altmin::Real,
      sat::SatData,
      sectracks::Vector,
      dataset::AbstractString,
      trackID::Union{Missing,AbstractString},
      maxtimediff::Int,
      stepwidth::Real,
      Xradius::Real,
      lidarprofile::NamedTuple,
      lidarrange::Tuple{Real,Real},
      primspan::Int,
      secspan::Int,
      expdist::Real,
      savesecondsattype::Bool,
      Float::DataType=Float32
    )

Using the interpolated flight or cloud `primtracks` and satellite `sectracks`,
add new spatial and temporal coordinates of all intersections of the current primary
track with satellite secondary track to `Xdata`, if the overpass of the aircraft/cloud
and the satellite at the intersection is within `maxtimediff` minutes and above
`altmin` in meters.
Additionally, save the measured flight/cloud `track` and `sat` track data near the intersection
(±`primspan`/±`secspan` datapoints of the intersection) in a `tracked` DataFrame
and information about the `accuracy` in another DataFrame.
MATLAB session `ms` is used to retrieve `CPro` and `CLay` satellite data.

When `savesecondsattype` is set to true, the additional satellite data type not
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
the identification of assicated data in the different dataframes of `Intersection`
and as flags in error messages.
"""
function find_intersections(
  ms::mat.MSession,
  track::Union{FlightTrack,CloudTrack},
  primtracks::Vector,
  altmin::Real,
  sat::SatData,
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
  savesecondsattype::Bool,
  Float::DataType=Float32
)
  # Initialise DataFrames for current flight
  Xdata = DataFrame(id=String[], lat=AbstractFloat[], lon=AbstractFloat[],
    alt=AbstractFloat[], tdiff=Dates.CompoundPeriod[], tflight = DateTime[],
    tsat = DateTime[], feature = Union{Missing,Symbol}[])
  tracked = DataFrame(id=String[], flight=FlightTrack[], CPro=CPro[], CLay=CLay[])
  accuracy = DataFrame(id=String[], intersection=AbstractFloat[], flightcoord=AbstractFloat[],
    satcoord=AbstractFloat[], flighttime=Dates.CompoundPeriod[], sattime=Dates.CompoundPeriod[])
  counter = 1 # for intersections within the same flight used in the id
  NA = Float(NaN) # set precision of NaNs according to Float

  # Loop over sat and flight tracks
  for st in sectracks, (n, pt) in enumerate(primtracks)
    # Continue only for sufficient overlap between flight/sat data
    pt.min < st.max && pt.max > st.min || continue
    # Find intersection coordinates
    Xf, Xs = findXcoords(pt, st, stepwidth, track.metadata.useLON, Float)
    for i = 1:length(Xf)
      id = string(dataset,-,trackID,-,counter)
      # Get precision of Intersection
      dx = dist.haversine(Xf[i], Xs[i], earthradius(Xf[i][1]))
      # Determine time difference between aircraf/satellite at intersection
      # Use only the current flight segment for the time interpolation
      # from the flex data in the flight metadata, which coincides with the flighttracks
      # For each flight segment between flex points, one interpolated flight-track exists
      tmf = interpolate_time(track.data[track.metadata.flex[n].range,:], Xf[i])
      tms = interpolate_time(sat.data, Xs[i])
      dt = Dates.canonicalize(Dates.CompoundPeriod(tms-tmf))
      # Skip intersections that exceed allowed time difference
      abs(tmf - tms) < Dates.Minute(maxtimediff) || continue
      # Extract the DataFrame rows of the sat/flight data near the intersection
      Xflight, ift = track isa FlightTrack ? (get_flightdata(track, Xf[i], primspan)) :
        (FlightTrack{Float}(), 0)
      alt = track isa FlightTrack ? Xflight.data.alt[ift] : NA
      cpro, clay, feature, ist = get_satdata(ms, sat, Xs[i], secspan, alt,
        altmin, trackID, lidarprofile, lidarrange, savesecondsattype, Float)
      Xsat = sat.metadata.type == :CPro ? cpro.data : clay.data
      # Calculate accuracies
      fxmeas = track isa FlightTrack ? dist.haversine(Xf[i],(Xflight.data.lat[ift],
        Xflight.data.lon[ift]), earthradius(Xf[i][1])) : NA
      ftmeas = track isa FlightTrack ? Dates.canonicalize(Dates.CompoundPeriod(tmf -
        Xflight.data.time[ift])) : Dates.canonicalize(Dates.CompoundPeriod())
      sxmeas = dist.haversine(Xs[i], (Xsat.lat[ist], Xsat.lon[ist]), earthradius(Xs[i][1]))
      stmeas = Dates.canonicalize(Dates.CompoundPeriod(tms - Xsat.time[ist]))
      # Exclude data with long distances to nearest flight measurement
      if fxmeas > expdist || sxmeas > expdist
        @info("maximum distance of intersection to next track point exceeded; data excluded",
          trackID)
        continue
      end
      # Save intersection data
      counter = addX!(Xdata, tracked, accuracy, counter, Xf[i], id, dx, dt, Xradius, Xflight,
        cpro, clay, tmf, tms, ift, feature, fxmeas, ftmeas, sxmeas, stmeas, NA)
    end #loop over intersections of current flight
  end #loop over flight and sat tracks
  # Return intersection data of current flight
  return Xdata, tracked, accuracy
end #function find_intersections


"""
    findoverlap(flight::FlightTrack, sat::SatDB, maxtimediff::Int) -> Vector{UnitRange}

From the data of the current `flight` and the `sat` data, calculate the data ranges
in the sat data that are in the vicinity of the flight track (min/max of lat/lon).
Consider only satellite data of ± `maxtimediff` minutes before the start and after
the end of the flight.
"""
function findoverlap(track::T where T<:PrimarySet, sat::SatData,
  maxtimediff::Int, ID::Union{Missing,Int,String})

  # Initialise
  overlap = NamedTuple{(:range, :min, :max),Tuple{UnitRange, Real, Real}}[]
  ## Retrieve sat data in the range ±maxtimediff minutes before and after the flight
  # Set time span
  t1 = findfirst(sat.data.time .≥ track.data.time[1] - Dates.Minute(maxtimediff))
  t2 = findlast(sat.data.time .≤ track.data.time[end] + Dates.Minute(maxtimediff))
  # return empty ranges, if no complete overlap is found
  if isnothing(t1) || isnothing(t2)
    @warn string("no sufficient satellite data for time index ",
      "$(track.data.time[1] - Dates.Minute(maxtimediff))...",
      "$(track.data.time[end] + Dates.Minute(maxtimediff))")
    return overlap
  elseif length(t1:t2) ≤ 1
    @warn string("no sufficient overlap between primary/secondary trajectory for track ",
      "$(ID) at $(sat.data.time[t1]) ... $(sat.data.time[t2])")
    return overlap
  end

  ## Find overlaps in flight and sat data
  satoverlap = (track.metadata.area.latmin .≤ sat.data.lat[t1:t2] .≤ track.metadata.area.latmax) .&
    ((track.metadata.area.elonmin .≤ sat.data.lon[t1:t2] .≤ track.metadata.area.elonmax) .|
    (track.metadata.area.wlonmin .≤ sat.data.lon[t1:t2] .≤ track.metadata.area.wlonmax))
  # Convert boolean vector of satellite overlapping data into ranges
  r = false # flag, whether index is part of a current range
  ind = 0   # index in the data array, when looping over data points
  for i = 1:length(satoverlap)
    if satoverlap[i] && !r #First data point of a range found
      r = true # flag as part of a range
      ind = i  # save start index
    elseif r && !satoverlap[i] && length(t1+ind-1:t1+i-2) > 1 # first index of non-overlapping data found
      r = false # flag as non-overlapping data
      # Define current track segment from saved first index to last index
      seg = t1+ind-1:t1+i-2
      # Find flex points at poles in sat tracks
      xdata = track.metadata.useLON ? sat.data.lon[seg] : sat.data.lat[seg]
      satsegments = findflex(xdata)
      # Save current range split into segments with monotonic latitude values
      for s in satsegments
        length(s.range) > 1 && push!(overlap, (range = seg[s.range], min = s.min, max = s.max))
      end
    end
  end
  # Return tuple with sat ranges, and the type of sat data that was used for the calculations
  return overlap
end#function findoverlap


"""
    interpolate_trackdata(flight::FlightTrack)

Using the `flight` data, construct a PCHIP polynomial and return it together with
the x data range (`min`/`max` values).
"""
function interpolate_trackdata(track::Union{FlightTrack,CloudTrack})

  # Define x and y data based on useLON
  x, y = track.metadata.useLON ?
    (track.data.lon, track.data.lat) : (track.data.lat, track.data.lon)
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
    function interpolate_satdata(
      sat::SatData,
      overlap::Vector{NamedTuple{(:range, :min, :max),Tuple{UnitRange, Real, Real}}},
      useLON::Bool
    ) -> Vector{Any}

Using the `sat` data and the stored `overlap` ranges, construct a PCHIP polynomial
and define the range of the x data (`min`/`max` values). X data is defined by the
prevailing flight direction from the `useLON` flag.
"""
function interpolate_satdata(
  sat::SatData,
  overlap::Vector{NamedTuple{(:range, :min, :max),Tuple{UnitRange, Real, Real}}},
  useLON::Bool
)
  # Define x and y data based on useLON
  x, y = useLON ?
    (sat.data.lon, sat.data.lat) : (sat.data.lat, sat.data.lon)

  # Interpolate satellite tracks and flight times for all segments of interest
  idata = []
  # Loop over satellite data
  for r in overlap
    # Interpolate track data with PCHIP and save all interpolated segments together with min/max
    ps = pchip(x[r.range], y[r.range])
    first(x[r.range]) < last(x[r.range]) ?
      push!(idata, (track = interpolate(ps), min = first(x[r.range]), max = last(x[r.range]))) :
      push!(idata, (track = interpolate(ps), min = last(x[r.range]), max = first(x[r.range])))
  end #loop over sat ranges

  # Return a vector with interpolation functions for each dataset
  return idata
end #function interpolate_satdata


"""
    function findXcoords(
      flight::NamedTuple,
      sat::NamedTuple,
      stepwidth::Real,
      useLON::Bool,
      Float::DataType=Float32
    ) -> Xf::, Xs::Tuple{Float,Float}[]

From the interpolated `flight` and `satellite` tracks (PCHIP polynomial), construct
a function `flight - sat` using the PCHIP method on interpolated data points constructed
from both PCHIP polynomials in the common track range with the defined `stepwidth`.
X data is defined by the flag `useLON` (true for longitude as x data) depending on
the prevailing flight direction.

Find all intersections between both tracks by finding all roots in the defined
function using the `IntervalRootFinding` package.

Floating point numbers are of the precision set by `Float`, by default `Float32`.
"""
function findXcoords(
  flight::NamedTuple,
  sat::NamedTuple,
  stepwidth::Real,
  useLON::Bool,
  Float::DataType=Float32
)
  # Define common interpolated x and y data
  xstart, xend = max(flight.min, sat.min), min(flight.max, sat.max)
  xdata = xstart < xend ? collect(Float, xstart:stepwidth:xend) : collect(Float, xend:stepwidth:xstart)
  ydata = flight.track(xdata) .- sat.track(xdata)
  length(xdata) > 1 || return Tuple{Float,Float}[], Tuple{Float,Float}[]

  # Define function to find minimum distance between both tracks
  coorddist(x) = interpolate(pchip(xdata, ydata), x)

  # Find minimum distance by solving flight track - sat track = 0
  rts = root.roots(coorddist, xdata[1] .. xdata[end], root.Krawczyk)
  X = Float.(root.mid.(root.interval.(rts)))

  # Return Vector with coordinate pairs
  Xf, Xs = Tuple{Float,Float}[], Tuple{Float,Float}[]
  for x in X
    isnan(x) && continue
    if useLON
      push!(Xf, (flight.track(x), x))
      push!(Xs, (sat.track(x), x))
    else
      push!(Xf, (x, flight.track(x)))
      push!(Xs, (x, sat.track(x)))
    end
  end

  return Xf, Xs
end #function findXcoords
