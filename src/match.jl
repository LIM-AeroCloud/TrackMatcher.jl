"""
    find_intersections(
      ms::mat.MSession,
      flight::FlightData,
      flighttracks::Vector,
      altmin::Real,
      sat::SatData,
      sattracks::Vector,
      overlap::Vector{UnitRange},
      maxtimediff::Int,
      dmin::Real,
      Xradius::Real,
      lidarprofile::NamedTuple,
      lidarrange::Tuple{Real,Real},
      flightspan::Int,
      satspan::Int,
      savesecondsattype::Bool,
      Float::DataType=Float32
    ) -> Xdata::DataFrame, track::DataFrame, accuracy::DataFrame

Using interpolated `flighttracks` and `sattracks`,
add new spatial and temporal coordinates of all intersections of the current flight
with satellite tracks to `Xdata`, if the overpass of the aircraft and the satellite
at the intersection is within `maxtimediff` minutes and above `altmin` in meters.
Additionally, save the measured `flight` and `sat` `tracked` data near the intersection
(±`flightspan`/±`satspan` datapoints of the intersection) and information about the
`accuracy`.

When `savesecondsattype` is set to true, the additional satellite data type not
used to derive the intersections from the `SatData` is stored as well in `Intersection`.
Satellite column data is stored over the `lidarrange` and defined by the `lidarprofile`.

The algorithm finds intersections, by finding the minimum distance between flight and
sat track points in the overlap region of both tracks. To be an intersection,
the minimum distance of the interpolated track points must be below a threshold `dmin`.
For duplicate intersection finds within an `Xradius`, only the one with the highest
accuracy (lowest `dmin`) is saved.

All floating point numbers are saved with single precision unless otherwise
specified by `Float`.
"""
function find_intersections(
  ms::mat.MSession,
  flight::FlightData,
  flighttracks::Vector,
  altmin::Real,
  sat::SatData,
  sattracks::Vector,
  maxtimediff::Int,
  stepwidth::Real,
  # Xradius::Real,
  # epsilon::Real,
  # tolerance::Real,
  lidarprofile::NamedTuple,
  lidarrange::Tuple{Real,Real},
  flightspan::Int,
  satspan::Int,
  expdist::Real,
  savesecondsattype::Bool,
  Float::DataType=Float32
)
  # Initialise DataFrames for current flight
  Xdata = DataFrame(id=String[], lat=AbstractFloat[], lon=AbstractFloat[],
    alt=AbstractFloat[], tdiff=Dates.CompoundPeriod[], tflight = DateTime[],
    tsat = DateTime[], feature = Union{Missing,Symbol}[])
  track = DataFrame(id=String[], flight=FlightData[], CPro=CPro[], CLay=CLay[])
  accuracy = DataFrame(id=String[], intersection=AbstractFloat[], flightcoord=AbstractFloat[],
    satcoord=AbstractFloat[], flighttime=Dates.CompoundPeriod[], sattime=Dates.CompoundPeriod[])
  counter = 0 # for intersections within the same flight used in the id
  # Loop over sat and flight tracks
  for (i, st) in enumerate(sattracks), ft in flighttracks
    # Continue only for sufficient overlap between flight/sat data
    ft.min < st.max && ft.max > st.min || continue
    # Find intersection coordinates
    Xs, Xf = findXcoords(ft, st, stepwidth, flight.metadata.useLON, Float)

    for i = 1:length(Xf)
      # Get ID for current intersection
      id = string(flight.metadata.source,-,flight.metadata.dbID,-,i)
      # Get precision of Intersection
      dx = dist.haversine(Xf[i], Xs[i], earthradius(Xf[i][1]))
      # Determine time difference between aircraf/satellite at intersection
      tmf = interpolate_time(flight.data, Xf[i])
      tms = interpolate_time(sat.data, Xs[i])
      dt = Dates.canonicalize(Dates.CompoundPeriod(tms-tmf))
      # Skip intersections that exceed allowed time difference
      abs(tmf - tms) < Dates.Minute(maxtimediff) || continue
      # Extract the DataFrame rows of the sat/flight data near the intersection
      Xflight, ift = get_flightdata(flight, Xf[i], flightspan)
      cpro, clay, feature, ist = get_satdata(ms, sat, Xs[i], satspan, Xflight.data.alt[ift],
        altmin, Xflight.metadata.dbID, lidarprofile, lidarrange, savesecondsattype, Float)
      Xsat = sat.metadata.type == :CPro ? cpro.data : clay.data
      # Calculate accuracies
      fxmeas = dist.haversine(Xf[i],(Xflight.data.lat[ift], Xflight.data.lon[ift]),
        earthradius(Xf[i][1]))
      ftmeas = Dates.canonicalize(Dates.CompoundPeriod(tmf - Xflight.data.time[ift]))
      sxmeas = dist.haversine(Xs[i], (Xsat.lat[ist], Xsat.lon[ist]), earthradius(Xs[i][1]))
      stmeas = Dates.canonicalize(Dates.CompoundPeriod(tms - Xsat.time[ist]))
      # Exclude data with long distances to nearest flight measurement
      if fxmeas > expdist || sxmeas > expdist
        @info("maximum distance of intersection to next track point exceeded; data excluded",
          flight.metadata.dbID)
        continue
      end
      # Save intersection data
      push!(Xdata, (id=id, lat=Xf[i][1], lon=Xf[i][2], alt=Xflight.data.alt[ift],
        tdiff=dt, tflight = tmf, tsat = tms, feature=feature))
      push!(track, (id = id, flight = Xflight, CPro = cpro, CLay = clay))
      # Save accuracies
      push!(accuracy, (id=id, intersection=dx, flightcoord=fxmeas,
        satcoord=sxmeas, flighttime=ftmeas, sattime=stmeas))
    end #loop over intersections of current flight
  end #loop over flight and sat tracks

  # Return intersection data of current flight
  return Xdata, track, accuracy
end #function find_intersections


"""
    findoverlap(flight::FlightData, sat::SatDB, maxtimediff::Int) -> Vector{UnitRange}

From the data of the current `flight` and the `sat` data, calculate the data ranges
in the sat data that are in the vicinity of the flight track (min/max of lat/lon).
Consider only satellite data of ± `maxtimediff` minutes before the start and after
the end of the flight.
"""
function findoverlap(flight::FlightData, sat::SatData, maxtimediff::Int)

  # Initialise
  overlap = NamedTuple{(:range, :min, :max),Tuple{UnitRange, Real, Real}}[]
  ## Retrieve sat data in the range ±maxtimediff minutes before and after the flight
  # Set time span
  t1 = findfirst(sat.data.time .≥ flight.data.time[1] - Dates.Minute(maxtimediff))
  t2 = findlast(sat.data.time .≤ flight.data.time[end] + Dates.Minute(maxtimediff))
  # return empty ranges, if no complete overlap is found
  if isnothing(t1) || isnothing(t2)
    @warn string("no sufficient satellite data for time index ",
      "$(flight.data.time[1] - Dates.Minute(maxtimediff))...",
      "$(flight.data.time[end] + Dates.Minute(maxtimediff))")
    return overlap
  elseif length(t1:t2) ≤ 1
    @warn string("no sufficient overlap between flight/satellite data for flight ",
      "$(flight.metadata.dbID) at ",
      "$(sat.data.time[t1]) ... $(sat.data.time[t2])")
    return overlap
  end

  ## Find overlaps in flight and sat data
  satoverlap = (flight.metadata.area.latmin .≤ sat.data.lat[t1:t2] .≤ flight.metadata.area.latmax) .&
    ((flight.metadata.area.elonmin .≤ sat.data.lon[t1:t2] .≤ flight.metadata.area.elonmax) .|
    (flight.metadata.area.wlonmin .≤ sat.data.lon[t1:t2] .≤ flight.metadata.area.wlonmax))
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
      xdata = flight.metadata.useLON ? sat.data.lon[seg] : sat.data.lat[seg]
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
    interpolate_flightdata(flight::FlightData)

Using the `flight` data, interpolate the data with the pchip method
using the defined `stepwidth` for interpolation.
"""
function interpolate_flightdata(flight::FlightData)

  # Define x and y data based on useLON
  x, y = flight.metadata.useLON ?
    (flight.data.lon, flight.data.lat) : (flight.data.lat, flight.data.lon)
  # Interpolate flight tracks and tims for all segments
  idata = []
  for f in flight.metadata.flex
    # Interpolate track data with PCHIP
    pf = pchip(x[f.range],y[f.range])
    # Save the interpolating polynomial to a vector
    first(x[f.range]) < last(x[f.range]) ?
      push!(idata, (track = interpolate(pf), min = first(x[f.range]), max = last(x[f.range]))) :
      push!(idata, (track = interpolate(pf), min = last(x[f.range]), max = first(x[f.range])))
  end

  # Return the interplated data
  return idata
end #function interpolate_flightdata


"""
    interpolate_satdata(sat::SatData, overlap::Vector{UnitRange})
      -> Vector{Any}

Using the `sat` data and the stored `overlap` ranges, interpolate the data
with the pchip method. Use the metadata in `flight` for error reports.
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


# coordinaterange(lat::Vector{<:AbstractFloat}, lon::Vector{<:AbstractFloat}) = lat[1] < lat[end] ?
#   (lat, lon) : (reverse(lat), reverse(lon))


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
  rts = root.roots(coorddist, xdata[1] .. xdata[end])
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
end

""" Generate a function to interpolate x vectors from the current PCHIP struct """
# interpolate(pc::Polynomial{<:AbstractFloat}) = x::Union{<:AbstractFloat,Vector{<:AbstractFloat}} -> interpolate(pc, x)
