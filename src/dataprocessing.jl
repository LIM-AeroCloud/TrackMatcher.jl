### Helper functions for data processing

# Overload abs function from base to get the absolute number of CompoundPeriod
Base.abs(dt::Dates.CompoundPeriod) = dt > Dates.CompoundPeriod(Dates.Millisecond(0)) ? dt : -dt

## Storage of intersection data

"""
    addX!(
      Xdata::DataFrame,
      track::DataFrame,
      accuracy::DataFrame,
      counter::Int,
      Xf::Tuple{T,T},
      id::String,
      dx::T,
      dt::Dates.CompoundPeriod,
      Xradius::Real,
      Xflight::FlightData{T},
      cpro::CPro,
      clay::CLay,
      tmf::DateTime,
      tms::DateTime,
      atmos::Union{Missing,Symbol},
      fxmeas::T,
      ftmeas::Dates.CompoundPeriod,
      sxmeas::T,
      stmeas::Dates.CompoundPeriod,
      alt::Real
    ) where T

Append DataFrames `Xdata`, `track`, and `accuracy` by data `Xf`, `id`,
`dx`, `dt`, `Xflight`, `cpro`, `clay`, `tmf`, `tms`, `atmos`, `fxmeas`,
`ftmeas`, `sxmeas`, `stmeas`, and `alt`. If an intersection already exists within
`Xradius` in `Xdata`, use the more accurate intersection with the lowest
`accuracy.intersection`. Increase the `counter` for new entries.
"""
function addX!(
  Xdata::DataFrame,
  track::DataFrame,
  accuracy::DataFrame,
  counter::Int,
  Xf::Tuple{T,T},
  id::String,
  dx::T,
  dt::Dates.CompoundPeriod,
  Xradius::Real,
  Xflight::FlightData{T},
  cpro::CPro,
  clay::CLay,
  tmf::DateTime,
  tms::DateTime,
  atmos::Union{Missing,Symbol},
  fxmeas::T,
  ftmeas::Dates.CompoundPeriod,
  sxmeas::T,
  stmeas::Dates.CompoundPeriod,
  alt::Real
) where T

  # Loop over previously found intersections
  for i = 1:size(Xdata, 1)
    # Use most accurate intersection, when duplicates are found within Xradius
    # or intersection with least decay between overpass times for equal accuracies
    if dist.haversine(Xf, (Xdata.lat[i], Xdata.lon[i]), earthradius(Xf[1])) ≤ Xradius
      dx ≤ accuracy.intersection[i] || return counter # previous intersection more accurate
      # previous intersection equally accurate, but smaller delay time:
      (dx == accuracy.intersection[i] && abs(dt) > abs(Xdata.tdiff[i])) && return counter

      # Save more accurate duplicate
      Xdata[i, 2:end] = (lat = Xf[1], lon = Xf[2], alt = alt,
        tdiff = dt, tflight = tmf, tsat = tms, atmos_state = atmos)
      track[i, 2:end] = (flight = Xflight, CPro = cpro, CLay = clay)
      accuracy[i, 2:end] = (intersection = dx, flightcoord = fxmeas,
        satcoord = sxmeas, flighttime = ftmeas, sattime = stmeas)
      return counter
    end # duplicate condition based on accuracy
  end #loop over already found intersection
  # Save new intersections that are not identified as duplicates and increase counter
  counter += 1
  push!(Xdata, (id = id, lat = Xf[1], lon = Xf[2], alt = alt,
    tdiff = dt, tflight = tmf, tsat = tms, atmos_state = atmos))
  push!(track, (id = id, flight = Xflight, CPro = cpro, CLay = clay))
  push!(accuracy, (id = id, intersection = dx, flightcoord = fxmeas,
    satcoord = sxmeas, flighttime = ftmeas, sattime = stmeas))

  return counter
end #function addX!


"""
    function add_intersections!(
      ms::mat.MSession,
      Xdata::DataFrame,
      tracked::DataFrame,
      accuracy::DataFrame,
      track::FlightTrack,
      sat::SatData,
      Xf::Tuple{T,T},
      Xs::Tuple{T,T},
      counter::Int,
      id::String,
      dx::T,
      dt::Dates.CompoundPeriod,
      tmf::DateTime,
      tms::DateTime,
      primspan::Int,
      secspan::Int,
      altmin::Real,
      trackID::Union{Missing,Int,AbstractString},
      Xradius::Real,
      expdist::Real,
      lidarprofile::NamedTuple,
      lidarrange::Tuple{Real,Real},
      savesecondsattype::Bool
    ) where T<:AbstractFloat

Add intersection data to the DataFrames `Xdata`, `tracked`, and `accuracy` using
the coordinates calculated from the primary (`Xf`) and secondary (`Xs`) track, the
primary flight `track` data, and the secondary `sat` track data. Each intersection
is identified by the `id`. New intersections or duplicates are identified according
to the parameters `dx`, `dt`, `tmf`, `tms`, `altmin`, `Xradius`, `expdist`,
`lidarprofile`, and `lidarrange`. Individual primary tracks are identified by their
`trackID`, which is used for error and warning messages. Experimental data is saved
for the `primspan`/`secspan` closest points to the intersection for primary/secondary
track points unless the closest measured point to the calculated intersection exceeds
the distance in meters of `expdist`. Only one type of satellite data is saved near
intersections (either profile or layer data) unless `savesecondsattype` is set to
`true`. Experimental data is read from the data files using the MATLAB sessions `ms`.
For new intersections, the counter is increased by 1, for duplicates the most accurate
calculation is saved.
"""
function add_intersections!(
  ms::mat.MSession,
  Xdata::DataFrame,
  tracked::DataFrame,
  accuracy::DataFrame,
  track::FlightTrack,
  sat::SatData,
  Xf::Tuple{T,T},
  Xs::Tuple{T,T},
  counter::Int,
  id::String,
  dx::T,
  dt::Dates.CompoundPeriod,
  tmf::DateTime,
  tms::DateTime,
  primspan::Int,
  secspan::Int,
  altmin::Real,
  trackID::Union{Missing,Int,AbstractString},
  Xradius::Real,
  expdist::Real,
  lidarprofile::NamedTuple,
  lidarrange::Tuple{Real,Real},
  saveobs::Union{String,Bool},
  savesecondsattype::Bool
) where T<:AbstractFloat
  # only CloudTrack: NA = T(NaN) # set precision of NaNs according to Float
  # Extract the DataFrame rows of the sat/flight data near the intersection
  Xflight, ift = get_flightdata(track, Xf, primspan, saveobs)
  alt = ift == 0 ? NaN : Xflight.data.alt[ift]
  cpro, clay, atmos, ist = get_satdata(ms, sat, Xs, secspan, alt, altmin,
    trackID, lidarprofile, lidarrange, saveobs, savesecondsattype, T)
  Xsat = sat.metadata.type == :CPro ? cpro.data : clay.data
  # Calculate accuracies (unless data is not saved, i.e. savedir = false)
  fxmeas, ftmeas, sxmeas, stmeas = ift == 0 ?
    (0, Dates.CompoundPeriod(), 0, Dates.CompoundPeriod()) :
    (dist.haversine(Xf,(Xflight.data.lat[ift], Xflight.data.lon[ift]),
    earthradius(Xf[1])),
    Dates.canonicalize(Dates.CompoundPeriod(tmf - Xflight.data.time[ift])),
    dist.haversine(Xs, (Xsat.lat[ist], Xsat.lon[ist]), earthradius(Xs[1])),
    Dates.canonicalize(Dates.CompoundPeriod(tms - Xsat.time[ist])))
  # Exclude data with long distances to nearest flight measurement
  if fxmeas > expdist || sxmeas > expdist
    @info("maximum distance of intersection to next track point exceeded; data excluded",
      trackID)
    return counter
  end
  # Save intersection data
  addX!(Xdata, tracked, accuracy, counter, Xf, id, dx, dt, Xradius, Xflight,
    cpro, clay, tmf, tms, atmos, fxmeas, ftmeas, sxmeas, stmeas, alt)
end #function add_intersections


"""
    function add_intersections!(
      ms::mat.MSession,
      Xdata::DataFrame,
      tracked::DataFrame,
      accuracy::DataFrame,
      sat::SatData,
      Xf::Tuple{T,T},
      Xs::Tuple{T,T},
      counter::Int,
      id::String,
      dx::T,
      dt::Dates.CompoundPeriod,
      tmf::DateTime,
      tms::DateTime,
      secspan::Int,
      altmin::Real,
      trackID::Union{Missing,Int,AbstractString},
      Xradius::Real,
      expdist::Real,
      lidarprofile::NamedTuple,
      lidarrange::Tuple{Real,Real},
      savesecondsattype::Bool
    ) where T<:AbstractFloat

Add intersection data to the DataFrames `Xdata`, `tracked`, and `accuracy` using
the coordinates calculated from the primary (`Xf`) and secondary (`Xs`) track, the
primary cloud `track` data, and the secondary `sat` track data. Each intersection
is identified by the `id`. New intersections or duplicates are identified according
to the parameters `dx`, `dt`, `tmf`, `tms`, `altmin`, `Xradius`, `expdist`,
`lidarprofile`, and `lidarrange`. Individual primary tracks are identified by their
`trackID`, which is used for error and warning messages. Experimental data is saved
for the `secspan` closest points to the intersection for secondary trajectories
unless the closest measured point to the calculated intersection exceeds
the distance in meters of `expdist`. Only one type of satellite data is saved near
intersections (either profile or layer data) unless `savesecondsattype` is set to
`true`. Experimental data is read from the data files using the MATLAB sessions `ms`.
For new intersections, the counter is increased by 1, for duplicates the most accurate
calculation is saved.
"""
function add_intersections!(
  ms::mat.MSession,
  Xdata::DataFrame,
  tracked::DataFrame,
  accuracy::DataFrame,
  sat::SatData,
  Xf::Tuple{T,T},
  Xs::Tuple{T,T},
  counter::Int,
  id::String,
  dx::T,
  dt::Dates.CompoundPeriod,
  tmf::DateTime,
  tms::DateTime,
  secspan::Int,
  altmin::Real,
  trackID::Union{Missing,Int,AbstractString},
  Xradius::Real,
  expdist::Real,
  lidarprofile::NamedTuple,
  lidarrange::Tuple{Real,Real},
  saveobs::Union{String,Bool},
  savesecondsattype::Bool
) where T<:AbstractFloat
  NA = T(NaN) # set precision of NaNs according to Float
  # Don't save additional cloud data near intersections at the moment
  Xcloud, ift = FlightTrack{T}(), 0
  cpro, clay, atmos, ist = get_satdata(ms, sat, Xs, secspan, NA, altmin,
    trackID, lidarprofile, lidarrange, saveobs, savesecondsattype, T)
  Xsat = sat.metadata.type == :CPro ? cpro.data : clay.data
  # Calculate accuracies
  fxmeas = NA
  ftmeas = Dates.canonicalize(Dates.CompoundPeriod())
  sxmeas = dist.haversine(Xs, (Xsat.lat[ist], Xsat.lon[ist]), earthradius(Xs[1]))
  stmeas = Dates.canonicalize(Dates.CompoundPeriod(tms - Xsat.time[ist]))
  # Exclude data with long distances to nearest flight measurement
  if fxmeas > expdist || sxmeas > expdist
    @info("maximum distance of intersection to next track point exceeded; data excluded",
      trackID)
    return counter
  end
  # Save intersection data
  addX!(Xdata, tracked, accuracy, counter, Xf, id, dx, dt, Xradius, Xcloud,
    cpro, clay, tmf, tms, atmos, fxmeas, ftmeas, sxmeas, stmeas, NA)
end #function add_intersections!


## Data extractions from raw data

"""
    find_timespan(sat::DataFrame, X::Tuple{<:AbstractFloat, <:AbstractFloat}, dataspan::Int=15)
      -> DataFrame, Vector{Int}

From the `sat` data in a `DataFrame` and the intersection `X` (as lat/lon pair),
find the time indices `t` ± `dataspan` for which the distance at `t` is minimal to `X`.

Return a `DataFrame` with `sat` data in the time span together with a `Vector{Int}`
holding the file indices of the corresponding granule file(s).
The `sat` data may be smaller than the `dataspan` at the edges of the `sat` `DataFrame`.
"""
function find_timespan(sat::DataFrame, X::Tuple{<:AbstractFloat, <:AbstractFloat},
  dataspan::Int=15)
  # Find index in sat data array with minimum distance to analytic intersection solution
  coords = ((sat.lat[i], sat.lon[i]) for i = 1:size(sat,1))
  imin = argmin(dist.haversine.(coords, [X], earthradius(X[1])))
  # Find first/last index of span acknowledging bounds of the data array
  t1 = max(1, min(imin-dataspan, length(sat.time)))
  t2 = min(length(sat.time), imin+dataspan)

  return (min=sat.time[t1], max=sat.time[t2]), unique(sat.fileindex[t1:t2])
end #function find_timespan


"""
    function get_flightdata(
      flight::FlightTrack{T},
      X::Tuple{<:AbstractFloat, <:AbstractFloat},
      primspan::Int
    ) where T -> track::FlightTrack, index::Int

From the measured `flight` data and lat/lon coordinates the intersection `X`,
save the closest measured value to the interpolated intersection ±`primspan` data points
to `track` and return it together with the `index` in `track` of the time with the coordinates
closest to `X`.
"""
function get_flightdata(
  flight::FlightTrack{T},
  X::Tuple{<:AbstractFloat, <:AbstractFloat},
  primspan::Int,
  saveobs::Union{String,Bool}
) where T
  # Retrieve observations, if saveobs is true
  (saveobs === false || isempty(saveobs)) && return FlightData{T}(), 0
  # Generate coordinate pairs from lat/lon columns
  coords = ((flight.data.lat[i], flight.data.lon[i]) for i = 1:size(flight.data,1))
  # Find the index (DataFrame row) of the intersection in the flight data
  imin = argmin(dist.haversine.(coords, [X], earthradius(X[1])))
  # Construct FlightTrack at Intersection
  t1 = max(1, min(imin-primspan, length(flight.data.time)))
  t2 = min(length(flight.data.time), imin+primspan)
  # t2 = t-timespan > length(data[:,1]) ? 0 : t2
  track = FlightData{T}(flight.data[t1:t2,:], flight.metadata)

  flightcoords = ((track.data.lat[i], track.data.lon[i])
    for i = 1:size(track.data, 1))
  return track, argmin(dist.haversine.(flightcoords, [X], earthradius(X[1])))
end #function get_flightdata


"""
    get_DateTimeRoute(filename::String, tzone::String)

From the `filename` and a custom time zone string (`tzone`), extract and return
the starting date, the standardised time zone, the flight ID, origin, and destination.
"""
function get_DateTimeRoute(filename::String, tzone::String)

    # Time is the first column and has to be addressed as flight[!,1] in the code
    # due to different column names, in which the timezone is included
    timezone = zonedict[tzone]
    # Retrieve date and metadata from filename
    flightID, datestr, course = try match(r"(.*?)_(.*?)_(.*)", filename).captures
    catch
      println()
      println()
      @warn "Flight ID, date, and course not found. Data skipped." file
      return missing, missing, missing, missing, missing
    end
    orig, dest = match(r"(.*)[-|_](.*)", course).captures
    date = try Dates.Date(datestr, "d-u-y", locale="english")
    catch
      println()
      println()
      @warn "Unable to parse date. Data skipped." file
      return missing, missing, missing, missing, missing
    end

    return date, timezone, flightID, orig, dest
end #function get_DateTimeRoute


"""
    get_satdata(
      ms::mat.MSession,
      sat::SatData,
      X::Tuple{<:AbstractFloat, <:AbstractFloat},
      secspan::Int,
      flightalt::Real,
      altmin::Real,
      flightid::Union{Int,String},
      lidarprofile::NamedTuple,
      lidarrange::Tuple{Real,Real},
      savesecondtype::Bool,
      Float::DataType=Float32
    ) -> cpro::CPro, clay::CLay, feature::Symbol, ts::Int

Using the `sat` data measurements within the overlap region and the MATLAB session
`ms`, extract CALIOP cloud profile (`cpro`) and/or layer data (`clay`) together with
the atmospheric conditions at flight level (`flightalt`) for the data point closest
to the calculated intersection `X` ± `secspan` timesteps and above the altitude
threshold `altmin`. In addition, return the index `ts` within `cpro`/`clay` of the
data point closest to `X`.
When `savesecondtype` is set to `false`, only the data type (`CLay`/`CPro`) in `sat`
is saved; if set to `true`, the corresponding data type is saved if available.
The lidar column data is saved for the height levels givin in the `lidarprofile` data
for the `lidarrange`. Floating point numbers are saved with single precision or
as defined by `Float`.
"""
function get_satdata(
  ms::mat.MSession,
  sat::SatData,
  X::Tuple{<:AbstractFloat, <:AbstractFloat},
  secspan::Int,
  flightalt::Real,
  altmin::Real,
  flightid::Union{Int,String},
  lidarprofile::NamedTuple,
  lidarrange::Tuple{Real,Real},
  saveobs::Union{String,Bool},
  savesecondtype::Bool,
  Float::DataType=Float32
)
  # Retrieve DataFrame at Intersection ± 15 time steps
  timespan, fileindex = find_timespan(sat.data, X, secspan)
  primfiles = map(f -> get(sat.metadata.files, f, 0), fileindex)
  secfiles = if sat.metadata.type == :CPro && savesecondtype
    replace.(primfiles, "CPro" => "CLay")
  elseif sat.metadata.type == :CLay && savesecondtype
    replace.(primfiles, "CLay" => "CPro")
  else
    String[]
  end

  # Get CPro/CLay data from near the intersection
  clay = if sat.metadata.type == :CLay
    CLay{Float}(ms, primfiles, timespan, lidarrange, altmin, saveobs)
  else
    try CLay{Float}(ms, secfiles, timespan, lidarrange, altmin)
    catch
      println(); @warn "could not load additional layer data" flightid
      CLay{Float}()
    end
  end
  cpro = if sat.metadata.type == :CPro
    CPro{Float}(ms, primfiles, timespan, lidarprofile, saveobs)
  else
    try CPro{Float}(ms, secfiles, timespan, lidarprofile, saveobs)
    catch
      println(); @warn "could not load additional profile data" flightid
      CPro{Float}()
    end
  end

  # Define primary data and index of intersection in primary data
  primdata = sat.metadata.type == :CPro ? cpro : clay
  coords = ((primdata.data.lat[i], primdata.data.lon[i]) for i = 1:size(primdata.data, 1))
  ts = argmin(dist.haversine.(coords, [X], earthradius(X[1])))

  # Get feature classification (atmospheric state)
  atmos = sat.metadata.type == :CPro ?
    atmosphericinfo(primdata, lidarprofile.fine, ts, flightalt, flightid) :
    atmosphericinfo(primdata, flightalt, ts)
  return cpro, clay, atmos, ts
end #function get_satdata


"""
    interpolate_time(data::DataFrame, X::Tuple{T,T}  where T<:AbstractFloat) -> DateTime

Return the linearly interpolated time at `X` (a lat/lon coordinate pair)
to the `data` in a DataFrame with a `time`, `lat`, and `lon` column.

Time is linearly interpolated between the 2 closest points to `X`.
"""
function interpolate_time(data::DataFrame, X::Tuple{T,T}  where T<:AbstractFloat)
  # Calculate distances for each coordinate pair to X
  d = dist.haversine.(((φ, λ) for (φ, λ) in zip(data.lat, data.lon)), [X], earthradius(X[1]))
  index = closest_points(d)
  d = dist.haversine((data.lat[index[1]], data.lon[index[1]]),
    (data.lat[index[2]], data.lon[index[2]]), earthradius(data.lat[index[1]]))
  ds = dist.haversine((data.lat[index[1]], data.lon[index[1]]), X, earthradius(data.lat[index[1]]))
  dt = data.time[index[2]] - data.time[index[1]]
  round(data.time[index[1]] + Dates.Millisecond(round(ds/d*dt.value)), Dates.Second)
end #function interpolate_time


## Helper functions for struct Instantiation

"""
    init_dict(data::Vector{Pair}, default)

From a vector of key, value pairs in `data`, return a DefaultDict with the given
`default` value. Ensure values are vectors by converting other data to vectors of
length 1.
"""
function init_dict(data::Vector{<:Pair}, default)
  dict = ds.DefaultDict(default)
  for (key, val) in data
    val isa Vector || (val = [val])
    dict[key] = val
  end

  return dict
end #function init_dict


"""
    trim_vec!(vec::Vector, cutoff::Int)

Trim `vec` to the first `cutoff` elements or return unchanged if `cutoff > length(vec)`.
"""
trim_vec!(vec::Vector, cutoff::Int) = deleteat!(vec, cutoff+1:length(vec))
