### Routines related to track/time interpolation and the calculation of track intersection points

"""
    intersection(flights::FlightDB, sat::SatDB, sattype::Symbol=:CLay;
      deltat::Int=30, flightspan::Int=0, satspan::Int=15, precision::Float64=0.01)
      -> Vector{Intersection}

From `flights` and `sat` data derive a `Vector{Intersection}` holding information
about all the intersection between the trajectories in both databases with a maximum
time difference of `deltat` in minutes between the flight and satellite overpass
at the intersection.

Flight and satellite data are interpolated with MATLAB's pchip routine to determine
the intersection. The stepwidth of the interpolation is performed with a `precision`
in degress (by default `0.001`). For the interpolation satellite data of `sattype`
is preferred (by default `:CLay`), only if no data is found the other satellite data
(`CPro` in default option) is used.

In each `Intersection`, the nearest `FlightData` measurement and `SatDB` measurements
are saved. Additionally, a variable time span (`flightspan`/`satspan`) can be saved including
`±n` timesteps additional to the nearest timestep (`±0` only the nearest point is saved;
default `flightspan=0` and `satspan=15`).
"""
function intersection(flights::FlightDB, sat::SatDB, sattype::Symbol=:CLay;
  deltat::Int=30, flightspan::Int=0, satspan::Int=15, precision::Float64=0.01)
  # Initialise Vector with Intersection data
  intersects = Intersection[]
  # New MATLAB session
  ms = mat.MSession()
  # Loop over data and interpolate track data and time, throw error on failure
  flight = nothing #Initialise variable flight
  try @pm.showprogress 1 "find intersections..." for outer flight in
    [flights.inventory; flights.archive; flights.onlineData]
      # Find sat tracks in the vicinity of flight tracks, where intersections are possible
      overlap = findoverlap(flight, sat, sattype, deltat)
      isempty(overlap.ranges) && continue
      # Interpolate trajectories using MATLAB's pchip routine
      sattracks = interpolate_satdata(ms, sat, overlap, flight.metadata)
      flighttracks = interpolate_flightdata(ms, flight, precision)
      # Calculate intersections and store in Vector
      intersects = find_intersections(intersects, flight, flighttracks,
        sat, overlap.type, sattracks, deltat, flightspan, satspan, precision)
    end #loop over flights
  catch
    # Issue warning on failure of interpolating track or time data
    printstyled(
      "\n\33[1mError:\33[0m Track data and/or time could not be interpolated in flight ",
      "$(flight.metadata.dbID) of $(flight.metadata.source) dataset\n",
      color=:red)
      println("Remaining data ignored.")
      # Return already calculate Intersections on failure
    return intersects
  finally #make sure MATLAB session is closed
    # Always close MATLAB session at the end or on failure
    mat.close(ms)
  end
  # Return Intersections after completion
  return intersects
end #function intersection


"""
    find_intersections(intersects::Vector{Intersection}, flight::FlightData,
      flighttracks::Vector, sat::SatDB, sattype::Symbol, sattracks::Vector,
      deltat::Int, flightspan::Int, satspan::Int, precision::Float64) -> Vector{Intersection}

To the vector `intersects` with intersections of sat and flight tracks, add new
`Intersection`s for the current `flight` and `sat` data using `sat` data of `sattype`
(`:CLay` or `:CPro`). Use the interpolated `flighttracks` and `sattracks` to find
the intersections with a `precision` (in degrees), where the maximum time difference
between the flight and satellite overpass at the intersection is `deltat` minutes.
In addition to the nearest messurement at the intersection save ± `flightspan`/`satspan`
data points of flight and sat data, respectively.
"""
function find_intersections(intersects::Vector{Intersection}, flight::FlightData,
  flighttracks::Vector, sat::SatDB, sattype::Symbol, sattracks::Vector,
  deltat::Int, flightspan::Int, satspan::Int, precision::Float64)

  # Get satellite data of the preferred type
  satdata = getfield(sat, sattype).data
  # Loop over sat and flight tracks
  for st in sattracks, ft in flighttracks
    if ft.min < st.max && ft.max > st.min
      # Use overlap of sat and flight data only and retrieve lat/lon values
      ilat = ft.lat[st.min .< ft.lat .< st.max]
      isempty(ilat) && continue  # skip data with no overlap region
      flon = ft.lon[st.min .< ft.lat .< st.max]
      # Interpolate sat data with same step width as flight data
      slon = st.track(ilat)

      # Calculate coordinate pairs of sat/flight data
      fcoord = Geodesy.LatLon[]
      for (lat, lon) in zip(ilat, flon)
        push!(fcoord, Geodesy.LatLon(lat, lon))
      end
      scoord = Geodesy.LatLon[]
      for (lat, lon) in zip(ilat, st.track(ilat))
        push!(scoord, Geodesy.LatLon(lat, lon))
      end

      # Calculate distances between each coordinate pair
      d = Geodesy.distance.(fcoord, scoord)
      # Find minimum in distance and check whether it is within precision
      m = argmin(d)
      # Calculate minimum accuracy at equator
      dprec = Geodesy.distance(Geodesy.LatLon(0,0), Geodesy.LatLon(0,0+precision))
      tm = Dates.unix2datetime(st.time(ilat[m]))
      if d[m] < dprec && Dates.Minute(-deltat) < ft.t[m] - tm < Dates.Minute(deltat)
        push!(intersects, Intersection(flight, sat, sattype, ft.t[m], tm, ilat[m], flon[m], flightspan, satspan, d[m]))
      end
    end
  end
  return intersects
end #function find_intersections


"""
    findoverlap(flight::FlightData, sat::SatDB, sattype::Symbol, deltat::Int)

From the data of the current `flight` and the `sat` data, calculate the data ranges
in the sat data that are in the vicinity of the flight track (min/max of lat/lon).
Prefer `sat` data of `sattype` (`:CLay`/`:CPro`) for the calculations and consider
only satellite data of ± `deltat` minutes before the start and after the end of the
flight.
"""
function findoverlap(flight::FlightData, sat::SatDB, sattype::Symbol, deltat::Int)

  # Initialise
  satranges = UnitRange[]; t1 = t2 = nothing; satdata = DataFrame();
  # Try to retrieve sat data in the range ±deltat minutes before and after the flight
  # Try prefer sattype, and on failure other sattype
  for i = 1:2
    # Get sat data of sattype
    satdata = getfield(sat, sattype).data
    # Determine time span of flight ± deltat minutes
    t1 = findfirst(satdata.time .≥ flight.data.time[1] - Dates.Minute(deltat))
    t2 = findlast(satdata.time .≤ flight.data.time[end] + Dates.Minute(deltat))
    t1 ≠ nothing && t2 ≠ nothing && break #exit on success, otherwise switch sattype
    sattype = swap_sattype(sattype)
  end
  # return empty ranges, if no complete overlap is found
  if t1 == nothing || t2 == nothing
    @warn string("no sufficient satellite data for time index ",
      "$(flight.data.time[1] - Dates.Minute(deltat))...",
      "$(flight.data.time[end] + Dates.Minute(deltat))")
    return (ranges=satranges, type=sattype)
  end
  # Find overlaps in flight and sat data
  satoverlap = (flight.metadata.area.latmin .≤ satdata.lat[t1:t2] .≤ flight.metadata.area.latmax) .&
    ((flight.metadata.area.elonmin .≤ satdata.lon[t1:t2] .≤ flight.metadata.area.elonmax) .|
    (flight.metadata.area.wlonmin .≤ satdata.lon[t1:t2] .≤ flight.metadata.area.wlonmax))
  # Return only overlaps with at least 2 data points
  length(satoverlap) > 1 || return (ranges=satranges, type=sattype)
  # Convert boolean vector of satellite overlapping data into ranges
  r = false # flag, whether index is part of a current range
  ind = 0   # index in the data array, when looping over data points
  for i = 1:length(satoverlap)
    if satoverlap[i] && !r #First data point of a range found
      r = true # flag as part of a range
      ind = i  # save start index
    elseif r && !satoverlap[i] # first index of non-overlapping data found
      r = false # flag as non-overlapping data
      # push range where satellite is in correct time and spatial range
      # t1 - 1 corrects indices from the satoverlap array to the sat.CLay/sat.CPro array
      push!(satranges, t1+ind-1:t1+i-2) # save range from saved index to last index
    end
  end
  # Return tuple with sat ranges, and the type of sat data that was used for the calculations
  return (ranges=satranges, type=sattype)
end#function findoverlap


"""
    interpolate_satdata(ms::mat.MSession, DB::SatDB, sat)

Using the satellite data in the `DB` database, and the stored `sat` ranges and types,
interpolate the data with the pchip method in the MATLAB session (`ms`).
"""
function interpolate_satdata(ms::mat.MSession, DB::SatDB, overlap, flight::FlightMetadata)

  # Get the satellite data of the correct type
  satdata = getfield(DB, overlap.type).data

  # Interpolate satellite tracks and flight times for all segments of interest
  interpolated_satdata = []
  # Loop over satellite data
  for r in overlap.ranges
    # Find possible flex points in satellite tracks
    satsegments = findFlex(satdata.lat[r])

    # Loop over satellite segments
    for seg in satsegments
      length(seg.range) > 1 || continue #ignore points or empty data in segments
      # Pass variables to MATLAB
      mat.put_variable(ms, :x, satdata.lat[r][seg.range])
      mat.put_variable(ms, :y, satdata.lon[r][seg.range])
      mat.put_variable(ms, :meta, flight)
      # Convert times to UNIX times for interpolation
      mat.put_variable(ms, :t, Dates.datetime2unix.(satdata.time[r][seg.range]))
      mat.eval_string(ms, string("try\nps = pchip(x, y);\ncatch\n",
        "disp('Error in sat track interpolation for')\n",
        "disp(['database/ID: ', string(meta.source), '/', string(meta.dbID)])\nend"))
      ps = mat.get_mvariable(ms, :ps)
      mat.eval_string(ms, string("try\npt = pchip(x, t);\ncatch\n",
        "disp('Error in sat time interpolation for')\n",
        "disp(['database/ID: ', string(meta.source), '/', string(meta.dbID)])\nend"))
      pt = mat.get_mvariable(ms, :pt)
      # Re-convert UNIX time values to DateTimes
      # Save the interpolation functions for track and time data for the current dataset
      # together with the min/max of the latitude (allowed x value range)
      push!(interpolated_satdata, (track=Minterpolate(ms, ps), time=Minterpolate(ms, pt),
        min=seg.min, max=seg.max))
    end #loop over sat segments
  end #loop over sat ranges

  # Return a vector with interpolation functions for each dataset
  return interpolated_satdata
end #function interpolate_satdata


"""
    interpolate_flightdata(ms::mat.MSession, flight::FlightData, precision::Real)

Using the `flight` data, interpolate the data with the pchip method in the MATLAB
session (`ms`) with a `precission` in degrees.
"""
function interpolate_flightdata(ms::mat.MSession, flight::FlightData, precision::Real)

  # Define x and y data based on useLON
  x, y = flight.metadata.useLON ?
    (flight.data.lon, flight.data.lat) : (flight.data.lat, flight.data.lon)
  # Interpolate flight tracks and tims for all segments
  flightdata = []
  for f in flight.metadata.flex
    # Hand over variables to MATLAB
    mat.put_variable(ms, :x, x[f.range])
    mat.put_variable(ms, :y, y[f.range])
    mat.put_variable(ms, :meta, flight.metadata)
    # Convert times to UNIX times for interpolation
    mat.put_variable(ms, :t, Dates.datetime2unix.(flight.data.time[f.range]))
    # Interpolate with MATLAB
    mat.eval_string(ms, string("try\npf = pchip(x, y);\ncatch\n",
      "disp('Error in flight track interpolation for')\n",
      "disp(['database/ID: ', string(meta.source), '/', string(meta.dbID)])\nend"))
    # Retrieve variables from MATLAB
    pf = mat.get_mvariable(ms, :pf)
    xr = interpolatedtrack(x[f.range], precision)
    # Consider only data segments with more than one data point
    length(xr) > 1 || continue
    mat.eval_string(ms, string("try\npt = pchip(x, t);\ncatch\n",
      "disp('Error in flight time interpolation for')\n",
      "disp(['database/ID: ', string(meta.source), '/', string(meta.dbID)])\nend"))
    pt = mat.get_mvariable(ms, :pt)
    # Save the interpolated data to a vector
    flight.metadata.useLON ? (push!(flightdata, (lat = Minterpolate(ms, pf)(xr),
      lon =  xr, t = Dates.unix2datetime.(Minterpolate(ms, pt)(xr)),
      min = minimum(flight.data.lat[f.range]), max = maximum(flight.data.lat[f.range])))) :
      (push!(flightdata, (lat = xr, lon =  Minterpolate(ms, pf)(xr),
      t = Dates.unix2datetime.(Minterpolate(ms, pt)(xr)), min = f.min, max = f.max)))
  end

  # Return the interplated data
  return flightdata
end #function interpolate_flightdata


"""
    interpolatedtrack(xdata::Vector{Float64}, precision::Real)

Return an interpolated vector of `xdata` with a step width of `precision`.
Data are ascending or descending depending on the first and last data point of
`xdata`.
"""
interpolatedtrack(xdata::Vector{Float64}, precision::Real) = xdata[1] < xdata[end] ?
  collect(Float64, xdata[1]:precision:xdata[end]) : collect(Float64, xdata[1]:-precision:xdata[end])
