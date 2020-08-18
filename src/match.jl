"""
    find_intersections(
      ms::mat.MSession,
      flight::FlightData,
      flighttracks::Vector,
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
      savesecondsattype::Bool
    ) -> Xdata::DataFrame, track::DataFrame, accuracy::DataFrame

Using interpolated `flighttracks` and `sattracks`,
add new spatial and temporal coordinates of all intersections of the current flight
with satellite tracks to `Xdata`, if the overpass of the aircraft and the satellite
at the intersection is within `maxtimediff` minutes. Additionally, save the measured
`flight` and `sat` `tracked` data near the intersection (±`flightspan`/±`satspan`
datapoints of the intersection) and information about the `accuracy`.

When `savesecondsattype` is set to true, the additional satellite data type not
used to derive the intersections from the `SatData` is stored as well in `Intersection`.
Satellite column data is stored over the `lidarrange` and defined by the `lidarprofile`.

The algorithm finds intersections, by finding the minimum distance between flight and
sat track points in the overlap region of both tracks. To be an intersection,
the minimum distance of the interpolated track points must be below a threshold `dmin`.
For duplicate intersection finds within an `Xradius`, only the one with the highest
accuracy (lowest `dmin`) is saved.
"""
function find_intersections(
  ms::mat.MSession,
  flight::FlightData,
  flighttracks::Vector,
  sat::SatData,
  sattracks::Vector,
  maxtimediff::Int,
  dmin::Real,
  Xradius::Real,
  lidarprofile::NamedTuple,
  lidarrange::Tuple{Real,Real},
  flightspan::Int,
  satspan::Int,
  savesecondsattype::Bool
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
    # Use overlap of sat and flight data only and retrieve lat/lon values
    ilat = ft.lat[st.min .< ft.lat .< st.max]
    length(ilat) > 1 || continue  # skip data with no overlap region
    flon = ft.lon[st.min .< ft.lat .< st.max]
    # Interpolate sat data with same step width as flight data
    slon = st.track(ilat)

    # Calculate coordinate pairs of sat/flight data
    fcoord = geo.LatLon[]
    for (lat, lon) in zip(ilat, flon)
      push!(fcoord, geo.LatLon(lat, lon))
    end
    scoord = geo.LatLon[]
    for (lat, lon) in zip(ilat, slon)
      push!(scoord, geo.LatLon(lat, lon))
    end

    # Calculate distances between each coordinate pair
    d = geo.distance.(fcoord, scoord)
    # Loop over track segment until all intersections within precision are found
    while length(fcoord) > 1
      # Find minimal distance and second least adjacent distance
      m1 = argmin(d)
      m2 = if m1 == 1
        2
      elseif m1 == length(d)
        m1 - 1
      elseif d[m1-1] < d[m1+1]
        m1-1
      else
        m1+1
      end
      # Determine intersection between the 2 point pairs analytically by linear interpolation
      X = (((fcoord[m1].lat*fcoord[m2].lon - fcoord[m1].lon*fcoord[m2].lat) *
        (scoord[m1].lat - scoord[m2].lat) - (fcoord[m1].lat - fcoord[m2].lat) *
        (scoord[m1].lat*scoord[m2].lon - scoord[m1].lon*scoord[m2].lat)) /
        ((fcoord[m1].lat - fcoord[m2].lat)*(scoord[m1].lon - scoord[m2].lon) -
        (fcoord[m1].lon - fcoord[m2].lon)*(scoord[m1].lat - scoord[m2].lat)),
        ((fcoord[m1].lat*fcoord[m2].lon - fcoord[m1].lon*fcoord[m2].lat) *
        (scoord[m1].lon - scoord[m2].lon) - (fcoord[m1].lon - fcoord[m2].lon) *
        (scoord[m1].lat*scoord[m2].lon - scoord[m1].lon*scoord[m2].lat)) /
        ((fcoord[m1].lat - fcoord[m2].lat)*(scoord[m1].lon - scoord[m2].lon) -
        (fcoord[m1].lon - fcoord[m2].lon)*(scoord[m1].lat - scoord[m2].lat)))
      # Filter mismatches and unreasonable analytic solutions
      any(isnan.(X)) && break
      -90 .≤ X[1] .≤ 90 || break
      -180 .≤ X[2] .≤ 180 || break
      # Determine intersection from flight/satellite point of view and deviation
      Xf = try flight.metadata.useLON ?
        geo.LatLon(ft.track(X[2]), X[2]) : geo.LatLon(X[1], ft.track(X[1]))
      catch
        break
      end
      Xs = try geo.LatLon(X[1], st.track(X[1]))
      catch
        break
      end
      dx = geo.distance(Xf, Xs)
      dx < dmin || break # ignore solutions outside tolerance
      # Determine time difference between aircraf/satellite at intersection
      tmf = flight.metadata.useLON ? timesec(ft.time(X[2])) : timesec(ft.time(X[1]))
      tms = timesec(st.time(X[1]))
      dt = Dates.canonicalize(Dates.CompoundPeriod(tms-tmf))
      # Consider only intersections within allowed time span
      abs(tmf - tms) < Dates.Minute(maxtimediff) || break
      # Look at previous intersections at the coordinates within Xradius
      dup = findfirst([geo.distance(Xf, geo.LatLon(Xdata[i,[:lat,:lon]]...))
        for i = 1:length(Xdata.id)] .< Xradius)
      # Only save the most accurate intersection calculation within an Xradius
      # of the current intersection, i.e. only continue, if current intersection
      # is more accurate or new (not within Xradius)
      if dup === nothing || dx < accuracy.intersection[dup]
        # Extract the DataFrame rows of the sat/flight data near the intersection
        Xflight, ift = get_flightdata(flight, X, flightspan)
        cpro, clay, feature, ist = get_satdata(ms, sat, X, satspan,
          Xflight.data.alt[ift], Xflight.metadata.dbID, lidarprofile, lidarrange, savesecondsattype)
        Xsat = sat.metadata.type == :CPro ? cpro.data : clay.data
        # Calculate accuracies
        fxmeas = geo.LatLon(Xflight.data.lat[ift], Xflight.data.lon[ift])
        ftmeas = Dates.canonicalize(Dates.CompoundPeriod(tmf - Xflight.data.time[ift]))
        sxmeas = geo.LatLon(Xsat.lat[ist], Xsat.lon[ist])
        stmeas = Dates.canonicalize(Dates.CompoundPeriod(tms - Xsat.time[ist]))

        if dup === nothing # new data
          # Construct ID of current Intersection
          counter += 1
          id = string(flight.metadata.source,-,flight.metadata.dbID,-,counter)
          # Save intersection data
          push!(Xdata, (id=id, lat=Xf.lat, lon=Xf.lon, alt=Xflight.data.alt[ift],
            tdiff=dt, tflight = tmf, tsat = tms, feature=feature))
          push!(track, (id=id, flight=Xflight, CPro=cpro, CLay = clay))
          # Save accuracies
          push!(accuracy, (id=id, intersection=dx, flightcoord=geo.distance(Xf,fxmeas),
            satcoord=geo.distance(Xs, sxmeas), flighttime=ftmeas, sattime=stmeas))
        else # more exact intersection calculations
          # Find most accurate solution
          dx < accuracy.intersection[dup] || break
          # Save intersection data
          Xdata[dup,:] = (id=Xdata.id[dup], lat=Xf.lat, lon=Xf.lon, alt=Xflight.data.alt[ift],
            tdiff=dt, tflight = tmf, tsat = tms, feature=feature)
          track[dup,:] = (id=Xdata.id[dup], flight=Xflight, CPro=cpro, CLay = clay)
          # Save accuracies
          accuracy[dup,:] = (id=Xdata.id[dup], intersection=dx, flightcoord=geo.distance(Xf,fxmeas),
            satcoord=geo.distance(Xs, sxmeas), flighttime=ftmeas, sattime=stmeas)
        end
      end
      # Reduce distance, flight, and sat coordinates by current indices, used for analytic solution of X
      deleteat!(d, sort([m1,m2]))
      deleteat!(scoord, sort([m1,m2]))
      deleteat!(fcoord, sort([m1,m2]))
      # Remove further entries within Xradius of the current intersection find
      d = d[geo.distance.(geo.LatLon(X...),fcoord) .> Xradius]
      scoord = scoord[geo.distance.(geo.LatLon(X...),fcoord) .> Xradius]
      fcoord = fcoord[geo.distance.(geo.LatLon(X...),fcoord) .> Xradius]
    end #loop finding intersections in track segment
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
  overlap = NamedTuple{(:range, :min, :max),Tuple{UnitRange, Real, Real}}[]; t1 = t2 = nothing
  ## Retrieve sat data in the range ±maxtimediff minutes before and after the flight
  # Set time span
  t1 = findfirst(sat.data.time .≥ flight.data.time[1] - Dates.Minute(maxtimediff))
  t2 = findlast(sat.data.time .≤ flight.data.time[end] + Dates.Minute(maxtimediff))
  # return empty ranges, if no complete overlap is found
  if t1 === nothing || t2 === nothing
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
      satsegments = findflex(sat.data.lat[seg])
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
    interpolate_satdata(sat::SatData, overlap::Vector{UnitRange})
      -> Vector{Any}

Using the `sat` data and the stored `overlap` ranges, interpolate the data
with the pchip method. Use the metadata in `flight` for error reports.
"""
function interpolate_satdata(
  sat::SatData,
  overlap::Vector{NamedTuple{(:range, :min, :max),Tuple{UnitRange, Real, Real}}}
)

  # Interpolate satellite tracks and flight times for all segments of interest
  idata = []
  # Loop over satellite data
  for r in overlap
    # Interpolate track data with PCHIP
    ps = pchip(sat.data.lat[r.range], sat.data.lon[r.range])
    # Convert times to UNIX times for interpolation
    unixtime = Dates.datetime2unix.(sat.data.time[r.range])
    pt = pchip(sat.data.lat[r.range], unixtime)
    # Re-convert UNIX time values to DateTimes
    # Save the interpolation functions for track and time data for the current dataset
    # together with the min/max of the latitude (allowed x value range)
    push!(idata, (track=interpolate(ps), time=interpolate(pt), min=r.min, max=r.max))
  end #loop over sat ranges

  # Return a vector with interpolation functions for each dataset
  return idata
end #function interpolate_satdata


"""
    interpolate_flightdata(flight::FlightData, stepwidth::Real)

Using the `flight` data, interpolate the data with the pchip method
using the defined `stepwidth` for interpolation.
"""
function interpolate_flightdata(flight::FlightData, stepwidth::Real)

  # Define x and y data based on useLON
  x, y = flight.metadata.useLON ?
    (flight.data.lon, flight.data.lat) : (flight.data.lat, flight.data.lon)
  # Interpolate flight tracks and tims for all segments
  flightdata = []
  for f in flight.metadata.flex
    # Interpolate track data with PCHIP
    pf = pchip(x[f.range],y[f.range])
    xi = interpolatedtrack(x[f.range], stepwidth)
    yi = interpolate.(pf, xi)
    # Consider only data segments with more than one data point
    length(xi) > 1 || continue
    # Convert times to UNIX times and interpolate with PCHIP
    unixtime = Dates.datetime2unix.(flight.data.time[f.range])
    pt = pchip(x[f.range], unixtime)
    ti = Dates.unix2datetime.(interpolate.(pt, xi))
    # Save the interpolated data to a vector
    if flight.metadata.useLON
      push!(flightdata, (lat = yi, lon =  xi, t = ti, track = interpolate(pf), time = interpolate(pt),
        min = minimum(flight.data.lat[f.range]), max = maximum(flight.data.lat[f.range])))
    else
      push!(flightdata, (lat = xi, lon =  yi, t = ti,
        track = interpolate(pf), time = interpolate(pt), min = f.min, max = f.max))
    end
  end

  # Return the interplated data
  return flightdata
end #function interpolate_flightdata


"""
    interpolatedtrack(xdata::Vector{<:AbstractFloat}, stepwidth::Real)

Return an interpolated vector of `xdata` with the defined `stepwidth`.
Data are ascending or descending depending on the first and last data point of
`xdata`.
"""
interpolatedtrack(xdata::Vector{<:AbstractFloat}, stepwidth::Real) = xdata[1] < xdata[end] ?
  collect(xdata[1]:stepwidth:xdata[end]) : collect(xdata[1]:-stepwidth:xdata[end])


  """ Generate a function to interpolate x vectors from the current PCHIP struct """
  interpolate(pc::PCHIPdata{<:AbstractFloat}) = x::Union{<:AbstractFloat,Vector{<:AbstractFloat}} -> interpolate(pc, x)
