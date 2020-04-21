"""
    find_intersections(flight::FlightData, flighttracks::Vector, sat::SatDB,
      sattype::Symbol, sattracks::Vector, maxtimediff::Int, stepwidth::Float64, Xradius::Real,
      flightspan::Int, satspan::Int) -> idata::DataFrame, track::DataFrame, accuracy::DataFrame

Using interpolated `flighttracks` and `sattracks`, add new spatial and temporal
coordinates of the current flight along with the measured `flight` and `sat` `track`
near the intersection (±`flightspan`/±`satspan` datapoints of the intersection).

For the calculation of the satellite data, prefer data of `sattype` (`CLay` by default,
otherwise `CPro`) and find intersections using the `stepwidth` for the step width
of the interpolated data to find intersections with a maximum time delay `maxtimediff`
between the aircraft and satellite overpass. Intersections within `Xradius` in meters
are viewed as duplicates and only the one with the highest stepwidth is stored.
"""
function find_intersections(flight::FlightData, flighttracks::Vector, sat::SatDB,
  sattype::Symbol, sattracks::Vector, maxtimediff::Int, stepwidth::Float64, Xradius::Real,
  flightspan::Int, satspan::Int)

  # Initialise DataFrames for current flight
  idata = DataFrame(id=String[], lat=Float64[], lon=Float64[],
    tdiff=Dates.CompoundPeriod[], tflight = DateTime[], tsat = DateTime[], feature = Symbol[])
  track = DataFrame(id=String[], flight=FlightData[], sat=SatDB[])
  accuracy = DataFrame(id=String[], intersection=Float64[], flightcoord=Float64[],
    satcoord=Float64[], flighttime=Dates.CompoundPeriod[], sattime=Dates.CompoundPeriod[])
  counter = 0 # for intersections within the same flight used in the id
  # Get satellite data of the preferred type
  satdata = getfield(sat, sattype).data
  # Loop over sat and flight tracks
  for st in sattracks, ft in flighttracks
    if ft.min < st.max && ft.max > st.min
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
      for (lat, lon) in zip(ilat, st.track(ilat))
        push!(scoord, geo.LatLon(lat, lon))
      end

      # Calculate distances between each coordinate pair
      d = geo.distance.(fcoord, scoord)
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
      # Calculate minimum accuracy at equator
      dprec = geo.distance(geo.LatLon(0,0), geo.LatLon(0,0+stepwidth))
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
      any(isnan.(X)) && continue
      # Determine intersection from flight/satellite point of view and deviation
      Xf = flight.metadata.useLON ? geo.LatLon(ft.track(X[2]), X[2]) : geo.LatLon(X[1], ft.track(X[1]))
      Xs = geo.LatLon(X[1], st.track(X[1]))
      dx = geo.distance(Xf, Xs)
      dx < dprec || continue # ignore solutions outside tolerance
      # Determine time difference between aircraf/satellite at intersection
      tmf = flight.metadata.useLON ? timesec(ft.time(X[2])) : timesec(ft.time(X[1]))
      tms = timesec(st.time(X[1]))
      dt = Dates.canonicalize(Dates.CompoundPeriod(tms-tmf))
      # Consider only intersections within allowed time span
      Dates.Minute(-maxtimediff) < tmf - tms < Dates.Minute(maxtimediff) || continue
      # Look at previous intersection coordinates within the current flight
      dup = findfirst([geo.distance(Xf, geo.LatLon(idata[i,[:lat,:lon]]...))
        for i = 1:length(idata.id)] .< dprec)
      # Only save the most accurate intersection calculation within an Xradius
      # of the current intersection, i.e. only continue, if current intersection
      # is more accurate or new (not within Xradius)
      if isnothing(dup) || dx < accuracy.intersection[dup]
        # Extract the DataFrame rows of the sat/flight data near the intersection
        flightdata, satdb = get_trackdata(flight, sat, sattype, tmf, tms,
          flightspan, satspan)
        # Calculate accuracies
        i = length(flightdata.data.time)÷2+1
        fxmeas = geo.LatLon(flightdata.data.lat[i], flightdata.data.lon[i])
        ftmeas = Dates.canonicalize(Dates.CompoundPeriod(tmf - flightdata.data.time[i]))
        satdbdata = getfield(satdb,sattype).data
        i = length(satdbdata.time)÷2+1
        sxmeas = geo.LatLon(satdbdata.lat[i], satdbdata.lon[i])
        stmeas = Dates.canonicalize(Dates.CompoundPeriod(tms - satdbdata.time[i]))
        if isnothing(dup) # new data
          # Construct ID of current Intersection
          counter += 1
          id = string(flight.metadata.source,-,flight.metadata.dbID,-,counter)
          # Save intersection data
          push!(idata, (id=id, lat=Xf.lat, lon=Xf.lon, tdiff=dt,
            tflight = tmf, tsat = tms, feature=:no_signal))
          push!(track, (id=id, flight=flightdata, sat=satdb))
          # Save accuracies
          push!(accuracy, (id=id, intersection=dx, flightcoord=geo.distance(Xf,fxmeas),
            satcoord=geo.distance(Xs, sxmeas), flighttime=ftmeas, sattime=stmeas))
        else # more exact intersection calculations
          # Save intersection data
          idata[dup,:] = (id=id, lat=Xf.lat, lon=Xf.lon, tdiff=dt,
            tflight = tmf, tsat = tms, feature=:no_signal)
          track[dup,:] = (id=id, flight=flightdata, sat=satdb)
          # Save accuracies
          accuracy[dup,:] = (id=id, intersection=dx, flightcoord=geo.distance(Xf,fxmeas),
            satcoord=geo.distance(Xs, sxmeas), flighttime=ftmeas, sattime=stmeas)
        end
      end
    end #condition for track overlap
  end #loop over flight and sat tracks

  # Return intersection data of current flight
  return idata, track, accuracy
end #function find_intersections


"""
    findoverlap(flight::FlightData, sat::SatDB, sattype::Symbol, maxtimediff::Int)

From the data of the current `flight` and the `sat` data, calculate the data ranges
in the sat data that are in the vicinity of the flight track (min/max of lat/lon).
Prefer `sat` data of `sattype` (`:CLay`/`:CPro`) for the calculations and consider
only satellite data of ± `maxtimediff` minutes before the start and after the end of the
flight.
"""
function findoverlap(flight::FlightData, sat::SatDB, sattype::Symbol, maxtimediff::Int)

  # Initialise
  satranges = UnitRange[]; t1 = t2 = nothing; satdata = DataFrame();
  # Try to retrieve sat data in the range ±maxtimediff minutes before and after the flight
  # Try prefer sattype, and on failure other sattype
  for i = 1:2
    # Get sat data of sattype
    satdata = getfield(sat, sattype).data
    # Determine time span of flight ± maxtimediff minutes
    t1 = findfirst(satdata.time .≥ flight.data.time[1] - Dates.Minute(maxtimediff))
    t2 = findlast(satdata.time .≤ flight.data.time[end] + Dates.Minute(maxtimediff))
    t1 ≠ nothing && t2 ≠ nothing && break #exit on success, otherwise switch sattype
    sattype = swap_sattype(sattype)
  end
  # return empty ranges, if no complete overlap is found
  if t1 == nothing || t2 == nothing
    @warn string("no sufficient satellite data for time index ",
      "$(flight.data.time[1] - Dates.Minute(maxtimediff))...",
      "$(flight.data.time[end] + Dates.Minute(maxtimediff))")
    return (ranges=satranges, type=sattype)
  elseif length(t1:t2) ≤ 1
    @warn string("no sufficient overlap between flight/satellite data for flight ",
      "$(flight.metadata.dbID) at ",
      "$(satdata.time[t1]) ... $(satdata.time[t2])")
    return (ranges=satranges, type=sattype)
  end
  # Find overlaps in flight and sat data
  satoverlap = (flight.metadata.area.latmin .≤ satdata.lat[t1:t2] .≤ flight.metadata.area.latmax) .&
    ((flight.metadata.area.elonmin .≤ satdata.lon[t1:t2] .≤ flight.metadata.area.elonmax) .|
    (flight.metadata.area.wlonmin .≤ satdata.lon[t1:t2] .≤ flight.metadata.area.wlonmax))
  # Convert boolean vector of satellite overlapping data into ranges
  r = false # flag, whether index is part of a current range
  ind = 0   # index in the data array, when looping over data points
  for i = 1:length(satoverlap)
    if satoverlap[i] && !r #First data point of a range found
      r = true # flag as part of a range
      ind = i  # save start index
    elseif r && !satoverlap[i] && length(t1+ind-1:t1+i-2) > 1 # first index of non-overlapping data found
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
    satsegments = findflex(satdata.lat[r])

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
    interpolate_flightdata(ms::mat.MSession, flight::FlightData, stepwidth::Real)

Using the `flight` data, interpolate the data with the pchip method in the MATLAB
session (`ms`) with a `precission` in degrees.
"""
function interpolate_flightdata(ms::mat.MSession, flight::FlightData, stepwidth::Real)

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
    xr = interpolatedtrack(x[f.range], stepwidth)
    # Consider only data segments with more than one data point
    length(xr) > 1 || continue
    mat.eval_string(ms, string("try\npt = pchip(x, t);\ncatch\n",
      "disp('Error in flight time interpolation for')\n",
      "disp(['database/ID: ', string(meta.source), '/', string(meta.dbID)])\nend"))
    pt = mat.get_mvariable(ms, :pt)
    # Save the interpolated data to a vector
    flight.metadata.useLON ? (push!(flightdata, (lat = Minterpolate(ms, pf)(xr),
      lon =  xr, t = Dates.unix2datetime.(Minterpolate(ms, pt)(xr)), track = Minterpolate(ms, pf),
      time = Minterpolate(ms, pt), min = minimum(flight.data.lat[f.range]), max = maximum(flight.data.lat[f.range])))) :
      (push!(flightdata, (lat = xr, lon =  Minterpolate(ms, pf)(xr),
      t = Dates.unix2datetime.(Minterpolate(ms, pt)(xr)), track = Minterpolate(ms, pf),
      time = Minterpolate(ms, pt), min = f.min, max = f.max)))
  end

  # Return the interplated data
  return flightdata
end #function interpolate_flightdata


"""
    interpolatedtrack(xdata::Vector{Float64}, stepwidth::Real)

Return an interpolated vector of `xdata` with a step width of `stepwidth`.
Data are ascending or descending depending on the first and last data point of
`xdata`.
"""
interpolatedtrack(xdata::Vector{Float64}, stepwidth::Real) = xdata[1] < xdata[end] ?
  collect(Float64, xdata[1]:stepwidth:xdata[end]) : collect(Float64, xdata[1]:-stepwidth:xdata[end])
