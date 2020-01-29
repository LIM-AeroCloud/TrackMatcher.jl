function intersection(flights::FlightDB, sat::SatDB, sattype::Symbol=:CLay; deltat::Int=30,
  precision::Real=0.001)
  intersects = Intersection[]
  # New MATLAB session
  ms = mat.MSession()
  # Loop over data and interpolate track data and time, throw error on failure
  flight = nothing #Initialise variable flight
  try @pm.showprogress 1 "find intersections..." for outer flight in flights.inventory[1:10]
  # try @pm.showprogress 1 "find intersections..." for outer flight in
  #   [flights.inventory; flights.archive; flights.onlineData]
      satranges = get_satranges(flight, sat, sattype, deltat)
      sattracks = interpolate_satdata(ms, sat, satranges)
      flighttracks = interpolate_flightdata(ms, flight, precision)
      intersects = find_intersections(intersects, flight, flighttracks,
        sat, satranges.type, sattracks, deltat, precision)
    end #loop over flights
  catch
    printstyled(
      "\n\33[1mError:\33[0m Track data and/or time could not be interpolated in flight ",
      "$(flight.metadata.dbID) of $(flight.metadata.source) dataset\n",
      color=:red)
      println("Remaining data ignored.")
    return intersects
  finally #make sure MATLAB session is closed
    mat.close(ms)
  end
  return intersects
end #function intersection


"""
    find_intersections(flighttracks::Vector, sattracks::Vector, satranges::Vector{UnitRange})

documentation
"""
function find_intersections(intersects::Vector{Intersection}, flight::FlightData,
  flighttracks::Vector, sat::SatDB, sattype::Symbol, sattracks::Vector,
  deltat::Real, precision::Real)

  satdata = getfield(sat, sattype).data
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
      dprec = Geodesy.distance(Geodesy.LatLon(0,0), Geodesy.LatLon(0,0+prec))
      tm = Dates.unix2datetime(st.time(ilat[m]))
      if d[m] < dprec && Dates.Minute(-deltat) < ft.t[m] - tm < Dates.Minute(deltat)
        push!(intersects, Intersection(flight, sat, sattype, ft.t[m], tm, ilat[m], flon[m], d[m]))
      end
    end
  end
  return intersects
end #function find_intersections


function get_satranges(flight::FlightData, sat::SatDB, sattype::Symbol, deltat::Int)

  # Initialise
  satranges = UnitRange[]; t1 = t2 = nothing; satdata = DataFrame();
  for i = 1:2
    satdata = getfield(sat, sattype).data
    t1 = findfirst(satdata.time .≥ flight.data.time[1] - Dates.Minute(deltat))
    t2 = findlast(satdata.time .≤ flight.data.time[end] + Dates.Minute(deltat))
    !isnothing(t1) && !isnothing(t2) && break
    sattype = swap_sattype(sattype)
  end
  (isnothing(t1) || isnothing(t2)) && return (ranges=satranges, type=sattype)
  satoverlap = (flight.metadata.area.latmin .≤ satdata.lat[t1:t2] .≤ flight.metadata.area.latmax) .&
    ((flight.metadata.area.plonmin .≤ satdata.lon[t1:t2] .≤ flight.metadata.area.plonmax) .|
    (flight.metadata.area.nlonmin .≤ satdata.lon[t1:t2] .≤ flight.metadata.area.nlonmax))
  length(satoverlap) > 1 || return (ranges=satranges, type=sattype)
  r = false; ind = 0
  for i = 1:length(satoverlap)
    if satoverlap[i] && !r
      r = true
      ind = i
    elseif r && !satoverlap[i]
      r = false
      # push range where satellite is in correct time and spatial range
      # t1 - 1 corrects indices from the satoverlap array to the sat.CLay/sat.CPro array
      push!(satranges, t1+ind-1:t1+i-2)
    end
  end
  return (ranges=satranges, type=sattype)
end#function get_satranges


function interpolate_satdata(ms::mat.MSession,
  DB::SatDB, sat)

  satdata = getfield(DB, sat.type).data

  # Interpolate satellite tracks and flight times for all segments of interest
  interpolated_satdata = []
  # Loop over satellite data
  for r in sat.ranges
    # Find possible flex points in satellite tracks
    satsegments = findFlex(satdata.lat[r])

    # Loop over satellite segments
    for seg in satsegments
      length(seg) > 1 || continue #ignore points or empty data in segments
      # Pass variables to MATLAB
      mat.put_variable(ms, :x, satdata.lat[r])
      mat.put_variable(ms, :y, satdata.lon[r])
      # Convert times to UNIX times for interpolation
      mat.put_variable(ms, :t, Dates.datetime2unix.(satdata.time[r]))
      mat.eval_string(ms, "ps = pchip(x, y);")
      ps = mat.get_mvariable(ms, :ps)
      mat.eval_string(ms, "pt = pchip(x, t);")
      pt = mat.get_mvariable(ms, :pt)
      # Re-convert UNIX time values to DateTimes
      push!(interpolated_satdata, (track=Minterpolate(ms, ps), time=Minterpolate(ms, pt),
        min=seg.min, max=seg.max))
    end #loop over sat segments

  end

  return interpolated_satdata
end #function interpolate_satdata


"""


"""
function interpolate_flightdata(ms::mat.MSession, flight::FlightData, precision::Real)

  # Define x and y data based on useLON
  x, y = flight.metadata.useLON ?
    (flight.data.lon, flight.data.lat) : (flight.data.lat, flight.data.lon)
  # Interpolate flight tracks and tims for all segments
  flightdata = []
  for f in flight.metadata.flex
    mat.put_variable(ms, :x, x[f.range])
    mat.put_variable(ms, :y, y[f.range])
    mat.put_variable(ms, :t, Dates.datetime2unix.(flight.data.time[f.range]))
    mat.eval_string(ms, "pf = pchip(x, y);")
    pf = mat.get_mvariable(ms, :pf)
    xr = interpolatedtrack(x[f.range], precision)
    length(xr) > 1 || continue
    mat.eval_string(ms, "pt = pchip(x, t);")
    pt = mat.get_mvariable(ms, :pt)
    flight.metadata.useLON ? (push!(flightdata, (lat = Minterpolate(ms, pf)(xr),
      lon =  xr, t = Dates.unix2datetime.(Minterpolate(ms, pt)(xr)),
      min = minimum(flight.data.lat[f.range]), max = maximum(flight.data.lat[f.range])))) :
      (push!(flightdata, (lat = xr, lon =  Minterpolate(ms, pf)(xr),
      t = Dates.unix2datetime.(Minterpolate(ms, pt)(xr)), min = f.min, max = f.max)))
  end

  return flightdata
end #function interpolate_flightdata

interpolatedtrack(xdata::Vector{Float64}, precision::Real) = xdata[1] < xdata[end] ?
  collect(Float64, xdata[1]:precision:xdata[end]) : collect(Float64, xdata[1]:-precision:xdata[end])
