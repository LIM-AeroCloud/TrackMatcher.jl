function intersection(flights::FlightDB, sat::SatDB, satdata::Symbol=:CLay; deltat::Int=30,
  precision::Real=0.001)
  intersects = Intersection[]
  # New MATLAB session
  ms = mat.MSession()
  # Loop over data and interpolate track data and time, throw error on failure
  # try
    @pm.showprogress 1 "find intersections..." for flight in flights.inventory[1:10]
    # try @pm.showprogress 1 "interpolate data..." for flight in [flights.inventory; flights.archive; flights.onlineData]
      satranges = get_satranges(flight, getfield(sat, satdata), deltat)
      sattracks = interpolate_satdata(ms, getfield(sat, satdata), satranges)
      flighttracks = interpolate_flightdata(ms, flight, precision)
      intersects = find_intersections(intersects, flight, flighttracks,
        getfield(sat, satdata), satdata, sattracks, deltat, precision)
    end #loop over flights
  # catch
  #   throw("Track data and/or time could not be interpolated")
  # finally #make sure MATLAB session is closed
    mat.close(ms)
  # end
  return intersects
end #function intersection


"""
    find_intersections(flighttracks::Vector, sattracks::Vector, satranges::Vector{UnitRange})

documentation
"""
function find_intersections(intersects::Vector{Intersection}, flight::FlightData,
  flighttracks::Vector, sat::Union{CLay,CPro}, sattype::Symbol, sattracks::Vector,
  deltat::Real, precision::Real)

  for st in sattracks, ft in flighttracks
    if ft.min < st.max && ft.max > st.min
      # Use overlap of sat and flight data only
      fi = ft.lat[st.min .< ft.lat .< st.max]
      # Interpolate sat data with same step width as flight data
      si = st.track(fi)

      # Find flight longitudes in overlapping sat/flight range
      fmin, fmax = fi[1] < fi[end] ? (fi[1], fi[end]) : (fi[end], fi[1])
      flon = ft.lon[fmin .<= ft.lat .<= fmax]
      # Calculate distances between each coordinate pair
      d = abs.(abs.(si) .- abs.(flon))
      # Find minimum in distance and check whether it is within precision
      m = argmin(d)
      tm = Dates.unix2datetime(st.time(fi[m]))
      if d[m] < precision && Dates.Minute(-deltat) < ft.t[m] - tm < Dates.Minute(deltat)
        push!(intersects, Intersection(flight, sat, sattype, ft.t[m], tm, fi[m], flon[m]))
      end
    end
  end
  return intersects
end #function find_intersections


function get_satranges(flight::FlightData, sat::Union{CLay,CPro}, deltat::Int)::Vector{UnitRange}
  satranges = UnitRange[]
  t1 = findfirst(sat.time .≥ flight.time[1] - Dates.Minute(deltat))
  t2 = findlast(sat.time .≤ flight.time[end] + Dates.Minute(deltat))
  (isnothing(t1) || isnothing(t2)) && return satranges
  satoverlap = (flight.metadata.area.latmin .≤ sat.lat[t1:t2] .≤ flight.metadata.area.latmax) .&
    ((flight.metadata.area.plonmin .≤ sat.lon[t1:t2] .≤ flight.metadata.area.plonmax) .|
    (flight.metadata.area.nlonmin .≤ sat.lon[t1:t2] .≤ flight.metadata.area.nlonmax))
  length(satoverlap) > 1 || return satranges
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
  return satranges
end#function get_satranges


function interpolate_satdata(ms::mat.MSession,
  sat::Union{CLay,CPro}, satranges::Vector{UnitRange})

  # Interpolate satellite tracks and flight times for all segments of interest
  interpolated_satdata = []
  # Loop over satellite data
  for r in satranges
    # Find possible flex points in satellite tracks
    satsegments = findFlex(sat.lat[r])

    # Loop over satellite segments
    for seg in satsegments
      length(seg) > 1 || continue #ignore points or empty data in segments
      # Pass variables to MATLAB
      mat.put_variable(ms, :x, sat.lat[r])
      mat.put_variable(ms, :y, sat.lon[r])
      # Convert times to UNIX times for interpolation
      mat.put_variable(ms, :t, Dates.datetime2unix.(sat.time[r]))
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
  x, y = flight.metadata.useLON ? (flight.lon, flight.lat) : (flight.lat, flight.lon)
  # Interpolate flight tracks and tims for all segments
  flightdata = []
  for f in flight.metadata.flex
    mat.put_variable(ms, :x, x[f.range])
    mat.put_variable(ms, :y, y[f.range])
    mat.put_variable(ms, :t, Dates.datetime2unix.(flight.time[f.range]))
    mat.eval_string(ms, "pf = pchip(x, y);")
    pf = mat.get_mvariable(ms, :pf)
    xr = interpolatedtrack(x[f.range], precision)
    length(xr) > 1 || continue
    mat.eval_string(ms, "pt = pchip(x, t);")
    pt = mat.get_mvariable(ms, :pt)
    flight.metadata.useLON ? (push!(flightdata, (lat = Minterpolate(ms, pf)(xr),
      lon =  xr, t = Dates.unix2datetime.(Minterpolate(ms, pt)(xr)),
      min = minimum(flight.lat[f.range]), max = maximum(flight.lat[f.range])))) :
      (push!(flightdata, (lat = xr, lon =  Minterpolate(ms, pf)(xr),
      t = Dates.unix2datetime.(Minterpolate(ms, pt)(xr)), min = f.min, max = f.max)))
  end

  return flightdata
end #function interpolate_flightdata

interpolatedtrack(xdata::Vector{Float64}, precision::Real) = xdata[1] < xdata[end] ?
  collect(Float64, xdata[1]:precision:xdata[end]) : collect(Float64, xdata[1]:-precision:xdata[end])
