function intersection(flights::FlightDB, sat::SatDB; deltat::Int=30,
  satdata::Symbol=:CPro, precision::Real=0.001)
  hits = []
  # New MATLAB session
  ms = mat.MSession()
  # Loop over data and interpolate track data and time, throw error on failure
  try
    @pm.showprogress 1 "interpolate data..." for flight in flights.inventory[1:10]
    # try @pm.showprogress 1 "interpolate data..." for flight in [flights.inventory; flights.archive; flights.onlineData]
      satranges = get_satranges(flight, getfield(sat, satdata), deltat)
      sattracks = interpolate_satdata(ms, getfield(sat, satdata), satranges)
      flighttracks = interpolate_flightdata(ms, flight, precision)
      push!(hits, (flight = flighttracks, sat = sattracks))
    end #loop over flights
  catch
    throw("Track data and/or time could not be interpolated")
  finally #make sure MATLAB session is closed
    mat.close(ms)
  end
  return hits
end #function intersection


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
  fx = Vector{Float64}[]; fy = Vector{Float64}[]; ft = Vector{DateTime}[]
  for f in flight.metadata.flex
    mat.put_variable(ms, :x, x[f.range])
    mat.put_variable(ms, :y, y[f.range])
    mat.put_variable(ms, :t, Dates.datetime2unix.(flight.time[f.range]))
    mat.eval_string(ms, "p = pchip(x, y);")
    p = mat.get_mvariable(ms, :p)
    xr = interpolatedtrack(x[f.range], precision)
    length(xr) > 1 || continue
    push!(fx, xr)
    push!(fy, Minterpolate(ms, p)(xr))
    mat.eval_string(ms, "p = pchip(x, t);")
    p = mat.get_mvariable(ms, :p)
    push!(ft, Dates.unix2datetime.(Minterpolate(ms, p)(xr)))
  end

  return flight.metadata.useLON ? (lat = fy, lon = fx, t = ft) : (lat = fx, lon = fy, t = ft)
end #function interpolate_flightdata

interpolatedtrack(xdata::Vector{Float64}, precision::Real) = xdata[1] < xdata[end] ?
  collect(Float64, xdata[1]:precision:xdata[end]) : collect(Float64, xdata[1]:-precision:xdata[end])
