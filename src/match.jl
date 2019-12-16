function intersection(flights::FlightDB, sat::SatDB; deltat::Int=30, satdata::Symbol=CPro)
  hits = FlightData[]
  for flight in [flights.inventory; flights.archive; flights.onlineData]
    satranges = get_satranges(flight, getfield(sat, satdata), deltat)
  end #loop over flights
end #function intersection

function get_satranges(flight::FlightData, sat::Union{CLay,CPro}, deltat::Int)::Vector{UnitRange}
  satranges = UnitRange[]
  t1 = findfirst(sat.time .≥ flight.time[1] - Dates.Minute(deltat))
  t2 = findlast(sat.time .≤ flight.time[end] + Dates.Minute(deltat))
  (isnothing(t1) || isnothing(t2)) && return satranges
  satoverlap = (flight.metadata.area.latmin .≤ sat.lat[t1:t2] .≤ flight.metadata.area.latmax) .&
    ((flight.metadata.area.plonmin .≤ sat.lon[t1:t2] .≤ flight.metadata.area.plonmax) .|
    (flight.metadata.area.nlonmin .≤ sat.lon[t1:t2] .≤ flight.metadata.area.nlonmax))
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
