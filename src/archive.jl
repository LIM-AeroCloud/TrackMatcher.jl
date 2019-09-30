"""
    loadArchive(files)

documentation
"""
function loadArchive(files)
  # Initialise inventory file array
  archive = FlightData[]
  @pm.showprogress 1 "load archive..." for file in files
    # Load data
    flights = CSV.read(file, datarow=2, header=[:id, :ident, :orig, :dest,
      :aircraft, :time, :lat, :lon, :speed, :alt, :climb, :heading, :direction,
      :facility, :description, :est])
    # Calculate time from individual columns and add as DateTime to DataFrame
    flights.time = ZonedDateTime.(DateTime.(flights.time, "m/d/y H:M:S"), tz.tz"UTC")

    # Initialise loop over file
    global iStart = 1
    # Loop over all data points
    for i = 1:length(flights.time)-1
      if flights.id[i] ≠ flights.id[i+1]
        # When flight ID changes save data as FlightData and set start index to next dataset
        push!(archive, FlightData(flights.time[iStart:i], flights.lat[iStart:i],
          flights.lon[iStart:i], flights.alt[iStart:i], flights.heading[iStart:i],
          flights.climb[iStart:i], flights.speed[iStart:i]))
        global iStart = i+1
      end
    end
    # Save last flight of the dataset
    i = length(flights.time)
    push!(archive, FlightData(flights.time[iStart:i], flights.lat[iStart:i],
      flights.lon[iStart:i], flights.alt[iStart:i], flights.heading[iStart:i],
      flights.climb[iStart:i], flights.speed[iStart:i]))
  end

  return archive
end #function loadArchive
