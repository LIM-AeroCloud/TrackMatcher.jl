"""
    loadArchive(files)

documentation
"""
function loadArchive(files)
  # Initialise inventory file array
  archive = flightData[]
  for file in files
    # Load data
    flights = DataFrame(csv.load(files[1], skiplines_begin=1,
      header_exists=false, colnames=[:id, :ident, :orig, :dest,
      :aircraft, :time, :lat, :lon, :speed, :alt, :climb, :heading, :direction,
      :facility, :description, :est]))
    # Calculate time from individual columns and add as DateTime to DataFrame
    flights.time = ZonedDateTime.(DateTime.(flights.time, "mm/dd/yyyy HH:MM:SS"), tz.tz"UTC")

    # Initialise loop over file
    global iStart = 1
    # Loop over all data points
    for i = 1:length(flights.time)-1
      if flights.id[i] ≠ flights.id[i+1]
        # When flight ID changes save data as flightData and set start index to next dataset
        push!(archive, flightData(flights.time[iStart:i], flights.lat[iStart:i],
          flights.lon[iStart:i], flights.alt[iStart:i], flights.heading[iStart:i],
          flights.climb[iStart:i], flights.speed[iStart:i]))
        global iStart = i+1
      end
    end
    # Save last flight of the dataset
    i = length(flights.time)
    push!(archive, flightData(flights.time[iStart:i], flights.lat[iStart:i],
      flights.lon[iStart:i], flights.alt[iStart:i], flights.heading[iStart:i],
      flights.climb[iStart:i], flights.speed[iStart:i]))
  end

  return archive
end #function loadArchive
