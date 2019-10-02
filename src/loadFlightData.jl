"""
    loadInventory(files::Vector{String}) -> inventory

From a list of `files`, return an `inventory` as `Vector{FlightData}` that can
be saved to the `inventory` field in `FlightDB`.
"""
function loadInventory(files::Vector{String})
  # Initialise inventory file array
  inventory = FlightData[]
  @pm.showprogress 1 "load inventory..." for file in files
    # Load data
    flights = CSV.read(file, datarow=3, copycols=true)
    # Delete additional line at end of file
    df.deleterows!(flights, length(flights.FLIGHT_ID)-1:length(flights.FLIGHT_ID))
    # Calculate time from individual columns and add as DateTime to DataFrame
    flights.time = ZonedDateTime.(flights.SEGMENT_YEAR, flights.SEGMENT_MONTH,
      flights.SEGMENT_DAY, flights.SEGMENT_HOUR, flights.SEGMENT_MIN,
      flights.SEGMENT_SEC, tz.tz"UTC")

    # Initialise loop over file
    global iStart = 1
    # Loop over all data points
    for i = 1:length(flights.time)-1
      if flights.FLIGHT_ID[i] ≠ flights.FLIGHT_ID[i+1]
        # When flight ID changes save data as FlightData and set start index to next dataset
        push!(inventory, FlightData(flights.time[iStart:i], flights.LATITUDE[iStart:i],
          flights.LONGITUDE[iStart:i], flights.ALTITUDE[iStart:i],
          [missing for j=iStart:i], [missing for j=iStart:i], flights.SPEED[iStart:i],
          parse(Int, flights.FLIGHT_ID[i]), missing, missing, missing, file))
        global iStart = i+1
      end
    end
    # Save last flight of the dataset
    i = length(flights.time)
    push!(inventory, FlightData(flights.time[iStart:i], flights.LATITUDE[iStart:i],
      flights.LONGITUDE[iStart:i], flights.ALTITUDE[iStart:i],
      [missing for j=iStart:i], [missing for j=iStart:i], flights.SPEED[iStart:i],
      parse(Int, flights.FLIGHT_ID[i]), missing, missing, missing, file))
  end

  return inventory
end #function loadInventory


"""
    loadArchive(files::Vector{String}) -> archive

From a list of `files`, return an `archive` as `Vector{FlightData}` that can
be saved to the `archive` field in `FlightDB`.
"""
function loadArchive(files::Vector{String})
  # Initialise inventory file array
  archive = FlightData[]
  @pm.showprogress 1 "load archive..." for file in files
    # Load data
    flights = CSV.read(file, datarow=2, header=[:id, :ident, :orig, :dest,
      :aircraft, :time, :lat, :lon, :speed, :alt, :climb, :heading, :direction,
      :facility, :description, :est])
    flights.speed = convert.(Union{Missing,AbstractFloat}, flights.speed)
    # Calculate time from individual columns and add as DateTime to DataFrame
    flights.time = ZonedDateTime.(DateTime.(flights.time, "m/d/y H:M:S"), tz.tz"UTC")

    # Initialise loop over file
    global iStart = 1
    # Loop over all data points
    for i = 1:length(flights.time)-1
      if flights.id[i] ≠ flights.id[i+1]
        # When flight ID changes save data as FlightData and set start index to next dataset
        push!(archive, FlightData(flights.time[iStart:i], flights.lat[iStart:i],
          flights.lon[iStart:i], flights.alt[iStart:i],
          flights.heading[iStart:i], flights.climb[iStart:i],
          flights.speed[iStart:i], flights.id[i],
          flights.ident[i], flights.aircraft[i],
          (orig=flights.orig[i], dest=flights.dest[i]), file))
        global iStart = i+1
      end
    end
    # Save last flight of the dataset
    i = length(flights.time)
    push!(archive, FlightData(flights.time[iStart:i], flights.lat[iStart:i],
      flights.lon[iStart:i], flights.alt[iStart:i],
      flights.heading[iStart:i], flights.climb[iStart:i], flights.speed[iStart:i],
      flights.id[i], flights.ident[i], flights.aircraft[i],
      (orig=flights.orig[i], dest=flights.dest[i]), file))
  end

  return archive
end #function loadArchive


"""
    loadOnlineData(files::Vector{String}) -> archive

From a list of `files`, return an `archive` as `Vector{FlightData}` that can
be saved to the `onlineData` field in `FlightDB`.
"""
function loadOnlineData(files::Vector{String})
  # Initialise inventory file array
  archive = Vector{FlightData}(undef,length(files))
  @pm.showprogress 1 "load online data..." for (n, file) in enumerate(files)
    # Read flight data
    flight = CSV.read(file, delim='\t', datarow=3,
      header=["time", "lat", "lon", "heading", "speed", "mph", "alt", "climb", "facility"])
    ### Get timezone of input data
    # Set default timezone to local
    timezone = tz.localzone()
    # Infer timezone from table header
    m = match(r"\((.*?)\)", readline(file))
    zone = m.captures[1]
    # Define timezones as UTC offset to avoid conflicts during
    # changes to/from daylight saving
    if zone == "CET"
      timezone = tz.tz"+0100"
    elseif zone == "CEST"
      timezone = tz.tz"+0200"
    end
    # Retrieve date and metadata from filename
    filename = splitext(basename(file))[1]
    flightID, datestr, course = match(r"(.*?)_(.*?)_(.*)", filename).captures
    orig, dest = match(r"(.*)[-|_](.*)", course).captures
    date = Dates.Date(datestr, "d-u-y", locale="english")
    # Set to 2 days prior to allow corrections for timezone diffences in the next step
    date -= Dates.Day(2)
    # Delete rows with invalid times
    flight = flight[[length(t) == 15 for t in flight.time],:]

    ### Convert times to datetime and extract heading and climbing rate as Int32
    # Initialise time and date vectors
    flightdate = Vector{Union{Missing,Date}}(undef, length(flight.time))
    flighttime = Vector{Union{Missing,Time}}(undef, length(flight.time))
    # Initialise vectors for altitude, heading and climb to convert from strings to Int32
    heading = Vector{Union{Missing,Int}}(undef, length(flight.time))
    climb = Vector{Union{Missing,Int}}(undef, length(flight.time))
    altitude = Vector{Union{Missing,AbstractFloat}}(undef, length(flight.time))
    # Loop over times
    for i=1:length(flight.time)
      # Derive date from day of week and filename
      while Dates.dayabbr(date) ≠ flight.time[i][1:3]
        date += Dates.Day(1)
      end
      # Save date for current time step
      flightdate[i] = date
      # Derive time from time string
      if VERSION > v"1.3"
        # Use AM/PM format for Julia > version 1.3
        t = Time(flight.time[i][5:end], "I:M:S p")
      else
        # Calculate time manually otherwise
        t = Time(flight.time[i][5:12], "H:M:S")
        if flight.time[i][end-1:end] == "PM" && !(Dates.hour(t)==12 && Dates.minute(t)==0 &&
          Dates.second(t)==0)
          t += Dates.Hour(12)
        elseif flight.time[i][end-1:end] == "AM" && Dates.hour(t)==12 &&
          Dates.minute(t)==0 && Dates.second(t)==0
          t -= Dates.Hour(12)
        end
      end
      # Save time of current time step
      flighttime[i] = t
      # Filter altitude, heading and climbing rate for numbers and convert to Int32
      altitude[i] = try
        alt=parse(Int, flight.alt[i][findall([isdigit.(i) for i ∈ flight.alt[i]])])
        convert(AbstractFloat,alt)
      catch; 0.
      end
      climb[i] = try
        parse(Int, flight.climb[i][findall([isdigit.(i) | (i == '-') for i ∈ flight.climb[i]])])
      catch; 0
      end
      try
        idx = [j for j in eachindex(flight.heading) if isdigit(flight.heading[j])]
        heading[i] = parse.(Int, flight.heading[idx])
      catch
        heading[i] = missing
      end
    end
    # Save revised DateTime
    flight.time = ZonedDateTime.(DateTime.(flightdate, flighttime), timezone)
    # Save revised heading and climbing rate
    flight.climb = climb
    flight.heading = heading
    flight.alt = altitude
    # Save data as FlightData
    archive[n] = FlightData(flight.time, flight.lat, flight.lon, flight.alt,
      flight.heading, flight.climb, convert.(Union{Missing,AbstractFloat},flight.speed),
      n, flightID, missing, (orig=orig, dest=dest), file)
  end #loop over files

  return archive
end #function loadOnlineData
