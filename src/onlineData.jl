"""
    loadOnlineData(files)

documentation
"""
function loadOnlineData(files)
  # Initialise inventory file array
  archive = Vector{flightData}(undef,length(files))
  for (n, file) in enumerate(files)
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
    # Retrieve date from filename
    m = match(r"_(.*?)_", basename(file))
    date = Dates.Date(m.captures[1], "d-u-y", locale="english")
    # Set to 2 days prior to allow corrections for timezone diffences in the next step
    date -= Dates.Day(2)
    # Delete rows with invalid times
    flight = flight[[length(t) == 15 for t in flight.time],:]

    ### Convert times to datetime and extract heading and climbing rate as Int32
    # Initialise time and date vectors
    flightdate = Vector{Union{Missing,Date}}(undef, length(flight.time))
    flighttime = Vector{Union{Missing,Time}}(undef, length(flight.time))
    # Initialise vectors for altitude, heading and climb to convert from strings to Int32
    heading = Vector{Union{Missing,Int32}}(undef, length(flight.time))
    climb = Vector{Union{Missing,Int32}}(undef, length(flight.time))
    altitude = Vector{Union{Missing,Int32}}(undef, length(flight.time))
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
        parse(Int32, flight.alt[i][findall([isdigit.(i) for i ∈ flight.alt[i]])])
      catch; 0
      end
      climb[i] = try
        parse(Int32, flight.climb[i][findall([isdigit.(i) | (i == '-') for i ∈ flight.climb[i]])])
      catch; 0
      end
      try
        idx = [j for j in eachindex(flight.heading) if isdigit(flight.heading[j])]
        heading[i] = parse.(Int32, flight.heading[idx])
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
    # Save data as flightData
    archive[n] = flightData(flight.time, flight.lat, flight.lon, flight.alt,
      flight.heading, flight.climb, flight.speed)
  end #loop over files

  return archive
end #function loadOnlineData


"""
    findtextfiles(inventory::Vector{String}, folder::String) -> inventory

Load `inventory` files saved in `folder` to the `inventory` holding the file names
and locations as a vector of strings.
"""
function findtextfiles(inventory::Vector{String}, folder::String)
  # Scan directory for files and folders and save directory
  dir = readdir(folder); path = abspath(folder)
  for file in dir
    cwd = joinpath(path, file)
    if endswith(file, ".txt") || endswith(file, ".dat")
      push!(inventory, cwd)
    elseif isdir(cwd)
      inventory = findtextfiles(inventory, cwd)
    end
  end

  return inventory
end # function findtextfiles
