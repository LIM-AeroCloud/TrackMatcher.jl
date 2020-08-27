### Routines related to loading FlightData

"""
    loadInventory(files::String...; altmin::Real=5_000) -> Vector{FlightData}

From a list of `files`, return a `Vector{FlightData}` that can
be saved to the `inventory` field in `FlightDB`.

When the `Vector{FlightData}` is constructed, data can be filtered by a minimum
altitude threshold in meters of the aircraft data (default: `altmin=5_000`).
"""
function loadInventory(files::String...; Float::DataType=Float32, altmin::Real=5_000)

  # Initialise inventory file array and start MATLAB for PCHIP fitting
  inventory = FlightData[]
  parallel = VERSION ≥ v"1.3" ? true : false

  # Loop over files
  prog = pm.Progress(2length(files), "load inventory...")
  for file in files
    # Load data
    flights = CSV.File(file, datarow=3, footerskip=2, ignoreemptylines=true,
      silencewarnings=true, threaded=parallel, dateformat="HH:MM:SS.sssm",
      types = Dict("LATITUDE" => Float, "LONGITUDE" => Float, "ALTITUDE" => Float,
      "SPEED" => Float), select = [1:11;16]) |> df.DataFrame!
    # Monitor progress for progress bar
    pm.next!(prog, showvalues = [(:file,splitext(basename(file))[1])])

    # Calculate time from individual columns and add as DateTime to DataFrame
    flights.time = [DateTime(flights.SEGMENT_YEAR[i], flights.SEGMENT_MONTH[i],
      flights.SEGMENT_DAY[i], flights.SEGMENT_HOUR[i], flights.SEGMENT_MIN[i],
      flights.SEGMENT_SEC[i]) for i = 1:length(flights.FLIGHT_ID)]
    # Unit conversions
    flights.ALTITUDE = ft2m.(flights.ALTITUDE)
    flights.SPEED = knot2mps.(flights.SPEED)

    # Initialise loop over file
    flightdata = DataFrame(time = DateTime[], lat = Float[], lon = Float[],
      alt = Float[], speed = Float[])
    FID = flights.FLIGHT_ID[1]

    # Loop over all data points
    for i = 1:length(flights.time)
      # If the next flight ID is found, save current flight
      if flights.FLIGHT_ID[i] ≠ FID || i == length(flights.time)
        # Ignore data with less than 2 data points
        if length(flightdata.time) ≤ 1
          FID = flights.FLIGHT_ID[i]
          # Empty possible entry
          flightdata = DataFrame(time = DateTime[], lat = Float[], lon = Float[],
            alt = Float[], speed = Float[])
          continue
        end
        # Determine predominant flight direction, inflection points, and remove duplicate entries
        flightdata, flex, useLON = preptrack(flightdata)
        # Save the FlightData in the inventory vector
        standardisecols!(flightdata)
        push!(inventory, FlightData(flightdata, FID,
          missing, missing, missing, flex, useLON, "VOLPE AEDT", file))
        # Empty data vectors
        flightdata = DataFrame(time = DateTime[], lat = Float[], lon = Float[],
          alt = Float[], speed = Float[])
        FID = flights.FLIGHT_ID[i]
      end # saving FlightData
      # Filter altitude threshold
      if ismissing(flights.ALTITUDE[i]) || flights.ALTITUDE[i] ≥ altmin
        push!(flightdata, [flights.time[i], flights.LATITUDE[i],
          flights.LONGITUDE[i], flights.ALTITUDE[i], flights.SPEED[i]])
      end
    end #loop over flights
    # Monitor progress for progress bar
    pm.next!(prog, showvalues = [(:file,splitext(basename(file))[1])])
  end #loop over files
  pm.finish!(prog)

  return inventory
end #function loadInventory


"""
    loadArchive(files::String...; altmin::Real=5_000)
      -> Vector{FlightData}

From a list of `files`, return a `Vector{FlightData}` that can
be saved to the `archive` field in `FlightDB`.

When the `Vector{FlightData}` is constructed, data can be filtered by a minimum
altitude threshold in meters of the aircraft data (default: `altmin=5_000`).
"""
function loadArchive(files::String...; Float::DataType=Float32, altmin::Real=5_000)
  # Initialise archive file array
  archive = FlightData[]
  parallel = VERSION ≥ v"1.3" ? true : false
  # Loop over database files
  prog = pm.Progress(length(files), "load archive...")
  for file in files
    # Load data from csv file into standardised DataFrame
    flights = readArchive(file, parallel, Float)
    # Unit conversions
    flights.alt = ft2m.(flights.alt)
    flights.speed = knot2mps.(flights.speed)
    flights.climb = ftpmin2mps.(flights.climb)


    # Initialise loop over file
    FID = flights.dbID[1]; n = 1
    flightdata = DataFrame(time = DateTime[], lat = Float[]; lon = Float[],
      alt = Union{Missing,Float}[], speed = Union{Missing,Float}[],
      climb = Union{Missing,Float}[], heading = Union{Missing,Int}[])

    # Loop over all data points
    for i = 1:length(flights.time)
      # Save flight, if flight ID changes
      if flights.dbID[i] ≠ FID || i == length(flights.time)
        # Ignore data with less than 2 data points
        if length(flightdata.time) ≤ 1
          n = i
          FID = flights.dbID[n]
          # Empty possible entry
          flightdata = DataFrame(time = DateTime[], lat = Float[]; lon = Float[],
            alt = Union{Missing,Float}[], speed = Union{Missing,Float}[],
            climb = Union{Missing,Float}[], heading = Union{Missing,Int}[])
          continue
        end
        # Determine predominant flight direction, inflection points, and remove duplicate entries
        flightdata, flex, useLON = preptrack(flightdata)
        # Save the FlightData in the archive vector
        standardisecols!(flightdata)
        push!(archive, FlightData(flightdata, FID, flights.flightID[n],
          flights.type[n], (orig=flights.orig[n],
          dest=flights.dest[n]), flex, useLON, "FlightAware", file))

        # Reset temporary data arrays
        flightdata = DataFrame(time = DateTime[], lat = Float[]; lon = Float[],
          alt = Union{Missing,Float}[], speed = Union{Missing,Float}[],
          climb = Union{Missing,Float}[], heading = Union{Missing,Int}[])
        # Set Flight ID and position to next flight
        n = i
        FID = flights.dbID[n]
      end
      # Filter data
      if !ismissing(flights.lat[i]) && !ismissing(flights.lon[i]) &&
        (ismissing(flights.alt[i]) || flights.alt[i] ≥ altmin)
        push!(flightdata, [flights.time[i], flights.lat[i],
          flights.lon[i], flights.alt[i], flights.speed[i],
          flights.climb[i], flights.heading[i]])
      end
    end #loop over current flight
    # Monitor progress for progress bar
    pm.next!(prog, showvalues = [(:file,splitext(basename(file))[1])])
  end #loop over files
  pm.finish!(prog)

  return archive
end #function loadArchive


"""
    loadOnlineData(files::String...; altmin::Real=5_000, delim::Union{Nothing,Char,String}=nothing)
      -> Vector{FlightData}

From a list of `files`, return a `Vector{FlightData}` that can
be saved to the `onlineData` field in `FlightDB`.

The `delim`iter of the data in the input file can be specified by a string or character.
Default is `nothing`, which means auto-detection of the delimiter is used.

When the `Vector{FlightData}` is constructed, data can be filtered by a minimum
altitude threshold in meters of the aircraft data (default: `altmin=5_000`).
"""
function loadOnlineData(files::String...; Float::DataType=Float32, altmin::Real=5_000,
  delim::Union{Nothing,Char,String}=nothing)
  # Initialise inventory file array
  # parallel = VERSION ≥ v"1.3" ? true : false
  archive = FlightData[]
  # Loop over files with online data
  prog = pm.Progress(length(files), "load online data...")
  for file in files
    # Read flight data
    flight = CSV.File(file, delim=delim, ignoreemptylines=true, normalizenames=true,
      silencewarnings=true, threaded=false, types=Dict(:Latitude => Float,
      :Longitude => Float, :feet => String, :kts => Float, :Course => String,
      :Rate => String), drop = ["mph", "Reporting_Facility"]) |> df.DataFrame!
    # Convert knots to m/s
    flight.kts = knot2mps.(flight.kts)
    if length(names(flight)) ≠ 7 || names(flight)[2:7] ≠
      ["Latitude", "Longitude", "Course", "kts", "feet", "Rate"]
      println()
      println()
      @warn "Unknown file format.\nTry to specify column delimiter. Data skipped." file
      continue
    else
      tzone = string(names(flight)[1])
      df.rename!(flight, 1 => :time)
      df.rename!(flight, :Latitude => :lat, :Longitude => :lon, :Course => :heading,
        :kts => :speed, :feet => :alt, :Rate => :climb)
    end

    ### Get timezone from input data or use local time for undefined timezones
    # Define timezones as UTC offset to avoid conflicts during
    # changes to/from daylight saving

    # Time is the first column and has to be addressed as flight[!,1] in the code
    # due to different column names, in which the timezone is included
    timezone = if occursin("_CET_", tzone)
      tz.tz"+0100"
    elseif occursin("_CEST_", tzone)
      tz.tz"+0200"
    else
      tz.localzone()
    end
    # Retrieve date and metadata from filename
    filename = splitext(basename(file))[1]
    flightID, datestr, course = try match(r"(.*?)_(.*?)_(.*)", filename).captures
    catch
      println()
      println()
      @warn "Flight ID, date, and course not found. Data skipped." file
      continue
    end
    orig, dest = match(r"(.*)[-|_](.*)", course).captures
    date = try Dates.Date(datestr, "d-u-y", locale="english")
    catch
      println()
      println()
      @warn "Unable to parse date. Data skipped." file
      continue
    end
    # Set to 2 days prior to allow corrections for timezone diffences in the next step
    date += Dates.Day(2)
    ### Convert times to datetime and extract heading and climbing rate as Int
    # Initialise time vector
    flighttime = ZonedDateTime[]
    # Initialise vectors for altitude, heading and climb to convert from strings to Int
    altitude = Union{Missing,Float}[]
    heading = Union{Missing,Int}[]
    climbingrate = Union{Missing,Float}[]
    # Loop over times
    for i=length(flight[!,1]):-1:1
      alt = try ft2m(parse(Float, join([n for n in flight.alt[i]
        if isnumeric(n) || n == '.'])))
      catch; missing;  end
      climb = try ftpmin2mps(parse(Float, join([n for n in flight.climb[i]
        if isnumeric(n) || n == '.' || n == '-'])))
      catch; missing;  end
      head = try parse(Int, join([n for n in flight.heading[i]
        if isnumeric(n) || n == '.']))
      catch; missing;  end
      if length(flight.time[i]) ≠ 15 || ismissing(flight.lat[i]) ||
          ismissing(flight.lon[i]) || (!ismissing(alt) && alt < altmin)
        delete!(flight, i)
        continue
      end
      # Derive date from day of week and filename
      while Dates.dayabbr(date) ≠ flight.time[i][1:3]
        date -= Dates.Day(1)
      end
      # Derive time from time string
      t = if VERSION ≥ v"1.3"
        # Use AM/PM format for Julia > version 1.3
        Time(flight.time[i][5:end], "I:M:S p")
      else
        # Calculate time manually otherwise
        t = Time(flight.time[i][5:12], "H:M:S")
        if flight.time[i][end-1:end] == "PM" && !(Dates.hour(t)==12)
          t += Dates.Hour(12)
        elseif flight.time[i][end-1:end] == "AM" && Dates.hour(t)==12
          t -= Dates.Hour(12)
        else
          t
        end
      end
      # Save data that needed tweaking of current time step
      if VERSION ≥ v"1.1"
        pushfirst!(flighttime, ZonedDateTime(DateTime(date, t), timezone))
      else
        pushfirst!(flighttime, ZonedDateTime(DateTime(Dates.yearmonthday(date)...,
          Dates.hour(t), Dates.minute(t), Dates.second(t)), timezone))
      end
      pushfirst!(altitude, alt); pushfirst!(climbingrate, climb); pushfirst!(heading, head)
    end #loop of flight

    # Skip data with all data points below the altitude threshold or missing
    isempty(altitude) && continue
    # Save revised data to DataFrame
    flight.time = flighttime
    flight.alt = altitude
    flight.climb = climbingrate
    flight.heading = heading
    # Convert lat/lon from type Union{Missing,Float64} to Float64 to be processed by findflex
    flight.lat = float.(flight.lat); flight.lon = float.(flight.lon)

    # Determine predominant flight direction, inflection points, and remove duplicate entries
    flight, flex, useLON = preptrack(flight)

    # Save data as FlightData
    standardisecols!(flight)
    push!(archive, FlightData(flight, replace(filename, "_" => "/"), flightID,
      missing, (orig=orig, dest=dest), flex, useLON, "flightaware.com", file))
    # Monitor progress for progress bar
    pm.next!(prog)
  end #loop over files
  pm.finish!(prog)

  return archive
end #function loadOnlineData


"""
    readArchive(file, parallel) -> DataFrame

Read FlightAware archived data from a csv `file` and return content as DataFrame.

Using Julia version 1.3 or higher, `parallel` can be set to `true` to speed up file reading.
The routine works for several FlightAware archive versions.
"""
function readArchive(file, parallel, Float=Float32)
  # Read file
  flightdata = CSV.File(file, datarow=2, normalizenames=true, ignoreemptylines=true,
    silencewarnings=true, threaded=parallel, dateformat="m/d/y H:M:S",
    types = Dict(:Latitude => Float, :Longitude => Float, :Altitude_feet_ => Float,
    :Altitude_ft_ => Float, :Groundspeed_knots_ => Float, :Groundspeed_kts_ => Float,
    :Rate => Float),
    drop = ["Direction", "Facility_Name", "Facility_Description", "Estimated"])
  ## Get column order and define column names
  # Get column names of current file
  datacols = propertynames(flightdata[1])
  # Define unique key phrases for each column found in every FlightAware version
  keys = [datacols[i] for i in findfirst.(occursin.(str, string.(datacols)) for str in
    ["Flight_ID", "Ident", "Orig", "Dest", "Type",
    "Time", "Lat", "Lon", "Alt", "speed", "Rate", "Course"])]
  # Define standard column names
  colnames = ["dbID", "flightID", "orig", "dest", "type",
    "time", "lat", "lon", "alt", "speed", "climb", "heading"]
  # Construct DataFrame
  DataFrame(colnames .=> [getproperty(flightdata, key) for key in keys])
end #function readArchive
