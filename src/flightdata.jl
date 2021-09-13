### Routines related to loading FlightTrack

"""
    loadVOLPE(files::String...; Float::DataType=Float32, altmin::Real=5_000)
      -> Vector{FlightTrack}

From a list of `files`, return a `Vector{FlightTrack}` that can
be saved to the `inventory` field in `FlightSet`.

When the `Vector{FlightTrack{T}}` is constructed, data can be filtered by a minimum
altitude threshold in meters of the aircraft data (default: `altmin=5_000`).
Floating point numbers in `FlightTrack` are of the precision set by `Float`,
by default `Float32`.
"""
function loadVOLPE(
  files::String...;
  Float::DataType=Float32,
  altmin::Real=5_000
)

  # Initialise inventory file array and start MATLAB for PCHIP fitting
  inventory = FlightData{Float}[]

  # Loop over files
  prog = pm.Progress(2length(files), "load inventory...")
  for file in files
    # Load data
    flights = CSV.read(file, DataFrame, datarow=3, footerskip=2, ignoreemptylines=true,
      silencewarnings=true, threaded=true, copycols = false, dateformat="HH:MM:SS.sssm",
      types = Dict("LATITUDE" => Float, "LONGITUDE" => Float, "ALTITUDE" => Float,
      "SPEED" => Float), select = [1:11;16])
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
    track = DataFrame(time = DateTime[], lat = Float[], lon = Float[],
      alt = Float[], speed = Float[])
    FID = flights.FLIGHT_ID[1]

    # Loop over all data points
    for i = 1:length(flights.time)
      # If the next flight ID is found, save current flight
      if flights.FLIGHT_ID[i] ≠ FID || i == length(flights.time)
        # Ignore data with less than 2 data points
        if length(track.time) ≤ 1
          FID = flights.FLIGHT_ID[i]
          # Empty possible entry
          track = DataFrame(time = DateTime[], lat = Float[], lon = Float[],
            alt = Float[], speed = Float[])
          continue
        end
        # Determine predominant flight direction, inflection points, and remove duplicate entries
        flex, useLON = preptrack!(track)
        # Save the FlightTrack in the inventory vector
        isempty(flex) || push!(inventory, FlightTrack{Float}(track, FID,
          missing, missing, missing, flex, useLON, "VOLPE", file))
        # Empty data vectors
        track = DataFrame(time = DateTime[], lat = Float[], lon = Float[],
          alt = Float[], speed = Float[])
        FID = flights.FLIGHT_ID[i]
      end # saving FlightTrack
      # Filter altitude threshold
      if ismissing(flights.ALTITUDE[i]) || flights.ALTITUDE[i] ≥ altmin
        push!(track, [flights.time[i], flights.LATITUDE[i],
          flights.LONGITUDE[i], flights.ALTITUDE[i], flights.SPEED[i]])
      end
    end #loop over flights
    # Monitor progress for progress bar
    pm.next!(prog, showvalues = [(:file,splitext(basename(file))[1])])
  end #loop over files
  pm.finish!(prog)

  return inventory
end #function loadVOLPE


"""
    loadFA(files::String...; Float::DataType=Float32, altmin::Real=5_000)
      -> Vector{FlightTrack}

From a list of `files`, return a `Vector{FlightTrack}` that can
be saved to the `archive` field in `FlightSet`.

When the `Vector{FlightTrack{T}}` is constructed, data can be filtered by a minimum
altitude threshold in meters of the aircraft data (default: `altmin=5_000`).
Floating point numbers in `FlightTrack` are of the precision set by `Float`,
by default `Float32`.
"""
function loadFA(
  files::String...;
  Float::DataType=Float32,
  altmin::Real=5_000
)
  # Initialise archive file array
  archive = FlightData{Float}[]
  # Loop over database files
  prog = pm.Progress(length(files), "load archive...")
  for file in files
    # Load data from csv file into standardised DataFrame
    flights = readArchive(file, Float)
    # Unit conversions
    flights.alt = ft2m.(flights.alt)
    flights.speed = knot2mps.(flights.speed)
    flights.climb = ftpmin2mps.(flights.climb)


    # Initialise loop over file
    FID = flights.dbID[1]; n = 1
    track = DataFrame(time = DateTime[], lat = Float[]; lon = Float[],
      alt = Union{Missing,Float}[], speed = Union{Missing,Float}[],
      climb = Union{Missing,Float}[], heading = Union{Missing,Int}[])

    # Loop over all data points
    for i = 1:length(flights.time)
      # Save flight, if flight ID changes
      if flights.dbID[i] ≠ FID || i == length(flights.time)
        # Ignore data with less than 2 data points
        if length(track.time) ≤ 1
          n = i
          FID = flights.dbID[n]
          # Empty possible entry
          track = DataFrame(time = DateTime[], lat = Float[]; lon = Float[],
            alt = Union{Missing,Float}[], speed = Union{Missing,Float}[],
            climb = Union{Missing,Float}[], heading = Union{Missing,Int}[])
          continue
        end
        # Determine predominant flight direction, inflection points, and remove duplicate entries
        flex, useLON = preptrack!(track)
        # Save the FlightTrack in the archive vector
        isempty(flex) || push!(archive, FlightTrack{Float}(track, FID, flights.flightID[n],
          flights.type[n], (orig=flights.orig[n], dest=flights.dest[n]), flex, useLON,
          "FlightAware", file))

        # Reset temporary data arrays
        track = DataFrame(time = DateTime[], lat = Float[]; lon = Float[],
          alt = Union{Missing,Float}[], speed = Union{Missing,Float}[],
          climb = Union{Missing,Float}[], heading = Union{Missing,Int}[])
        # Set Flight ID and position to next flight
        n = i
        FID = flights.dbID[n]
      end
      # Filter data
      if !ismissing(flights.lat[i]) && !ismissing(flights.lon[i]) &&
        (ismissing(flights.alt[i]) || flights.alt[i] ≥ altmin)
        push!(track, [flights.time[i], flights.lat[i],
          flights.lon[i], flights.alt[i], flights.speed[i],
          flights.climb[i], flights.heading[i]])
      end
    end #loop over current flight
    # Monitor progress for progress bar
    pm.next!(prog, showvalues = [(:file,splitext(basename(file))[1])])
  end #loop over files
  pm.finish!(prog)

  return archive
end #function loadFA


"""
    loadWD(files::String...; Float::DataType=Float32, altmin::Real=5_000, delim::Union{Nothing,Char,String}=nothing)
      -> Vector{FlightTrack}

From a list of `files`, return a `Vector{FlightTrack}` that can
be saved to the `onlineData` field in `FlightSet`.

The `delim`iter of the data in the input file can be specified by a string or character.
Default is `nothing`, which means auto-detection of the delimiter is used.

When the `Vector{FlightTrack{T}}` is constructed, data can be filtered by a minimum
altitude threshold in meters of the aircraft data (default: `altmin=5_000`).
Floating point numbers in `FlightTrack` are of the precision set by `Float`,
by default `Float32`.
"""
function loadWD(
  files::String...;
  Float::DataType=Float32,
  altmin::Real=5_000,
  delim::Union{Nothing,Char,String}=nothing
)
  # Initialise inventory file array
  archive = FlightData{Float}[]
  # Loop over files with online data
  prog = pm.Progress(length(files), "load online data...")
  for file in files
    # Read flight data
    flight = CSV.read(file, DataFrame, delim=delim, ignoreemptylines=true, normalizenames=true,
      silencewarnings=true, threaded=false, copycols=false, types=Dict(:Latitude => Float,
      :Longitude => Float, :feet => String, :kts => Float, :Course => String,
      :Rate => String), drop = ["mph", "Reporting_Facility"])
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

    # Get time zone, data and flight metadata from file name and header of time column
    filename = splitext(basename(file))[1]
    date, timezone, flightID, orig, dest = get_DateTimeRoute(filename, tzone)

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
      t = Time(flight.time[i][5:end], "I:M:S p")
      # Save data that needed tweaking of current time step
      pushfirst!(flighttime, ZonedDateTime(DateTime(date, t), timezone))
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
    flex, useLON = preptrack!(flight)

    # Save data as FlightTrack
    isempty(flex) || push!(archive, FlightTrack{Float}(flight, replace(filename, "_" => "/"),
      flightID, missing, (orig=orig, dest=dest), flex, useLON, "web", file))
    # Monitor progress for progress bar
    pm.next!(prog)
  end #loop over files
  pm.finish!(prog)

  return archive
end #function loadWD


"""
    readArchive(file, Float=Float32) -> DataFrame

Read FlightAware archived data from a csv `file` and return content as DataFrame.

The routine works for several FlightAware archive versions. Floating point numbers
are read with single precision or as defined by kwarg `Float`.
"""
function readArchive(file, Float=Float32)
  # Read header and check version
  header = readline(file)
  old = contains(header, "(kts)")

  # Load correct data version and return as DataFrame
  old ?
    CSV.read(file, DataFrame, datarow=2, normalizenames=true, ignoreemptyrows=true,
    silencewarnings=true, threaded=true, dateformat="m/d/y H:M:S",
    header = ["dbID", "flightID", "type", "orig", "dest", "time", "lat", "lon",
      "speed", "alt", "climb", "heading", "direction", "estimated"],
    types = Dict(:lat => Float, :lon => Float, :alt => Float, :speed => Float, :climb => Float),
    drop = ["direction", "estimated"]) :
    CSV.read(file, DataFrame, datarow=2, normalizenames=true, ignoreemptyrows=true,
    silencewarnings=true, threaded=true, dateformat="m/d/y H:M:S",
    header = ["dbID", "flightID", "type", "orig", "dest", "time", "lat", "lon",
      "speed", "alt", "climb", "heading", "direction", "fac_name", "fac_descr", "estimated"],
    types = Dict(:lat => Float, :lon => Float, :alt => Float, :speed => Float, :climb => Float),
    drop = ["direction", "fac_name", "fac_descr", "estimated"])
end #function readArchive
