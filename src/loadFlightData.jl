### Routines related to loading FlightData

"""
    loadInventory(files::Vector{String}; altmin=15_000, filterCloudfree::Bool=true) -> Vector{FlightData}

From a list of `files`, return a `Vector{FlightData}` that can
be saved to the `inventory` field in `FlightDB`.

When the `Vector{FlightData}` is constructed, data can be filtered by a minimum
altitude threshold of the aircraft data (default: `altmin=15_000`) and by the
existance of cirrus clouds at flight level (default: `filterCloudfree=true`;
currently only place holder, still needs to be implemented).
"""
function loadInventory(files::Vector{String}; altmin=15_000, filterCloudfree::Bool=true)

  # Initialise inventory file array and start MATLAB for PCHIP fitting
  inventory = FlightData[]

  # Loop over files
  prog = pm.Progress(2length(files), "load inventory...")
  for file in files

    # Load data
    parallel = VERSION ≥ v"1.3" ? true : false
    flights = CSV.read(file, datarow=3, footerskip=2, ignoreemptylines=true,
      silencewarnings=true, threaded=parallel, dateformat="HH:MM:SS.sssm")
    # Monitor progress for progress bar
    pm.next!(prog, showvalues = [(:file,splitext(basename(file))[1])])

    # Calculate time from individual columns and add as DateTime to DataFrame
    flights.time = [ZonedDateTime(flights.SEGMENT_YEAR[i], flights.SEGMENT_MONTH[i],
      flights.SEGMENT_DAY[i], flights.SEGMENT_HOUR[i], flights.SEGMENT_MIN[i],
      flights.SEGMENT_SEC[i], tz.tz"UTC") for i = 1:length(flights.FLIGHT_ID)]

    # Initialise loop over file
    FID = flights.FLIGHT_ID[1]
    lat = Float64[]; lon = Float64[];
    alt = Float64[]; t = ZonedDateTime[]; speed = Float64[]

    # Loop over all data points
    for i = 1:length(flights.time)
      # If the next flight ID is found, save current flight
      if flights.FLIGHT_ID[i] ≠ FID || i == length(flights.time)
        # Ignore data with less than 2 data points
        if length(t) ≤ 1  FID = flights.FLIGHT_ID[i]; continue  end
        # calculate area covered by flight
        lp = any(lon .> 0) ? maximum(lon[lon.≥0]) - minimum(lon[lon.≥0]) : 0
        ln = any(lon .< 0) ? maximum(lon[lon.<0]) - minimum(lon[lon.<0]) : 0
        # Determine main direction of flight (N<>S, E<>W) and use it as x values
        # for flight interpolation (info stored as bool useLON)
        useLON = maximum(lat) - minimum(lat) ≤ (lp + ln) * cosd(stats.mean(lat)) ? true : false
        useLON ? (x = lon; y = lat) : (x = lat; y = lon)
        # Remove duplicate points in data
        x, y, alt, speed, t = remdup(x, y, alt, speed, t)
        # find flex points to cut data in segments needed for the interpolation
        flex = findFlex(x)
        # Save the FlightData in the inventory vector
        push!(inventory, FlightData(t, lat, lon, alt, [missing for i = 1:length(t)],
          [missing for i = 1:length(t)], speed, FID, missing,
          missing, missing, flex, useLON, "VOLPE AEDT", file))
        # Empty data vectors
        lat = Float64[]; lon = Float64[];
        alt = Float64[]; t = ZonedDateTime[]; speed = Float64[]
        FID = flights.FLIGHT_ID[i]
        # segtime = DateTime[]; segdist = Float64[]
      end
      # Filter altitude threshold
      if flights.ALTITUDE[i] ≥ altmin
        push!(lat, flights.LATITUDE[i]); push!(lon, flights.LONGITUDE[i])
        push!(alt, flights.ALTITUDE[i]); push!(speed, flights.SPEED[i])
        push!(t, flights.time[i])
      end
    end #loop over flights
    # Monitor progress for progress bar
    pm.next!(prog, showvalues = [(:file,splitext(basename(file))[1])])
  end #loop over files
  pm.finish!(prog)

  return inventory
end #function loadInventory


"""
    loadArchive(files::Vector{String}; altmin::Int=15_000, filterCloudfree::Bool=true) -> Vector{FlightData}

From a list of `files`, return a `Vector{FlightData}` that can
be saved to the `archive` field in `FlightDB`.

When the `Vector{FlightData}` is constructed, data can be filtered by a minimum
altitude threshold of the aircraft data (default: `altmin=15_000`) and by the
existance of cirrus clouds at flight level (default: `filterCloudfree=true`;
currently only place holder, still needs to be implemented).
"""
function loadArchive(files::Vector{String}; altmin::Int=15_000, filterCloudfree::Bool=true)
  # Initialise archive file array
  archive = FlightData[]
  # Loop over database files
  prog = pm.Progress(length(files), "load archive...")
  for file in files
    # Load data
    parallel = VERSION ≥ v"1.3" ? true : false
    flights = CSV.read(file, datarow=2, normalizenames=true, ignoreemptylines=true,
      silencewarnings=true, threaded=parallel, copycols=true, dateformat="m/d/y H:M:S",
      types = Dict(:Altitude_feet_ => Float64, :Groundspeed_knots_ => Float64))
    # Calculate time from individual columns and add as DateTime to DataFrame
    flights.Time_UTC_ = ZonedDateTime.(flights.Time_UTC_, tz.tz"UTC")

    # Initialise loop over file
    FID = flights.Flight_ID[1]; n = 1
    lat = Float64[]; lon = Float64[]
    alt = Float64[]; t = ZonedDateTime[]; speed = Union{Missing,Float64}[]
    climb = Union{Missing,Int}[]; head = Union{Missing,Int}[]

    # Initialise loop over file
    # Loop over all data points
    for i = 1:length(flights.Time_UTC_)
      # Save flight, if flight ID changes
      if flights.Flight_ID[i] ≠ FID || i == length(flights.Time_UTC_)
        # Ignore data with less than 2 data points
        if length(t) ≤ 1
          n = i
          FID = flights.Flight_ID[n]
          continue
        end
        # calculate area covered by flight
        lp = any(lon .> 0) ? maximum(lon[lon.≥0]) - minimum(lon[lon.≥0]) : 0
        ln = any(lon .< 0) ? maximum(lon[lon.<0]) - minimum(lon[lon.<0]) : 0
        # Determine main direction of flight (N<>S, E<>W) and use it as x values
        # for flight interpolation (info stored as bool useLON)
        useLON = maximum(lat) - minimum(lat) ≤ (lp + ln) * cosd(stats.mean(lat)) ? true : false
        # find flex points to cut data in segments needed for the interpolation
        flex = useLON ? findFlex(lon) : findFlex(lat)
        # Save the FlightData in the archive vector
        push!(archive, FlightData(t, lat, lon, alt, head, climb, speed,
        FID, flights.Ident[n], flights.Aircraft_Type[n],
        (orig=flights.Origin[n], dest=flights.Destination[n]), flex, useLON,
        "FlightAware", file))

        # Reset temporary data arrays
        lat = Float64[]; lon = Float64[]
        alt = Float64[]; t = ZonedDateTime[]; speed = Union{Missing,Float64}[]
        climb = Union{Missing,Int}[]; head = Union{Missing,Int}[]
        # Set Flight ID and position to next flight
        n = i
        FID = flights.Flight_ID[n]
      end
      # Filter data
      if !ismissing(flights.Latitude[i]) && !ismissing(flights.Longitude[i]) &&
        !ismissing(flights.Altitude_feet_[i]) && flights.Altitude_feet_[i] ≥ altmin
        push!(t, flights.Time_UTC_[i])
        push!(lat, flights.Latitude[i]); push!(lon, flights.Longitude[i])
        push!(alt, flights.Altitude_feet_[i]); push!(speed, flights.Groundspeed_knots_[i])
        push!(climb, flights.Rate[i]); push!(head, flights.Course[i])
      end
      # Monitor progress for progress bar
      pm.next!(prog, showvalues = [(:file,splitext(basename(file))[1])])
    end #loop over files
    pm.finish!(prog)

    return archive
  end

  return archive
end #function loadArchive


"""
    loadOnlineData(files::Vector{String}; altmin::Int=15_000, filterCloudfree::Bool=true,
      delim::Union{Nothing,Char,String}=nothing) -> Vector{FlightData}

From a list of `files`, return a `Vector{FlightData}` that can
be saved to the `onlineData` field in `FlightDB`.

The delimiter of the data in the input file can be specified by a string or character.
Default is `nothing`, which means auto-detection of the delimiter is used.

When the `Vector{FlightData}` is constructed, data can be filtered by a minimum
altitude threshold of the aircraft data (default: `altmin=15_000`) and by the
existance of cirrus clouds at flight level (default: `filterCloudfree=true`;
currently only place holder, still needs to be implemented).
"""
function loadOnlineData(files::Vector{String}; altmin::Int=15_000, filterCloudfree::Bool=true,
  delim::Union{Nothing,Char,String}=nothing)
  # Initialise inventory file array
  archive = FlightData[]
  # Loop over files with online data
  prog = pm.Progress(length(files), "load online data...")
  for (n, file) in enumerate(files)
    # Read flight data
    parallel = VERSION ≥ v"1.3" ? true : false
    flight = CSV.read(file, delim=delim, ignoreemptylines=true, normalizenames=true, copycols=true,
      silencewarnings=true, threaded=parallel, types=Dict(:Latitude => Float64,
      :Longitude => Float64, :feet => String, :kts => Float64, :Course => String,
      :Rate => String))
    if df.names(flight)[2:9] ≠ [:Latitude, :Longitude, :Course, :kts, :mph, :feet, :Rate, :Reporting_Facility]
      println()
      println()
      @warn "Unknown file format in $file.\nData skipped."
      continue
    end

    ### Get timezone from input data or use local time for undefined timezones
    # Define timezones as UTC offset to avoid conflicts during
    # changes to/from daylight saving

    # Time is the first column and has to be addressed as flight[!,1] in the code
    # due to different column names, in which the timezone is included
    if occursin("_CET_", string(names(flight)[1]))
      timezone = tz.tz"+0100"
    elseif occursin("_CEST_", string(names(flight)[1]))
      timezone = tz.tz"+0200"
    else
      timezone = tz.localzone()
    end
    # Retrieve date and metadata from filename
    filename = splitext(basename(file))[1]
    flightID, datestr, course = try match(r"(.*?)_(.*?)_(.*)", filename).captures
    catch
      println()
      println()
      @warn "Flight ID, date, and course not found in $file. Data skipped."
      continue
    end
    orig, dest = match(r"(.*)[-|_](.*)", course).captures
    date = try Dates.Date(datestr, "d-u-y", locale="english")
    catch
      println()
      println()
      @warn "Unable to parse date in $file. Data skipped."
      continue
    end
    # Set to 2 days prior to allow corrections for timezone diffences in the next step
    date -= Dates.Day(2)
    ### Convert times to datetime and extract heading and climbing rate as Int32
    # Initialise time vector
    flighttime = ZonedDateTime[]
    # Initialise vectors for altitude, heading and climb to convert from strings to Int32
    altitude = Union{Missing,Float64}[]
    heading = Union{Missing,Int}[]
    climbingrate = Union{Missing,Int}[]
    # Loop over times
    for i=length(flight[!,1]):-1:1
      alt = try parse(Float64, join([n for n in flight.feet[i] if isnumeric(n)]))
      catch; missing;  end
      climb = try parse(Int, join([n for n in flight.Rate[i] if isnumeric(n) || n == '-']))
      catch; missing;  end
      head = try parse(Int, join([n for n in flight.Course[i] if isnumeric(n)]))
      catch; missing;  end
      if length(flight[i,1]) ≠ 15 || ismissing(flight.Latitude[i]) ||
          ismissing(flight.Longitude[i]) || ismissing(alt) || alt < altmin
        df.deleterows!(flight, i)
        continue
      end
      # Derive date from day of week and filename
      while Dates.dayabbr(date) ≠ flight[i,1][1:3]
        date += Dates.Day(1)
      end
      # Derive time from time string
      t = if VERSION ≥ v"1.3"
        # Use AM/PM format for Julia > version 1.3
        Time(flight[i,1][5:end], "I:M:S p")
      else
        # Calculate time manually otherwise
        t = Time(flight[i,1][5:12], "H:M:S")
        if flight[i,1][end-1:end] == "PM" && !(Dates.hour(t)==12)
          t += Dates.Hour(12)
        elseif flight[i,1][end-1:end] == "AM" && Dates.hour(t)==12
          t -= Dates.Hour(12)
        else
          t
        end
      end
      # Save data that needed tweaking of current time step
      if VERSION ≥ v"1.1"
        push!(flighttime, ZonedDateTime(DateTime(date, t), timezone))
      else
        push!(flighttime, ZonedDateTime(DateTime(Dates.yearmonthday(date)...,
          Dates.hour(t), Dates.minute(t), Dates.second(t)), timezone))
      end
      push!(altitude, alt); push!(climbingrate, climb); push!(heading, head)
    end #loop of flight

    # Skip data with all data points below the altitude threshold
    isempty(altitude) && continue
    # Save revised data to DataFrame
    flight.time = flighttime
    flight.feet = altitude
    flight.Rate = climbingrate
    flight.Course = heading
    # Convert lat/lon from type Union{Missing,Float64} to Float64 to be processed by findFlex
    flight.Latitude = float.(flight.Latitude); flight.Longitude = float.(flight.Longitude)

    # Determine, whether to use lat or lon as x data and find flex points
    lp = any(flight.Longitude .> 0) ?
      maximum(flight.Longitude[flight.Longitude.≥0]) - minimum(flight.Longitude[flight.Longitude.≥0]) : 0
    ln = any(flight.Longitude .< 0) ?
      maximum(flight.Longitude[flight.Longitude.<0]) - minimum(flight.Longitude[flight.Longitude.<0]) : 0
    useLON = maximum(flight.Latitude) - minimum(flight.Latitude) ≤ (lp + ln) *
      cosd(stats.mean(flight.Latitude)) ? true : false
    flex = useLON ? findFlex(flight.Longitude) : findFlex(flight.Latitude)

    # Save data as FlightData
    push!(archive, FlightData(flight.time, flight.Latitude, flight.Longitude,
      flight.feet, flight.Course, flight.Rate, flight.kts, n, flightID, missing,
      (orig=orig, dest=dest), flex, useLON, "flightaware.com", file))
    # Monitor progress for progress bar
    pm.next!(prog, showvalues = [(:file,filename)])
  end #loop over files
  pm.finish!(prog)

  return archive
end #function loadOnlineData
