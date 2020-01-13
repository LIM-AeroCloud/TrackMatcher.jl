"""
    loadFlightDB(DBtype::String, folder::Union{String, Vector{String}}...; remarks=nothing) -> struct `FlightDB`

Construct an instance of `FlightDB` with relevant flight data.


# DBtype

Specifies the database type; up to 3 types can be selected:
- `1` or `"i"`: Flight inventory with flight tracks saved in csv files and saved
  to property `inventory`.
- `2` or `"a"`: Commercially available flight track data by FlightAware saved as
  csv files and saved to property `archive`.
- `3` or `"o"`: Flight tracks available online from the FlightAware website in
  text (`.txt`/`.dat`) files and saved to property `onlineData`.


# folder

Directory holding the files of the databases specified by `DBtype`.
Folders must be given as vararg in the order given by `DBtype`.

e.g.:
```julia
flight = flightDB("i3", "data/inventory", "data/online_data")
```

# remarks
Any data can be attached to `FlightDB` with the keyword argument `remarks`.
"""
function loadFlightDB(DBtype::String, folder::Union{String, Vector{String}}...;
  altmin::Int=15_000, remarks=nothing)
  # Save time of database creation
  tc = Dates.now()
  # Find database types
  if occursin('i', DBtype)  i1 = findfirst(isequal('i'), DBtype)
  else  i1 = findfirst(isequal('1'), DBtype);  end
  if occursin('a', DBtype)  i2 = findfirst(isequal('a'), DBtype)
  else  i2 = findfirst(isequal('2'), DBtype);  end
  if occursin('o', DBtype)  i3 = findfirst(isequal('o'), DBtype)
  else  i3 = findfirst(isequal('3'), DBtype);  end

  # Load databases for each type
  if !isnothing(i1)
    ifiles = String[]
    ifiles = findFiles(ifiles, folder[i1], ".csv")
    inventory = try loadInventory(ifiles, altmin=altmin)
    catch
      @warn "Flight inventory couldn't be loaded."
      FlightData[]
    end
  else inventory = FlightData[];  end
  if !isnothing(i2)
    ifiles = String[]
    ifiles = findFiles(ifiles, folder[i2], ".csv")
    archive = try loadArchive(ifiles, altmin=altmin)
    catch
      @warn "FlightAware archive couldn't be loaded."
      FlightData[]
    end
  else archive = FlightData[];  end
  if !isnothing(i3)
    ifiles = String[]
    ifiles = findFiles(ifiles, folder[i3], ".txt", ".dat")
    onlineData = try loadOnlineData(ifiles, altmin=altmin)
    catch
      @warn "FlightAware online data couldn't be loaded."
      FlightData[]
    end
  else onlineData = FlightData[];  end

  println("\ndone loading data to properties\n- inventory\n- archive\n- onlineData\n", "")

  return FlightDB(inventory, archive, onlineData, tc, remarks)
end # function loadFlightDB


"""
    loadInventory(files::Vector{String}) -> inventory

From a list of `files`, return an `inventory` as `Vector{FlightData}` that can
be saved to the `inventory` field in `FlightDB`.
"""
function loadInventory(files::Vector{String}; altmin=15_000, filterCloudfree=true)

  # Initialise inventory file array and start MATLAB for PCHIP fitting
  inventory = FlightData[]

  # Loop over files
  for file in files

    # Load data
    flights = CSV.read(file, datarow=3, footerskip=2, ignoreemptylines=true,
      silencewarnings=true, threaded=true, dateformat="HH:MM:SS.sssm")

    # Calculate time from individual columns and add as DateTime to DataFrame
    flights.time = [ZonedDateTime(flights.SEGMENT_YEAR[i], flights.SEGMENT_MONTH[i],
      flights.SEGMENT_DAY[i], flights.SEGMENT_HOUR[i], flights.SEGMENT_MIN[i],
      flights.SEGMENT_SEC[i], tz.tz"UTC") for i = 1:length(flights.FLIGHT_ID)]

    # Initialise loop over file
    FID = flights.FLIGHT_ID[1]
    lat = Float64[]; lon = Float64[];
    alt = Float64[]; t = ZonedDateTime[]; speed = Float64[]

    # Loop over all data points
    @pm.showprogress 1 "load inventory from $(basename(splitext(file)[1]))..." for i = 1:length(flights.time)
      if flights.FLIGHT_ID[i] ≠ FID || i == length(flights.time)
        if length(t) ≤ 1  FID = flights.FLIGHT_ID[i]; continue  end
        lp = any(lon .> 0) ? maximum(lon[lon.≥0]) - minimum(lon[lon.≥0]) : 0
        ln = any(lon .< 0) ? maximum(lon[lon.<0]) - minimum(lon[lon.<0]) : 0
        useLON = maximum(lat) - minimum(lat) ≤ (lp + ln) * cosd(stats.mean(lat)) ? true : false
        useLON ? (x = lon; y = lat) : (x = lat; y = lon)
        x, y, alt, speed, t = remdup(x, y, alt, speed, t)
        flex = findFlex(x)
        push!(inventory, FlightData(t, lat, lon, alt, [missing for i = 1:length(t)],
          [missing for i = 1:length(t)], speed, FID, missing,
          missing, missing, flex, useLON, file))
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
  end #loop over files

  return inventory
end #function loadInventory


"""
    loadArchive(files::Vector{String}) -> archive

From a list of `files`, return an `archive` as `Vector{FlightData}` that can
be saved to the `archive` field in `FlightDB`.
"""
function loadArchive(files::Vector{String}; altmin::Int=15_000)
  # Initialise inventory file array
  archive = FlightData[]
  @pm.showprogress 1 "load archive..." for file in files
    # Load data
    flights = CSV.read(file, datarow=2, normalizenames=true, ignoreemptylines=true,
      silencewarnings=true, threaded=false, copycols=true, dateformat="m/d/y H:M:S",
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
      if flights.Flight_ID[i] ≠ FID || i == length(flights.Time_UTC_)
        if length(t) ≤ 1
          n = i
          FID = flights.Flight_ID[n]
          continue
        end
        lp = any(lon .> 0) ? maximum(lon[lon.≥0]) - minimum(lon[lon.≥0]) : 0
        ln = any(lon .< 0) ? maximum(lon[lon.<0]) - minimum(lon[lon.<0]) : 0
        useLON = maximum(lat) - minimum(lat) ≤ (lp + ln) * cosd(stats.mean(lat)) ? true : false
        flex = useLON ? findFlex(lon) : findFlex(lat)
        push!(archive, FlightData(t, lat, lon, alt, head, climb, speed,
        FID, flights.Ident[n], flights.Aircraft_Type[n],
        (orig=flights.Origin[n], dest=flights.Destination[n]), flex, useLON, file))

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
    end

    return archive
  end

  return archive
end #function loadArchive


"""
    loadOnlineData(files::Vector{String}) -> archive

From a list of `files`, return an `archive` as `Vector{FlightData}` that can
be saved to the `onlineData` field in `FlightDB`.
"""
function loadOnlineData(files::Vector{String}; altmin::Int=15_000)
  # Initialise inventory file array
  archive = FlightData[]
  @pm.showprogress 1 "load online data..." for (n, file) in enumerate(files)
    # Read flight data
    flight = CSV.read(file, ignoreemptylines=true, normalizenames=true, copycols=true,
      silencewarnings=true, threaded=false, types=Dict(:Latitude => Float64,
      :Longitude => Float64, :feet => String, :kts => Float64, :Course => String,
      :Rate => String))

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
    flightID, datestr, course = match(r"(.*?)_(.*?)_(.*)", filename).captures
    orig, dest = match(r"(.*)[-|_](.*)", course).captures
    date = Dates.Date(datestr, "d-u-y", locale="english")
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
        if flight[i,1][end-1:end] == "PM" && !(Dates.hour(t)==12 && Dates.minute(t)==0 &&
          Dates.second(t)==0)
          t += Dates.Hour(12)
        elseif flight[i,1][end-1:end] == "AM" && Dates.hour(t)==12 &&
          Dates.minute(t)==0 && Dates.second(t)==0
          t -= Dates.Hour(12)
        end
      end
      # Save data that needed tweaking of current time step
      push!(flighttime, ZonedDateTime(DateTime(date, t), timezone))
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
      (orig=orig, dest=dest), flex, useLON, file))
  end #loop over files

  return archive
end #function loadOnlineData
