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
function loadFlightDB(DBtype::String, folder::Union{String, Vector{String}}...; remarks=nothing)
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
    inventory = try loadInventory(ifiles)
    catch
      @warn "Flight inventory couldn't be loaded."
      FlightData[]
    end
  else inventory = FlightData[];  end
  if !isnothing(i2)
    ifiles = String[]
    ifiles = findFiles(ifiles, folder[i2], ".csv")
    archive = try loadArchive(ifiles)
    catch
      @warn "FlightAware archive couldn't be loaded."
      FlightData[]
    end
  else archive = FlightData[];  end
  if !isnothing(i3)
    ifiles = String[]
    ifiles = findFiles(ifiles, folder[i3], ".txt", ".dat")
    onlineData = try loadOnlineData(ifiles)
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
function loadInventory(files::Vector{String}; filterAlt=15000, filterCloudfree=true)

  # Initialise inventory file array and start MATLAB for PCHIP fitting
  inventory = FlightData[]
  ms = mat.MSession() # rather hand over as function argument?

  # Loop over files
  for file in files

    # Load data
    flights = CSV.read(file, datarow=3, footerskip=2,
      ignoreemptylines=true, silencewarnings=true)

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
        pchip = PCHIP(x, y, flex, FID, useLON, ms)
        push!(inventory, FlightData(t, lat, lon, alt, [missing for i = 1:length(t)],
          [missing for i = 1:length(t)], speed, FID, missing,
          missing, missing, pchip, file))
        lat = Float64[]; lon = Float64[];
        alt = Float64[]; t = ZonedDateTime[]; speed = Float64[]
        FID = flights.FLIGHT_ID[i]
        # segtime = DateTime[]; segdist = Float64[]
      end
      # Filter altitude threshold
      if flights.ALTITUDE[i] ≥ filterAlt
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
function loadArchive(files::Vector{String})
  # Initialise inventory file array
  archive = FlightData[]
  @pm.showprogress 1 "load archive..." for file in files
    # Load data
    flights = CSV.read(file, datarow=2, header=[:id, :ident, :orig, :dest,
      :aircraft, :time, :lat, :lon, :speed, :alt, :climb, :heading, :direction,
      :facility, :description, :est])
    flights.speed = convert.(Union{Missing,Float64}, flights.speed)
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
    altitude = Vector{Union{Missing,Float64}}(undef, length(flight.time))
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
        alt=parse(Float64, flight.alt[i][findall([isdigit.(i) for i ∈ flight.alt[i]])])
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
      flight.heading, flight.climb, convert.(Union{Missing,Float64},flight.speed),
      n, flightID, missing, (orig=orig, dest=dest), file)
  end #loop over files

  return archive
end #function loadOnlineData
