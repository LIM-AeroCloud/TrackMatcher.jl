### Routines related to loading FlightTrack

"""
    load_volpe(
        paths::Vector{@NamedTuple{root::String,files::Vector{String}}},
        roots::ds.OrderedDict{String,UInt16},
        altmin::Real=5_000,
        Float::DataType=Float32
    ) -> StructArray{FlightData{Float}}

From a list of `files`, return a `StructArray{FlightData{Float}}` that can
be saved to the `inventory` field in `FlightSet`.

When the `StructArray{FlightData{Float}}` is constructed, data can be filtered by a minimum
altitude threshold in meters of the aircraft data (default: `altmin=5_000`).
Floating point numbers in `FlightData{Float}` are of the precision set by `Float`,
by default `Float32`.
"""
function load_volpe(
    paths::Vector{@NamedTuple{root::String,files::Vector{String}}},
    roots::ds.OrderedDict{String,UInt16},
    altmin::Real=5_000,
    Float::DataType=Float32
)::StructArray{FlightData{Float}}
    # Initialise loop over file
    inventory = StructArray{FlightData{Float}}(undef, 0)
    isempty(paths) && return inventory
    # Loop over files
    n = sum(length.(getfield.(paths, :files)))
    prog = pm.Progress(n, desc = "load VOLPE...")
    for path in paths, file in path.files
        # Load data file
        flightdata = CSV.read(joinpath(path.root, file), DataFrame, skipto=3, footerskip=2,
            ignoreemptyrows=true, silencewarnings=true, dateformat="HH:MM:SS.sssm",
            types = Dict("FLIGHT_ID" => Int, "LATITUDE" => Float, "LONGITUDE" => Float,
            "ALTITUDE" => Float, "SPEED" => Float, "SEGMENT_YEAR" => Int,
            "SEGMENT_MONTH" => Int, "SEGMENT_DAY" => Int, "SEGMENT_HOUR" => Int,
            "SEGMENT_MIN" => Int, "SEGMENT_SEC" => Int), select = [1:11;16])
        # Prepare and filter data
        df.transform!(flightdata,
            [:SEGMENT_YEAR, :SEGMENT_MONTH, :SEGMENT_DAY,
            :SEGMENT_HOUR, :SEGMENT_MIN, :SEGMENT_SEC] =>
            df.ByRow((y,m,d,h,min,s) -> DateTime(y,m,d,h,min,s)) => :time)
        df.select!(flightdata,
            :FLIGHT_ID => :FID,
            :time,
            :LATITUDE => :lat,
            :LONGITUDE => :lon,
            :ALTITUDE => :alt,
            :SPEED => :speed)
        flightdata.alt = ft2m.(flightdata.alt)
        flightdata.speed = knot2mps.(flightdata.speed)
        flightdata = flightdata[(ismissing.(flightdata.alt) .| (flightdata.alt .≥ altmin)), :]
        # Group by individual flights
        groups = df.groupby(flightdata, :FID)
        sizehint!(inventory, length(inventory) + length(groups))
        # Loop over flights
        for group in groups
            df.nrow(group) ≤ 1 && continue
            track = copy(group) # ℹ avoid problems with views when modifying data
            df.select!(track, df.Not(:FID))
            flex, useLON = preptrack!(track)
            isempty(flex) || push!(inventory, FlightData{Float}(track, group.FID[1],
                missing, missing, missing, flex, useLON, 0x01, roots[path.root], file))
        end #loop over flights
        # Monitor progress for progress bar
        pm.next!(prog, showvalues = [(:file,splitext(basename(file))[1])])
    end #loop over files
    pm.finish!(prog)

    return inventory
end #function loadVOLPE


"""
    load_flightaware(files::String...; Float::DataType=Float32, altmin::Real=5_000)
        -> StructArray{FlightData{Float}}

From a list of `files`, return a `Vector{FlightTrack}` that can
be saved to the `archive` field in `FlightSet`.

When the `Vector{FlightTrack{T}}` is constructed, data can be filtered by a minimum
altitude threshold in meters of the aircraft data (default: `altmin=5_000`).
Floating point numbers in `FlightTrack` are of the precision set by `Float`,
by default `Float32`.
"""
function load_flightaware(
    paths::Vector{@NamedTuple{root::String,files::Vector{String}}},
    roots::ds.OrderedDict{String,UInt16},
    altmin::Real=5_000,
    Float::DataType=Float32
)::StructArray{FlightData{Float}}
    # Initialise archive file array
    archive = StructArray{FlightData{Float}}(undef, 0)
    isempty(paths) && return archive
    # Loop over database files
    prog = pm.Progress(sum(length.(getfield.(paths, :files))), desc = "load archive...")
    for path in paths, file in path.files
        # Load data from csv file into standardised DataFrame
        filepath = joinpath(path.root, file)
        flights = read_archive(filepath, Float)
        # Unit conversions
        flights.alt = ft2m.(flights.alt)
        flights.speed = knot2mps.(flights.speed)
        flights.climb = ftpmin2mps.(flights.climb)

        # Filter data: keep if lat/lon not missing and alt is missing or alt >= altmin
        flights = flights[(.!ismissing.(flights.lat)) .& (.!ismissing.(flights.lon)) .&
            (ismissing.(flights.alt) .| (flights.alt .≥ altmin)), :]

        # Group by flight ID
        groups = df.groupby(flights, :dbID)
        sizehint!(archive, length(archive) + length(groups))

        # Process each flight
        for group in groups
            df.nrow(group) ≤ 1 && continue

            # Copy to DataFrame and select track columns
            track = df.select(copy(group), :time, :lat, :lon, :alt, :speed, :climb, :heading)

            # Determine predominant flight direction, inflection points, and remove duplicates
            flex, useLON = preptrack!(track)

            # Save the FlightTrack in the archive vector
            isempty(flex) || push!(archive, FlightData{Float}(track, group.dbID[1],
                group.flightID[1], group.type[1],
                (orig=group.orig[1], dest=group.dest[1]), flex, useLON,
                0x02, roots[path.root], file))
        end #loop over flights

        # Monitor progress for progress bar
        pm.next!(prog, showvalues = [(:file,splitext(basename(file))[1])])
    end #loop over files
    pm.finish!(prog)

    return archive
end #function load_flightaware


"""
    load_webdata(
        paths::Vector{@NamedTuple{root::String,files::Vector{String}}},
        roots::ds.OrderedDict{String,UInt16},
        altmin::Real=5_000,
        Float::DataType=Float32,
        delim::Union{Nothing,Char,String}=nothing
    ) -> StructArray{FlightData{Float}}

From a list of `files`, return a `Vector{FlightTrack}` that can
be saved to the `onlineData` field in `FlightSet`.

The `delim`iter of the data in the input file can be specified by a string or character.
Default is `nothing`, which means auto-detection of the delimiter is used.

When the `Vector{FlightTrack{T}}` is constructed, data can be filtered by a minimum
altitude threshold in meters of the aircraft data (default: `altmin=5_000`).
Floating point numbers in `FlightTrack` are of the precision set by `Float`,
by default `Float32`.
"""
function load_webdata(
    paths::Vector{@NamedTuple{root::String,files::Vector{String}}},
    roots::ds.OrderedDict{String,UInt16},
    altmin::Real=5_000,
    Float::DataType=Float32,
    delim::Union{Nothing,Char,String}=nothing
)::StructArray{FlightData{Float}}
    # Initialise inventory file array
    archive = StructArray{FlightData{Float}}(undef, 0)
    isempty(paths) && return archive
    # Loop over files with online data
    prog = pm.Progress(sum(length.(paths.files)), desc = "load online data...")
    for path in paths, file in path.files
        # Read flight data
        flight = CSV.read(file, DataFrame, delim=delim, ignoreemptylines=true, normalizenames=true,
        silencewarnings=true, ntasks=1, copycols=false, types=Dict(:Latitude => Float,
        :Longitude => Float, :feet => String, :kts => Float, :Course => String,
        :Rate => String), drop = ["mph", "Reporting_Facility"], stringtype=String)
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
        date, timezone, flightID, orig, dest = get_date_time_route(filename, tzone)

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
        isempty(flex) || push!(archive, FlightData{Float}(flight, replace(filename, "_" => "/"),
        flightID, missing, (orig=orig, dest=dest), flex, useLON, 0x03, roots[path.root], file))
        # Monitor progress for progress bar
        pm.next!(prog)
    end #loop over files
    pm.finish!(prog)

    return archive
end #function loadWD


"""
    read_archive(file::AbstractString, Float::DataType=Float32) -> DataFrame

Read FlightAware archived data from a csv `file` and return content as DataFrame.

The routine works for several FlightAware archive versions. Floating point numbers
are read with single precision or as defined by kwarg `Float`.
"""
function read_archive(file::AbstractString, Float::DataType=Float32)::DataFrame
    # Read header and check version
    header = readline(file)
    old = contains(header, "(kts)")
    # Load correct data version and return as DataFrame
    if old
        CSV.read(file, DataFrame, skipto=2, normalizenames=true, ignoreemptyrows=true,
            silencewarnings=true, dateformat="m/d/y H:M:S",
            header = ["dbID", "flightID", "type", "orig", "dest", "time", "lat", "lon",
            "speed", "alt", "climb", "heading", "direction", "estimated"],
            types = Dict(:lat => Float, :lon => Float, :alt => Float, :speed => Float, :climb => Float),
            drop = ["direction", "estimated"], stringtype=String)
    else
        CSV.read(file, DataFrame, skipto=2, normalizenames=true, ignoreemptyrows=true,
            silencewarnings=true, dateformat="m/d/y H:M:S",
            header = ["dbID", "flightID", "type", "orig", "dest", "time", "lat", "lon",
            "speed", "alt", "climb", "heading", "direction", "fac_name", "fac_descr", "estimated"],
            types = Dict(:lat => Float, :lon => Float, :alt => Float, :speed => Float, :climb => Float),
            drop = ["direction", "fac_name", "fac_descr", "estimated"], stringtype=String)
    end
end #function read_archive


"""
    get_date_time_route(filename::String, tzone::String)

From the `filename` and a custom time zone string (`tzone`), extract and return
the starting date, the standardised time zone, the flight ID, origin, and destination.
"""
function get_date_time_route(filename::String, tzone::String)

    # Time is the first column and has to be addressed as flight[!,1] in the code
    # due to different column names, in which the timezone is included
    timezone = zonedict[tzone]
    # Retrieve date and metadata from filename
    flightID, datestr, course = try match(r"(.*?)_(.*?)_(.*)", filename).captures
    catch
        println()
        println()
        @warn "Flight ID, date, and course not found. Data skipped." file
        return missing, missing, missing, missing, missing
    end
    orig, dest = match(r"(.*)[-|_](.*)", course).captures
    date = try Dates.Date(datestr, "d-u-y", locale="english")
    catch
        println()
        println()
        @warn "Unable to parse date. Data skipped." file
        return missing, missing, missing, missing, missing
    end

    return date, timezone, flightID, orig, dest
end #function get_date_time_route
