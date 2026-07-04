### Routines related to loading FlightTrack

"""
    load_volpe(
        paths::Vector{@NamedTuple{root::String,files::Vector{String}}},
        pathdict::ds.OrderedDict{String,ds.OrderedDict},
        altmin::Real=5_000,
        Float::DataType=Float32
    ) -> StructArray{FlightData{Float}}

From a list of `paths`, return a `StructArray{FlightData{Float}}` that can
be saved to the `volpe` field in `FlightSet`. The `pathdict` connects file names with integer
indices used in the metadata of `FlightData` and `FlightSet` to save memory.

When the `StructArray{FlightData{Float}}` is constructed, data can be filtered by a minimum
altitude threshold in meters of the aircraft data (default: `altmin=5_000`).
Floating point numbers in `FlightData{Float}` are of the precision set by `Float`,
by default `Float32`.
"""
function load_volpe(
    paths::Vector{@NamedTuple{root::String,files::Vector{String}}},
    pathdict::ds.OrderedDict{String,ds.OrderedDict},
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
            "SEGMENT_MIN" => Int, "SEGMENT_SEC" => Int), select = [1;3:11;16])
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
            flex, use_lon = preptrack!(track)
            isempty(flex) || push!(inventory, FlightData{Float}(track, group.FID[1],
                missing, missing, missing, flex, use_lon, "V", pathdict["roots"][path.root],
                pathdict["files"][file]))
        end #loop over flights
        # Monitor progress for progress bar
        pm.next!(prog, showvalues = [(:file,splitext(basename(file))[1])])
    end #loop over files
    pm.finish!(prog)

    return inventory
end #function loadVOLPE


"""
    load_flightaware(
        paths::Vector{@NamedTuple{root::String,files::Vector{String}}},
        pathdict::ds.OrderedDict{String,ds.OrderedDict},
        altmin::Real=5_000,
        Float::DataType=Float32
    ) -> StructArray{FlightData{Float}}

From a list of `paths`, return a `StructArray{FlightData{Float}}` that can
be saved to the `flightaware` field in `FlightSet`. The `pathdict` connects file names with
integer  indices used in the metadata of `FlightData` and `FlightSet` to save memory.

When the `StructArray{FlightData{Float}}` is constructed, data can be filtered by a minimum
altitude threshold in meters of the aircraft data (default: `altmin=5_000`).
Floating point numbers in `FlightTrack` are of the precision set by `Float`,
by default `Float32`.
"""
function load_flightaware(
    paths::Vector{@NamedTuple{root::String,files::Vector{String}}},
    pathdict::ds.OrderedDict{String,ds.OrderedDict},
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
        groups = df.groupby(flights, :id)
        sizehint!(archive, length(archive) + length(groups))
        # Process each flight
        for group in groups
            df.nrow(group) ≤ 1 && continue
            # Copy to DataFrame and select track columns
            track = df.select(copy(group), :time, :lat, :lon, :alt, :speed, :climb, :heading)
            # Determine predominant flight direction, inflection points, and remove duplicates
            flex, use_lon = preptrack!(track)
            # Save the FlightTrack in the archive vector
            isempty(flex) || push!(archive, FlightData{Float}(track, group.id[1],
                group.flight_num[1], group.type[1],
                (orig=group.orig[1], dest=group.dest[1]), flex, use_lon,
                "FA", pathdict["roots"][path.root], pathdict["files"][file]))
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
        pathdict::ds.OrderedDict{String,ds.OrderedDict},
        altmin::Real=5_000,
        Float::DataType=Float32,
        delim::Union{Nothing,Char,String}=nothing
    ) -> StructArray{FlightData{Float}}

From a list of `paths`, return a `StructArray{FlightData{Float}}` that can
be saved to the `webdata` field in `FlightSet`. The `pathdict` connects file names with
integer  indices used in the metadata of `FlightData` and `FlightSet` to save memory.

The `delim`iter of the data in the input file can be specified by a string or character.
Default is `nothing`, which means auto-detection of the delimiter is used.

When the `Vector{FlightTrack{T}}` is constructed, data can be filtered by a minimum
altitude threshold in meters of the aircraft data (default: `altmin=5_000`).
Floating point numbers in `FlightTrack` are of the precision set by `Float`,
by default `Float32`.
"""
function load_webdata(
    paths::Vector{@NamedTuple{root::String,files::Vector{String}}},
    pathdict::ds.OrderedDict{String,ds.OrderedDict},
    altmin::Real=5_000,
    Float::DataType=Float32,
    delim::Union{Nothing,Char,String}=nothing
)::StructArray{FlightData{Float}}
    # Initialise inventory file array
    archive = StructArray{FlightData{Float}}(undef, 0)
    isempty(paths) && return archive
    nlength = sum(length.(getfield.(paths, :files)))
    sizehint!(archive, nlength)
    # Loop over files with online data
    prog = pm.Progress(nlength, desc = "load web data...")
    for path in paths, file in path.files
        # Read flight data
        flight = try CSV.read(joinpath(path.root, file), DataFrame, delim=delim, ignoreemptyrows=true,
            normalizenames=true, silencewarnings=true, ntasks=1, copycols=false,
            types=Dict(:Latitude => Float, :Longitude => Float, :feet => String, :kts => Float,
            :Course => String, :Rate => String), drop = ["mph", "Reporting_Facility"], stringtype=String)
        catch err
            println()
            println()
            @error "Error reading file. Try to specify column delimiter. Data skipped." file exception=(err, catch_backtrace())
            continue
        end
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
        date, timezone, flight_num, orig, dest = get_date_time_route(filename, tzone)

        # Set to 2 days prior to allow corrections for timezone diffences in the next step
        date += Dates.Day(2)
        ### Convert times to datetime and extract heading and climbing rate
        # Initialise vectors
        lat = Float[]
        lon = Float[]
        flighttime = ZonedDateTime[]
        altitude = Union{Missing,Float}[]
        heading = Union{Missing,Int}[]
        climbingrate = Union{Missing,Float}[]
        speed = Union{Missing,Float}[]
        # Loop over times (backwards to preserve day matching logic)
        for i = length(flight[!,1]):-1:1
            alt = try ft2m(parse(Float, join([n for n in flight.alt[i]
                if isnumeric(n) || n == '.'])))
            catch; missing; end
            climb = try ftpmin2mps(parse(Float, join([n for n in flight.climb[i]
                if isnumeric(n) || n == '.' || n == '-'])))
            catch; missing; end
            head = try parse(Int, join([n for n in flight.heading[i]
                if isnumeric(n) || n == '.']))
            catch; missing; end
            if length(flight.time[i]) ≠ 15 || ismissing(flight.lat[i]) ||
                ismissing(flight.lon[i]) || (!ismissing(alt) && alt < altmin)
                continue
            end
            # Derive date from day of week and filename
            while Dates.dayabbr(date) ≠ flight.time[i][1:3]
                date -= Dates.Day(1)
            end
            # Derive time from time string
            t = Time(flight.time[i][5:end], "I:M:S p")
            # Save data (reverse order for now)
            push!(flighttime, ZonedDateTime(DateTime(date, t), timezone))
            push!(altitude, ismissing(alt) ? missing : Float(alt))
            push!(climbingrate, ismissing(climb) ? missing : Float(climb))
            push!(heading, head)
            push!(lat, Float(flight.lat[i]))
            push!(lon, Float(flight.lon[i]))
            push!(speed, ismissing(flight.speed[i]) ? missing : Float(flight.speed[i]))
        end #loop of flight

        # Skip data with all data points below the altitude threshold or missing
        isempty(altitude) && continue
        # Restore original order
        reverse!(flighttime)
        reverse!(altitude)
        reverse!(climbingrate)
        reverse!(heading)
        reverse!(lat)
        reverse!(lon)
        reverse!(speed)
        # Build track DataFrame
        track = DataFrame(time = flighttime, lat = lat, lon = lon, alt = altitude,
            speed = speed, climb = climbingrate, heading = heading)

        # Determine predominant flight direction, inflection points, and remove duplicate entries
        flex, use_lon = preptrack!(track)

        # Save data as FlightTrack
        isempty(flex) || push!(archive, FlightData{Float}(track, replace(filename, "_" => "/"),
            flight_num, missing, (orig=orig, dest=dest), flex, use_lon, "W",
            pathdict["roots"][path.root], pathdict["files"][file]))
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
            header = ["id", "flight_num", "type", "orig", "dest", "time", "lat", "lon",
            "speed", "alt", "climb", "heading", "direction", "estimated"],
            types = Dict(:lat => Float, :lon => Float, :alt => Float, :speed => Float, :climb => Float),
            drop = ["direction", "estimated"], stringtype=String)
    else
        CSV.read(file, DataFrame, skipto=2, normalizenames=true, ignoreemptyrows=true,
            silencewarnings=true, dateformat="m/d/y H:M:S",
            header = ["id", "flight_num", "type", "orig", "dest", "time", "lat", "lon",
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
    flight_num, datestr, course = try match(r"(.*?)_(.*?)_(.*)", filename).captures
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

    return date, timezone, flight_num, orig, dest
end #function get_date_time_route
