## Metadata structs

"""
# struct FlightMetadata{T<:AbstractFloat} <: FlightTrack{T}

Immutable struct holding metadata for an individual `FlightTrack` of the `FlightSet` with fields

- `dbID::Union{Int,AbstractString}`
- `flightID::Union{Missing,AbstractString}`
- `route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}}`
- `aircraft::Union{Missing,AbstractString}`
- `date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}`
- `area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,AbstractFloat}}`
- `flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,AbstractFloat,AbstractFloat}}}}`
- `useLON::Bool`
- `source::UInt8`
- `root::UInt16`
- `file::UInt16`

By default, `Float32` is used for `T`.

See also [`PrimaryMetadata`](@ref), [`FlightData`](@ref), [`FlightTrack`](@ref), [`FlightSet`](@ref), and [`PrimarySet`](@ref).

## dbID
Database ID – integer counter for `volpe`,
String with information about `FlightID`, `route`, and/or scheduled arrival for
FlightAware data.

## FlightID and aircraft
Strings with aircraft identification and type.

## route
`NamedTuple` with fields for `orig`in and `dest`ination holding ICAO airport codes.

## area
`NamedTuple` with fields for latitude and Longitude range. For the longitude range,
it is distinguished between positive and negative ranges to avoid problems with
flights passing the date line.

Fields:
- `latmin`
- `latmax`
- `elonmin`
- `elonmax`
- `wlonmin`
- `wlonmax`

## date
`NamedTuple` with fields `start` and `stop` for start and end time of the current
flight.

## flex
`Tuple` of `NamedTuple`s with entries
- `range` (`UnitRange`): flight segment between inflection points of the current flight track
- `min` (`AbstractFloat`): minimum x value in the flight segment
- `max` (`AbstractFloat`): maximum x value in the flight segment

## useLON
Flag (`Bool`) whether to use longitude as x data for track interpolation.

## source
String describing the database source of the current flight:
- `"VOLPE"`
- `"FlightAware"`
- `"flightaware.com"`

## file
String holding the file name and path relative to the root directory.

## root
String holding the root directory path.


## Instantiation

`FlightMetadata` is constructed automatically, when `FlightData` is instantiated using
a modified constructor and `dbID`, `flightID`, `route`, `aircraft` type, the `date` and `lat`/`lon`
vectors, `useLON`, `flex`, `source`, `root`, and `file`.
Fields `area` and `date` are calculated from `lat`/`lon`, and `date` vectors.

    function FlightMetadata{T}(
        dbID::Union{Int,AbstractString},
        flightID::Union{Missing,AbstractString},
        route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
        aircraft::Union{Missing,AbstractString},
        date::Vector{DateTime},
        lat::Vector{<:Union{Missing,T}},
        lon::Vector{<:Union{Missing,T}},
        useLON::Bool,
        flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}},
        source::UInt8,
        root::UInt16,
        file::UInt16
    ) where T<:AbstractFloat -> struct FlightMetadata

Or construct `FlightMetadata` by directly handing over every field:

    function FlightMetadata{T}(
        dbID::Union{Int,AbstractString},
        flightID::Union{Missing,AbstractString},
        route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
        aircraft::Union{Missing,AbstractString},
        date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
        area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,T}},
        flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}},
        useLON::Bool,
        source::UInt8,
        root::UInt16,
        file::UInt16
    ) where T<:AbstractFloat -> struct FlightMetadata
"""
struct FlightMetadata{T<:AbstractFloat} <: FlightTrack{T}
    dbID::Union{Int,AbstractString}
    flightID::Union{Missing,AbstractString}
    route::Union{Missing,NamedTuple{(:orig,:dest),Tuple{AbstractString,AbstractString}}}
    aircraft::Union{Missing,AbstractString}
    date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
    area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,T}}
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}}
    useLON::Bool
    source::UInt8
    root::UInt16
    file::UInt16
end

#* Main constructor with limited tests and limited automated data construction
function FlightMetadata{T}(
    dbID::Union{Int,AbstractString},
    flightID::Union{Missing,AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
    aircraft::Union{Missing,AbstractString},
    date::Vector{DateTime},
    lat::Vector{<:Union{Missing,T}},
    lon::Vector{<:Union{Missing,T}},
    useLON::Bool,
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}},
    source::UInt8,
    root::UInt16,
    file::UInt16
) where T<:AbstractFloat
    elonmin, elonmax = lonextrema(lon, ≥)
    wlonmin, wlonmax = lonextrema(lon, <)
    area = (latmin=minimum(lat), latmax=maximum(lat), elonmin, elonmax, wlonmin, wlonmax)
    FlightMetadata{T}(dbID, flightID, route, aircraft, (start=date[1], stop=date[end]), area,
        flex, useLON, source, root, file)
end #constructor 2 FlightMetadata

#* Constructor for empty FlightMetadata
FlightMetadata{T}() where T<:AbstractFloat = FlightMetadata{T}(
    "", missing, missing, missing, (start=Dates.now(), stop=Dates.now()),
    (latmin=T(NaN), latmax=T(NaN), elonmin=T(NaN), elonmax=T(NaN), wlonmin=T(NaN), wlonmax=T(NaN)),
    ((range=0:0, min=T(NaN), max=T(NaN)),), false, 0x00, 0x0000, 0x0000
)

#* Constructor for default Float32 FlightMetadata
FlightMetadata(args...) = FlightMetadata{Float32}(args...)

#* Constructor for floating point type promotion
FlightMetadata{T}(meta::FlightMetadata) where T<:AbstractFloat = FlightMetadata{T}(
    meta.dbID, meta.flightID, meta.route, meta.aircraft, meta.date,
    (latmin = T(meta.area.latmin), latmax = T(meta.area.latmax), elonmin = T(meta.area.elonmin),
    elonmax = T(meta.area.elonmax), wlonmin = T(meta.area.wlonmin), wlonmax = T(meta.area.wlonmax)),
    Tuple([(range = m.range, min = T.(m.min), max = T.(m.max)) for m in meta.flex]),
    meta.useLON, meta.source, meta.root, meta.file
)


"""
# struct PrimaryMetadata{T<:AbstractFloat} <: PrimarySet{T}

Immutable struct with additional information of databases:

- `altmin::Real`: Minimum altitude threshold for which flight data is considered
- `date`: NamedTuple with entries `start`/`stop` giving the time range of the database
- `sources`: OrderedDict with a lookup table for database type
- `pathlookup`: OrderedDict with a lookup table for root directories of the given files and
  file names for the  flight data
- `created`: time of creation of database
- `loadtime`: time it took to read data files and load it to the struct
- `attachments`: any additional data or comments that can be attached to the database

See also [`FlightMetadata`](@ref), [`FlightData`](@ref), [`FlightTrack`](@ref), [`FlightSet`](@ref), and [`PrimarySet`](@ref).
"""
struct PrimaryMetadata{T<:AbstractFloat} <: PrimarySet{T}
    altmin::T
    date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
    sources::ds.OrderedDict{UInt8,String}
    pathlookup::ds.OrderedDict{String,ds.OrderedDict}
    created::Union{DateTime,ZonedDateTime}
    loadtime::Dates.CompoundPeriod
    attachments
end #struct PrimaryMetadata

#* Constructor for empty PrimaryMetadata
PrimaryMetadata{T}() where T<:AbstractFloat = PrimaryMetadata{T}(NaN,
    (start=Dates.now(), stop=Dates.now()), ds.OrderedDict{UInt8,String}(0x00 => "undefined"),
        ds.OrderedDict("roots" => ds.OrderedDict{UInt16,String}(0x0000 => "/"),
        "files" => ds.OrderedDict{UInt16,String}(0x0000 => "")),
        Dates.now(), Dates.CompoundPeriod(), nothing)

#* Constructor for default Float32 PrimaryMetadata
PrimaryMetadata(args...) = PrimaryMetadata{Float32}(args...)

#* Constructor for floating point type promotion
PrimaryMetadata{T}(meta::PrimaryMetadata) where T<:AbstractFloat = PrimaryMetadata{T}(T(meta.altmin),
    meta.date, meta.sources, meta.pathlookup, meta.created, meta.loadtime, meta.attachments)


## Struct for single flight tracks

"""
# struct FlightData{T<:AbstractFloat} <: FlightTrack{T}

Immutable struct holding data for an individual flight track of the `FlightSet` with fields

- `time::Vector{DateTime}` for UTC datetimes of the flight track points
- `lat::Vector{T}` for latitudes of the flight track points
- `lon::Vector{T}` for longitudes of the flight track points
- `alt::Vector{<:Union{Missing,T}}` for altitudes of the flight track points in meters
- `heading::Vector{<:Union{Missing,Int}}` for headings of the flight track points in degrees
- `climb::Vector{<:Union{Missing,T}}` for climb rates of the flight track points in meters per second
- `speed::Vector{<:Union{Missing,T}}` for speeds at each flight track point in meters per second
- `metadata::FlightMetadata{T}` for metadata associated with the flight track

By default, `Float32` is used for `T`.

See also [`FlightTrack`](@ref), [`FlightMetadata`](@ref), [`FlightSet`](@ref), and [`PrimarySet`](@ref).

## Instantiation

`FlightData` can be instantiated by directly handing over all data vectors and metadata Fields
or using a modified constructor that takes a `DataFrame` with the flight track data read from
a csv file and the metadata fields:

    function FlightData{T}(
        track::DataFrame,
        dbID::Union{Int,AbstractString},
        flightID::Union{Missing,AbstractString},
        aircraft::Union{Missing,AbstractString},
        route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{<:AbstractString,<:AbstractString}}},
        flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}},
        useLON::Bool,
        source::UInt8,
        root::UInt16,
        file::UInt16
    ) where T<:AbstractFloat -> struct FlightData
"""
struct FlightData{T<:AbstractFloat} <: FlightTrack{T}
    time::Vector{DateTime}
    lat::Vector{T}
    lon::Vector{T}
    alt::Vector{<:Union{Missing,T}}
    heading::Vector{<:Union{Missing,Int}}
    climb::Vector{<:Union{Missing,T}}
    speed::Vector{<:Union{Missing,T}}
    metadata::FlightMetadata{T}

    #* Internal constructor with data checks
    function FlightData{T}(time::Vector{DateTime}, lat::Vector{T}, lon::Vector{T},
        alt::Vector{<:Union{Missing,T}}, heading::Vector{<:Union{Missing,Int}},
        climb::Vector{<:Union{Missing,T}}, speed::Vector{<:Union{Missing,T}},
        metadata::FlightMetadata{T}) where T
        length(time) == length(lat) == length(lon) == length(alt) == length(heading) ==
            length(climb) == length(speed) ||
            throw(DimensionMismatch("all data vectors must have the same length in FlightData"))
        checklimits(time, DateTime(2000), Dates.now(), "time")
        checklimits(lat, -90, 90, "latitude")
        checklimits(lon, -180, 180, "longitude")
        checklimits(alt, 0, 40000, "altitude")
        checklimits(heading, 0, 360, "heading")
        checklimits(speed, 0, 2500, "speed")
        new{T}(time, lat, lon, alt, heading, climb, speed, metadata)
    end
end #struct FlightData


#* Main constructor parsing a DataFrame from file input and ensuring UTC time
function FlightData{T}(
    track::DataFrame,
    dbID::Union{Int,<:AbstractString},
    flightID::Union{Missing,<:AbstractString},
    aircraft::Union{Missing,<:AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{<:AbstractString,<:AbstractString}}},
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}},
    useLON::Bool,
    source::UInt8,
    root::UInt16,
    file::UInt16
) where T<:AbstractFloat
    # Check dataframe columns of flight data; fill missing columns with missing values
    t = track.time
    t = t isa Vector{ZonedDateTime} ? [zt.utc_datetime for zt in t] : t
    lat = track.lat
    lon = track.lon
    alt = track.alt
    heading = hasproperty(track, :heading) ? track.heading : [missing for _ in t]
    climb = hasproperty(track, :climb) ? track.climb : [missing for _ in t]
    speed = hasproperty(track, :speed) ? track.speed : [missing for _ in t]
    metadata = FlightMetadata{T}(dbID, flightID, route, aircraft, t, lat, lon, useLON, flex,
        source, root, file)

    # Instatiate new FlightTrack
    FlightData{T}(t, lat, lon, alt, heading, climb, speed, metadata)
end #constructor 2 FlightData

#* Constructor for empty FlightData
FlightData{T}() where T<:AbstractFloat = FlightData{T}(DateTime[], T[], T[],
    Union{Missing,T}[], Union{Missing,Int}[], Union{Missing,T}[], Union{Missing,T}[],
    FlightMetadata{T}())

#* Constructor for type promotion of FlightData
FlightData{T}(flight::FlightData) where T<:AbstractFloat =
  FlightData{T}(flight.time, T.(flight.lat), T.(flight.lon),
    [ismissing(x) ? missing : T(x) for x in flight.alt], flight.heading,
    [ismissing(x) ? missing : T(x) for x in flight.climb],
    [ismissing(x) ? missing : T(x) for x in flight.speed],
    FlightMetadata{T}(flight.metadata))

#* Constructor for default Float32 FlightData
FlightData(args...; kwargs...) = FlightData{Float32}(args...; kwargs...)


"""
    FlightTrack(args...; kwargs...)
    FlightTrack{T}(args...; kwargs...) where T<:AbstractFloat

Alias constructors for `FlightData` with default `Float32` precision and type promotion.

See also [`FlightData`](@ref), [`FlightMetadata`](@ref), [`FlightSet`](@ref), and [`PrimarySet`](@ref).
"""
function FlightTrack end

#* Alias constructor for FlightData with default Float32 precision
FlightTrack(args...; kwargs...) = FlightData{Float32}(args...; kwargs...)

#* Alias constructor for FlightData with type promotion
FlightTrack{T}(args...; kwargs...) where T<:AbstractFloat = FlightData{T}(args...; kwargs...)


## Struct for sets of flight tracks

"""
# struct FlightSet{T<:AbstractFloat} <: PrimarySet{T}

Immutable struct holding data for a set of flight tracks of the `FlightSet` with fields

- `volpe::StructArray{FlightData{T}}` for flight track data from the VOLPE AEDT inventory
- `flightaware::StructArray{FlightData{T}}` for flight track data from the commercially available database by FlightAware
- `webdata::StructArray{FlightData{T}}` for flight track data from free online data by FlightAware
- `metadata::PrimaryMetadata{T}` for metadata associated with the flight track set

By default, `Float32` is used for `T`.

`FlightSet` is used to store primary data and as source for TrackMatcher for intersection
calculations.

See also [`PrimarySet`](@ref), [`FlightData`](@ref), [`FlightTrack`](@ref), and [`PrimaryMetadata`](@ref).

## Instantiation

`FlightSet` can be instantiated by directly handing over the `StructArray`s for the different
inventories and the `PrimaryMetadata` or using a main constructor that takes folder paths for
the different inventories as `AbstractString` or `Vector{AbstractString}` and other parameters
to load the data and construct the `StructArray`s and `PrimaryMetadata` automatically.

    function FlightSet{T}(;
        volpe::Union{String,Vector{String}}=String[],
        flightaware::Union{String,Vector{String}}=String[],
        webdata::Union{String,Vector{String}}=String[],
        altmin::Real=5000,
        delim::Union{Nothing,Char,String}=nothing,
        attachments=nothing
    ) where T<:AbstractFloat -> struct FlightSet

### Keyword arguments for main constructor

- `volpe`: folder path(s) for the VOLPE AEDT inventory data as `AbstractString` or `Vector{AbstractString}`
- `flightaware`: folder path(s) for the commercially available database by FlightAware as `AbstractString` or `Vector{AbstractString}`
- `webdata`: folder path(s) for the free online data by FlightAware as `AbstractString` or `Vector{AbstractString}`
- `altmin`: minimum altitude threshold for which flight data is considered in meters (default: `5000`)
- `delim`: delimiter for the free online data by FlightAware (default: `nothing`, i.e. auto-detection)
- `attachments`: any additional data or comments that can be attached to the database (default: `nothing`)
"""
struct FlightSet{T<:AbstractFloat} <: PrimarySet{T}
    volpe::StructArray{FlightData{T}}
    flightaware::StructArray{FlightData{T}}
    webdata::StructArray{FlightData{T}}
    metadata::PrimaryMetadata{T}

    #* Internal constructor with data checks
    function FlightSet{T}(volpe::StructArray{FlightData{T}}, flightaware::StructArray{FlightData{T}},
        webdata::StructArray{FlightData{T}}, metadata::PrimaryMetadata{T}) where T<:AbstractFloat

        # Check for correct dataset type in each vector and for correct floating point precision
        deleteat!(volpe, findall(data -> data.metadata.source ≠ 0x01, volpe))
        deleteat!(flightaware, findall(data -> data.metadata.source ≠ 0x02, flightaware))
        deleteat!(webdata, findall(data -> data.metadata.source ≠ 0x03, webdata))

        # Instantiate new struct
        new{T}(volpe, flightaware, webdata, metadata)
    end #constructor 1 FlightSet
end #struct FlightSet

#* Main constructor
function FlightSet{T}(;
    volpe::Union{String,Vector{String}}=String[],
    flightaware::Union{String,Vector{String}}=String[],
    webdata::Union{String,Vector{String}}=String[],
    altmin::Real=5000,
    delim::Union{Nothing,Char,String}=nothing,
    attachments=nothing
) where T<:AbstractFloat

    # Return empty FlightSet, if no folders are passed to constructor
    all(isempty.([volpe, flightaware, webdata])) && return FlightSet{T}()
    # Save time of database creation
    tstart = Dates.now()

    ## Load databases for each type
    pathdict = ds.OrderedDict{String,ds.OrderedDict}(
        "roots" => ds.OrderedDict{String,UInt16}(),
        "files" => ds.OrderedDict{String,UInt16}()
    )
    # VOLPE AEDT dataset
    files = findfiles!(pathdict, volpe, ".csv")
    volpe = load_volpe(files, pathdict, altmin, T)
    # FlightAware commercial dataset
    files = findfiles!(pathdict, flightaware, ".csv")
    flightaware = load_flightaware(files, pathdict, altmin, T)
    # FlightAware web content
    files = findfiles!(pathdict, webdata, [".tsv", ".txt", ".dat"])
    webdata = load_webdata(files, pathdict, altmin, T, delim)
    # Determine overall time range of database from all flight tracks
    tmin, tmax = if isempty([volpe; flightaware; webdata])
        tstart, tstart
    else
        minimum([[f.metadata.date.start for f in volpe];
            [f.metadata.date.start for f in flightaware];
            [f.metadata.date.start for f in webdata]]),
        maximum([[f.metadata.date.stop for f in volpe];
            [f.metadata.date.stop for f in flightaware];
            [f.metadata.date.stop for f in webdata]])
    end

    ## Finalise and construct flight set
    # wrap up and calculate load time
    tend = Dates.now()
    tc = tz.ZonedDateTime(tend, tz.localzone())
    loadtime = Dates.canonicalize(Dates.CompoundPeriod(tend - tstart))

    @info string("FlightSet loaded in ",
        "$(join(loadtime.periods[1:min(2,length(loadtime.periods))], ", ")) to",
        "\n▪ volpe ($(length(volpe)) entries)\n▪ flightaware ($(length(flightaware)) entries)\n",
        "▪ webdata ($(length(webdata)) entries)\n▪ metadata")

    # Create lookup dictionaries for source and root IDs
    sources = ds.OrderedDict{UInt8,String}(
        0x01 => "VOLPE",
        0x02 => "FlightAware",
        0x03 => "web (flightaware.com)"
    )
    roots = ds.OrderedDict(v => k for (k,v) in pathdict["roots"])
    files = ds.OrderedDict(v => k for (k,v) in pathdict["files"])
    # Instantiate
    FlightSet{T}(volpe, flightaware, webdata,
        PrimaryMetadata{T}(altmin, (start=tmin, stop=tmax), sources, ds.OrderedDict(
            "roots" => roots, "files" => files
        ), tc, loadtime, attachments))
end # Main constructor FlightSet

#* Constructor for default Float32 FlightSet
FlightSet(args...; kwargs...) = FlightSet{Float32}(args...; kwargs...)

#* Constructor for empty FlightSet
FlightSet{T}() where T<:AbstractFloat = FlightSet{T}(FlightData{T}[], FlightData{T}[], FlightData{T}[], PrimaryMetadata{T}())

#* Constructor for floating point type promotion
FlightSet{T}(flights::FlightSet) where T<:AbstractFloat = FlightSet{T}(
    FlightData{T}.(flights.volpe),
    FlightData{T}.(flights.flightaware),
    FlightData{T}.(flights.webdata),
    PrimaryMetadata{T}(flights.metadata)
)

"""
    PrimarySet{T}(
        volpe::Vector{FlightData{T}},
        flightaware::Vector{FlightData{T}},
        webdata::Vector{FlightData{T}},
        metadata::PrimaryMetadata{T}
    ) where T<:AbstractFloat
    PrimarySet{T}(; kwargs...) where T<:AbstractFloat
    PrimarySet{T}(tracks::FlightSet) where T<:AbstractFloat

Alias constructors for `FlightSet{T}` instantiation.

See also [`FlightSet`](@ref), [`FlightData`](@ref), [`FlightTrack`](@ref), and [`PrimaryMetadata`](@ref).
"""
function PrimarySet end

PrimarySet{T}(
  volpe::Vector{FlightData{T}},
  flightaware::Vector{FlightData{T}},
  webdata::Vector{FlightData{T}},
  metadata::PrimaryMetadata{T}
) where T<:AbstractFloat = FlightSet{T}(volpe, flightaware, webdata, metadata)

#* Alias constructor for default Float32 PrimarySet
PrimarySet{T}(; kwargs...) where T<:AbstractFloat =
  FlightSet{T}(; kwargs...)

#* Alias constructor for type promotion of FlightSet
PrimarySet{T}(tracks::FlightSet) where T<:AbstractFloat = FlightSet{T}(tracks)
