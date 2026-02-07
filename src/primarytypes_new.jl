## Metadata structs

"""
# struct FlightMetadata{T} <: FlightTrack{T}

Immutable struct holding metadata for an individual `FlightTrack` of the `FlightSet` with fields

- `dbID::Union{Int,AbstractString}`
- `flightID::Union{Missing,AbstractString}`
- `route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}}`
- `aircraft::Union{Missing,AbstractString}`
- `date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}`
- `area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,AbstractFloat}}`
- `flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,AbstractFloat,AbstractFloat}}}}`
- `useLON::Bool`
- `source::AbstractString`
- `file::AbstractString`

By default, `Float32` is used for `T`.

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


# Instantiation
# TODO Update description

`FlightMetadata` is constructed automatically, when `FlightTrack` is instantiated using
a modified constructor and `dbID`, `flightID`, `aircraft` type, `route`, `useLON`,
`flex`, `source`, and `file`.
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
      source::AbstractString,
      file::AbstractString
    ) where T -> struct FlightMetadata

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
      source::AbstractString,
      file::AbstractString
    ) where T -> struct FlightMetadata
"""
struct FlightMetadata{T} <: FlightTrack{T}
    dbID::Union{Int,AbstractString}
    flightID::Union{Missing,AbstractString}
    route::Union{Missing,NamedTuple{(:orig,:dest),Tuple{AbstractString,AbstractString}}}
    aircraft::Union{Missing,AbstractString}
    date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
    area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,T}}
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}}
    useLON::Bool
    source::AbstractString
    file::AbstractString
    root::AbstractString
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
    source::AbstractString,
    file::AbstractString,
    root::AbstractString
) where T
    elonmin, elonmax = lonextrema(lon, ≥)
    wlonmin, wlonmax = lonextrema(lon, <)
    area = (latmin=minimum(lat), latmax=maximum(lat), elonmin, elonmax, wlonmin, wlonmax)
    new{T}(dbID, flightID, route, aircraft, (start=date[1], stop=date[end]), area,
        flex, useLON, source, file, root)
end #constructor 2 FlightMetadata

#* Constructor for empty FlightMetadata
FlightMetadata{T}() where T = FlightMetadata{T}(
    "", missing, missing, missing, (start=Dates.now(), stop=Dates.now()),
    (latmin=T(NaN), latmax=T(NaN), elonmin=T(NaN), elonmax=T(NaN), wlonmin=T(NaN), wlonmax=T(NaN)),
    ((range=0:0, min=T(NaN), max=T(NaN)),), false, "", "", ""
)

#* Constructor for default Float32 FlightMetadata
FlightMetadata(args...) = FlightMetadata{Float32}(args...)

#* Constructor for floating point type promotion
FlightMetadata{T}(meta::FlightMetadata) where T = FlightMetadata{T}(
    meta.dbID, meta.flightID, meta.route, meta.aircraft, meta.date,
    (latmin = T(meta.area.latmin), latmax = T(meta.area.latmax), elonmin = T(meta.area.elonmin),
    elonmax = T(meta.area.elonmax), wlonmin = T(meta.area.wlonmin), wlonmax = T(meta.area.wlonmax)),
    Tuple([(range = m.range, min = T.(m.min), max = T.(m.max)) for m in meta.flex]),
    meta.useLON, meta.source, meta.file, meta.root
)


"""
# struct PrimaryMetadata{T} <: PrimarySet{T}

Immutable struct with additional information of databases:

- `altmin::Real`: Minimum altitude threshold for which flight data is considered
- `date`: NamedTuple with entries `start`/`stop` giving the time range of the database
- `created`: time of creation of database
- `loadtime`: time it took to read data files and load it to the struct
- `remarks`: any additional data or comments that can be attached to the database
"""
struct PrimaryMetadata{T} <: PrimarySet{T}
    altmin::T
    date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
    created::Union{DateTime,ZonedDateTime}
    loadtime::Dates.CompoundPeriod
    remarks
end #struct PrimaryMetadata

#* Constructor for empty PrimaryMetadata
PrimaryMetadata{T}() where T = PrimaryMetadata{T}(NaN,
    (start=Dates.now(), stop=Dates.now()), Dates.now(), Dates.CompoundPeriod(), nothing)

#* Constructor for default Float32 PrimaryMetadata
PrimaryMetadata(args...) = PrimaryMetadata{Float32}(args...)

#* Constructor for floating point type promotion
PrimaryMetadata{T}(meta::PrimaryMetadata) where T = PrimaryMetadata{T}(T(meta.altmin),
    meta.date, meta.created, meta.loadtime, meta.remarks)


## Struct for single flight tracks

struct FlightData{T} <: FlightTrack{T}
    time::Vector{DateTime}
    lat::Vector{T}
    lon::Vector{T}
    alt::Vector{T}
    heading::Vector{<:Union{Missing,Int}}
    climb::Vector{<:Union{Missing,Int}}
    speed::Vector{T}
    metadata::FlightMetadata{T}

  """ Unmodified constructor for `FlightData` with basic checks for correct `data`"""
  function FlightData{T}(data::DataFrame, metadata::FlightMetadata{T}) where T
    # Ensure floats of correct precision
    convert_floats!(data, T)
    # Column checks and warnings
    standardnames = ["time", "lat", "lon", "alt", "heading", "climb", "speed"]
    standardtypes = [Union{DateTime,Vector{DateTime}},
      Vector{T}, Vector{T},
      Vector{<:Union{Missing,<:T}},
      Vector{<:Union{Missing,Int}},
      Vector{<:Union{Missing,<:T}},
      Vector{<:Union{Missing,<:T}}]
    bounds = (:lat => (-90, 90), :lon => (-180, 180), :alt => (0,Inf),
      :heading => (0, 360), :speed => (0, Inf))
    checkcols!(data, standardnames, standardtypes, bounds,
      metadata.source, metadata.dbID)
    new{T}(data,metadata)
  end #constructor 1 FlightData

  """ Modified constructor with variable checks and some automated calculation of fields """
  function FlightData{T}(
    track::DataFrame,
    dbID::Union{Int,AbstractString},
    flightID::Union{Missing,AbstractString},
    aircraft::Union{Missing,AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{<:AbstractString,<:AbstractString}}},
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}},
    useLON::Bool,
    source::String,
    file::AbstractString
  ) where T
    # Check dataframe columns of flight data; fill missing columns with missing values
    t = getproperty(track, :time)
    t = t isa Vector{ZonedDateTime} ? [zt.utc_datetime for zt in t] : t
    lat = getproperty(track, :lat)
    lon = getproperty(track, :lon)
    alt = getproperty(track, :alt)
    heading = df.hasproperty(track, :heading) ? df.getproperty(track, :heading) :
      [missing for i in t]
    climb = df.hasproperty(track, :climb) ? df.getproperty(track, :climb) :
      [missing for i in t]
    speed = df.hasproperty(track, :speed) ? df.getproperty(track, :speed) :
      [missing for i in t]
    metadata = FlightMetadata{T}(dbID,flightID,route,aircraft,t,lat,lon,useLON,flex,
      source,file)

    # Instatiate new FlightTrack
    new{T}(DataFrame(time=t,lat=lat,lon=lon,alt=alt,heading=heading,climb=climb,speed=speed),metadata)
  end #constructor 2 FlightData
end #struct FlightData


"""
    FlightData{T}() where T

External constructor for emtpy FlightData.
"""
FlightData{T}() where T = FlightData{T}(DataFrame(time = DateTime[],
  lat = T[], lon = T[], alt=T[], heading = Int[], climb = T[],
  speed = T[]), FlightMetadata{T}()
)

"""
    FlightData{T}(flight::FlightData) where T

External FlightData constructor for conversion of floating point precision.
"""
function FlightData{T}(flight::FlightData) where T
  convert_floats!(flight.data, T)
  FlightData{T}(flight.data, FlightMetadata{T}(flight.metadata))
end

"""
    FlightData(args...; kwargs...)

Default `FlightData` constructor for `Float32`.
"""
FlightData(args...; kwargs...) = FlightData{Float32}(args...; kwargs...)

"""
    FlightTrack(args...; kwargs...)

Alias constructor for default `Float32` `FlightData`.
"""
FlightTrack(args...; kwargs...) = FlightData{Float32}(args...; kwargs...)

"""
    FlightTrack{T}(args...; kwargs...) where T

Alias constructor for `FlightData{T}`.
"""
FlightTrack{T}(args...; kwargs...) where T = FlightData{T}(args...; kwargs...)


## Struct for sets of flight tracks

struct FlightSet{T} <: PrimarySet{T}
  volpe::Vector{FlightData{T}}
  flightaware::Vector{FlightData{T}}
  webdata::Vector{FlightData{T}}
  metadata::PrimaryMetadata{T}

    #* Internal constructor with data checks
    function FlightSet{T}(volpe::Vector{FlightData{T}}, flightaware::Vector{FlightData{T}},
        webdata::Vector{FlightData{T}}, metadata::PrimaryMetadata{T}) where T

        # Check for correct dataset type in each vector and for correct floating point precision
        volpe = checkDBtype(volpe, "VOLPE")
        flightaware = checkDBtype(flightaware, "FlightAware")
        webdata = checkDBtype(webdata, "flightaware.com")

        # Instantiate new struct
        new{T}(volpe, flightaware, webdata, metadata)
    end #constructor 1 FlightSet
end #struct FlightSet

  #* Modified constructor creating the database from an identifer of the
  #* database type and the respective folder path for that database.
  function FlightSet{T}(;
    volpe::Union{String,Vector{String}}=String[],
    flightaware::Union{String,Vector{String}}=String[],
    webdata::Union{String,Vector{String}}=String[],
    altmin::Real=5000,
    odelim::Union{Nothing,Char,String}=nothing,
    savedir::Union{String,Bool}="abs",
    remarks=nothing
  ) where T

    # Return empty FlightSet, if no folders are passed to constructor
    all(isempty.([volpe, flightaware, webdata])) &&
        return FlightSet{T}()
    # Save time of database creation
    tstart = Dates.now()

    ## Load databases for each type
    # VOLPE AEDT dataset
    volpe isa Vector || (volpe = [volpe])
    files = String[]
    for dir in volpe
        findfiles!(files, dir, ".csv")
    end
    files = convertdir.(files, savedir)
    volpe = loadVOLPE(files...; Float=T, altmin)
    # FlightAware commercial dataset
    flightaware isa Vector || (flightaware = [flightaware])
    files = String[]
    for dir in flightaware
        findfiles!(files, dir, ".csv")
    end
    files = convertdir.(files, savedir)
    flightaware = loadFA(files...; Float=T, altmin)
    # FlightAware web content
    webdata isa Vector || (webdata = [webdata])
    files = String[]
    for dir in webdata
        findfiles!(files, dir, ".tsv", ".txt", ".dat")
    end
    files = convertdir.(files, savedir)
    webdata = loadWD(files...; Float=T, altmin, delim=odelim)
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

    # wrap up and calculate load time
    tend = Dates.now()
    tc = tz.ZonedDateTime(tend, tz.localzone())
    loadtime = Dates.canonicalize(Dates.CompoundPeriod(tend - tstart))

    @info string("FlightSet loaded in ",
        "$(join(loadtime.periods[1:min(2,length(loadtime.periods))], ", ")) to",
        "\n▪ volpe ($(length(volpe)) entries)\n▪ flightaware ($(length(flightaware)) entries)\n",
        "▪ webdata ($(length(webdata)) entries)\n▪ metadata")

    # Instantiate
    new{T}(volpe, flightaware, webdata,
        PrimaryMetadata{T}(altmin, (start=tmin, stop=tmax), tc, loadtime, remarks))
  end # constructor 2 FlightSet



#* Main constructor
function FlightSet{T}(;
    volpe::Union{String,Vector{String}}=String[],
    flightaware::Union{String,Vector{String}}=String[],
    webdata::Union{String,Vector{String}}=String[],
    altmin::Real=5000,
    odelim::Union{Nothing,Char,String}=nothing,
    savedir::Union{String,Bool}="abs",
    remarks=nothing
) where T

    # Return empty FlightSet, if no folders are passed to constructor
    all(isempty.([volpe, flightaware, webdata])) &&
        return FlightSet{T}(FlightData{T}[], FlightData{T}[], FlightData{T}[], PrimaryMetadata{T}())
    # Save time of database creation
    tstart = Dates.now()

    ## Load databases for each type
    # VOLPE AEDT dataset
    volpe isa Vector || (volpe = [volpe])
    files = String[]
    for dir in volpe
        findfiles!(files, dir, ".csv")
    end
    files = convertdir.(files, savedir)
    volpe = loadVOLPE(files...; Float=T, altmin)
    # FlightAware commercial dataset
    flightaware isa Vector || (flightaware = [flightaware])
    files = String[]
    for dir in flightaware
        findfiles!(files, dir, ".csv")
    end
    files = convertdir.(files, savedir)
    flightaware = loadFA(files...; Float=T, altmin)
    # FlightAware web content
    webdata isa Vector || (webdata = [webdata])
    files = String[]
    for dir in webdata
        findfiles!(files, dir, ".tsv", ".txt", ".dat")
    end
    files = convertdir.(files, savedir)
    webdata = loadWD(files...; Float=T, altmin, delim=odelim)
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

    # wrap up and calculate load time
    tend = Dates.now()
    tc = tz.ZonedDateTime(tend, tz.localzone())
    loadtime = Dates.canonicalize(Dates.CompoundPeriod(tend - tstart))

    @info string("FlightSet loaded in ",
        "$(join(loadtime.periods[1:min(2,length(loadtime.periods))], ", ")) to",
        "\n▪ volpe ($(length(volpe)) entries)\n▪ flightaware ($(length(flightaware)) entries)\n",
        "▪ webdata ($(length(webdata)) entries)\n▪ metadata")

    # Instantiate
    new{T}(volpe, flightaware, webdata,
        PrimaryMetadata{T}(altmin, (start=tmin, stop=tmax), tc, loadtime, remarks))
end # Main constructor FlightSet

#* Constructor for default Float32 FlightSet
FlightSet(args...; kwargs...) = FlightSet{Float32}(args...; kwargs...)

#* Constructor for empty FlightSet
FlightSet{T}() where T = FlightSet{T}(FlightData{T}[], FlightData{T}[], FlightData{T}[], PrimaryMetadata{T}())

#* Constructor for floating point type promotion
FlightSet{T}(flights::FlightSet) where T = FlightSet{T}(
    FlightData{T}.(flights.volpe),
    FlightData{T}.(flights.flightaware),
    FlightData{T}.(flights.webdata),
    PrimaryMetadata{T}(flights.metadata)
)
