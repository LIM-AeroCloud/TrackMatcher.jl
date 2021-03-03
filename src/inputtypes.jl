## Define own Metadata structs
"""
# struct FlightMetadata{T}

Immutable struct to hold metadata for `FlightTrack` of the `FlightSet` with fields

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

By default, all `AbstractFloat` are set to Float32 but can be set to any other precision
in `FlightSet` with the kwarg `Float`.

## dbID
Database ID – integer counter for `inventory`,
String with information about `FlightID`, `route`, and/or scheduled arrival for
FlightAware data.

## FlightID and aircraft
Strings with aircraft identification and type.

## route
`NamedTuple` with fields for `orig`in and `dest`ination holding the airport codes.

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
- `"VOLPE AEDT"`
- `"FlightAware"`
- `"flightaware.com"`

## file
String holding the absolute folder path and file name.


# Instantiation

`FlightMetadata` is constructed automatically, when `FlightTrack` is instatiated using
a modified constructor and `dbID`, `flightID`, `aircraft` type, `route`, `useLON`,
`flex`, `source`, and `file`.
Fields `area` and `date` are calculated from `lat`/`lon`, and `date` vectors.

    FlightMetadata(
      dbID::Union{Int,AbstractString},
      flightID::Union{Missing,AbstractString},
      route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
      aircraft::Union{Missing,AbstractString},
      date::Vector{DateTime},
      lat::Vector{<:Union{Missing,<:AbstractFloat}},
      lon::Vector{<:Union{Missing,<:AbstractFloat}},
      useLON::Bool,
      flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,AbstractFloat,AbstractFloat}}}},
      source::AbstractString,
      file::AbstractString
    ) -> struct FlightMetadata

Or construct `FlightMetadata` by directly handing over every field:

    FlightMetadata(
      dbID::Union{Int,AbstractString},
      flightID::Union{Missing,AbstractString},
      route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
      aircraft::Union{Missing,AbstractString},
      date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
      area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,AbstractFloat}},
      flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,AbstractFloat,AbstractFloat}}}},
      useLON::Bool,
      source::AbstractString,
      file::AbstractString
    ) -> struct FlightMetadata
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

  """ Unmodified constructor for `FlightMetadata` """
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
  ) where T
    # T = typeof(area.latmin)
    new{T}(dbID, flightID, route, aircraft, date, area, flex, useLON, source, file)
  end #constructor 1 FlightMetadata

  """
  Modified constructor for FlightMetadata with some automated construction of fields
  and variable checks.
  """
  function FlightMetadata{T}(
    dbID::Union{Int,AbstractString},
    flightID::Union{Missing,AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
    aircraft::Union{Missing,AbstractString},
    date::Vector{DateTime},
    lat::Vector{<:Union{Missing,<:AbstractFloat}},
    lon::Vector{<:Union{Missing,<:AbstractFloat}},
    useLON::Bool,
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}},
    source::AbstractString,
    file::AbstractString
  ) where T
    # T = promote_type(eltype(lat), eltype(lon))
    elonmax = isempty(lon[lon.≥0]) ? T(NaN) : maximum(lon[lon.≥0])
    elonmin = isempty(lon[lon.≥0]) ? T(NaN) : minimum(lon[lon.≥0])
    wlonmax = isempty(lon[lon.<0]) ? T(NaN) : maximum(lon[lon.<0])
    wlonmin = isempty(lon[lon.<0]) ? T(NaN) : minimum(lon[lon.<0])
    area = (latmin=minimum(lat), latmax=maximum(lat),
      elonmin=elonmin, elonmax=elonmax, wlonmin=wlonmin, wlonmax=wlonmax)
    new{T}(dbID, flightID, route, aircraft, (start=date[1], stop=date[end]), area,
      flex, useLON, source, file)
  end #constructor 2 FlightMetadata
end #struct FlightMetadata

""" External constructor for emtpy FlightMetadata struct """
FlightMetadata{T}() where T = FlightMetadata("", missing, missing, missing,
  (start=Dates.now(), stop=Dates.now()), (latmin=T(NaN), latmax=T(NaN),
  elonmin=T(NaN), elonmax=T(NaN), wlonmin=T(NaN), wlonmax=T(NaN)),
  ((range=0:0, min=T(NaN), max=T(NaN)),), true, "","")


"""
# struct CloudMetadata

Currently only place holder for remarks available during test phase.
"""
struct CloudMetadata{T} <: CloudTrack{T}
  ID::String
  date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,T}}
  flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}}
  useLON::Bool
  file::String

  """ Unmodified constructor for `CloudMetadata` """
  function CloudMetadata(
    ID::String,
    date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,T}},
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}},
    useLON::Bool,
    file::String
  ) where T<:AbstractFloat
    # T = typeof(area.latmin)
    new{T}(ID, date, area, flex, useLON, file)
  end #constructor 1 CloudMetadata

  """
  Modified constructor for CloudMetadata with some automated construction of fields
  and variable checks.
  """
  function CloudMetadata(
    ID::Union{Int,AbstractString},
    data::DataFrame,
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}},
    useLON::Bool,
    file::AbstractString
  ) where T<:AbstractFloat
    # T = promote_type(eltype(data.lat), eltype(data.lon))
    elonmax = isempty(data.lon[data.lon.≥0]) ? T(NaN) : maximum(data.lon[data.lon.≥0])
    elonmin = isempty(data.lon[data.lon.≥0]) ? T(NaN) : minimum(data.lon[data.lon.≥0])
    wlonmax = isempty(data.lon[data.lon.<0]) ? T(NaN) : maximum(data.lon[data.lon.<0])
    wlonmin = isempty(data.lon[data.lon.<0]) ? T(NaN) : minimum(data.lon[data.lon.<0])
    area = (latmin=minimum(data.lat), latmax=maximum(data.lat),
      elonmin=elonmin, elonmax=elonmax, wlonmin=wlonmin, wlonmax=wlonmax)
    new{T}(ID, (start=data.time[1], stop=data.time[end]), area, flex, useLON, file)
  end #constructor 2 CloudMetadata
end #struct CloudMetadata


"""
# struct SatMetadata

Immutable struct to hold metadata for `SatData` with fields

- `files::Dict{Int,String}`
- `type::Symbol`
- `date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}`
- `created::Union{DateTime,ZonedDateTime}`
- `loadtime::Dates.CompoundPeriod`
- `remarks`

## files
Dictionary with indices of `fileindex` column in `SatData.data` pointing to the
full file names.

## type
Symbol indicating, whether profile or layer data is stored.

## date
`NamedTuple` with fields `start` and `stop` for start and end time of the monitored
satellite period.

## created
time of creation of database

## loadtime
time it took to read data files and load it to the struct

##remarks
any additional data or comments that can be attached to the database


# Instantiation

`SatMetadata` is constructed automatically, when `SatData` is instatiated using
a modified constructor and `files`, `date`, `loadtime`, and `remarks`.

    function SatMetadata(
      files::Vector{String},
      date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
      loadtime::Dates.CompoundPeriod=Dates.canonicalize(Dates.CompoundPeriod());
      remarks=nothing
    ) -> struct SattMetadata
"""
struct SatMetadata
  files::Dict{Int,String}
  type::Symbol
  date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  created::Union{DateTime,ZonedDateTime}
  loadtime::Dates.CompoundPeriod
  remarks

  function SatMetadata(
    files::Dict{Int,String},
    type::Symbol,
    date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    created::Union{DateTime,ZonedDateTime},
    loadtime::Dates.CompoundPeriod,
    remarks=nothing
  )
    new(files, type, date, created, loadtime, remarks)
  end #constructor 1 SatMetadata

  function SatMetadata(
    files::Vector{String},
    date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    loadtime::Dates.CompoundPeriod=Dates.canonicalize(Dates.CompoundPeriod());
    remarks=nothing
  )
    # Find type of satellite data based on first 50 files (~2 days)
    type = occursin("CLay", files[1]) ≥
      count(occursin.("CPro", files[1:min(length(files), 50)])) ? :CLay : :CPro
    # Create a new instance of SatMetadata
    new(Dict(enumerate(files)), type, date, tz.now(tz.localzone()), loadtime, remarks)
  end #constructor 2 SatMetadata
end #struct SatMetadata


"""
# struct SetMetadata

Immutable struct with additional information of databases:

- `altmin::Real`: Minimum altitude threshold for which flight data is considered
- `date`: NamedTuple with entries `start`/`stop` giving the time range of the database
- `created`: time of creation of database
- `loadtime`: time it took to read data files and load it to the struct
- `remarks`: any additional data or comments that can be attached to the database
"""
struct SetMetadata{T} <: PrimarySet{T}
  altmin::T
  date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  created::Union{DateTime,ZonedDateTime}
  loadtime::Dates.CompoundPeriod
  remarks
end #struct SetMetadata


## Define structs related to flight data

"""
# struct FlightTrack{T<:AbstractFloat}

Aircraft data with fields
- `data::DataFrame`
- `metadata::FlightMetadata`

The `DataFrame` of `data` has columns in the following order with the respective types:

- `time::Vector{DateTime}`                        (time stamp of measurement)
- `lat::Vector{<:Union{Missing,<:AbstractFloat}}`   (latitude in deg)
- `lon::Vector{<:Union{Missing,<:AbstractFloat}}`   (longitude in deg)
- `alt::Vector{<:Union{Missing,<:AbstractFloat}}`   (altitude in meters)
- `heading::Vector{<:Union{Missing,Int}}`         (heading/direction of the aircraft in deg)
- `climb::Vector{<:Union{Missing,Int}}`           (climbing (positive)/sinking (negative values) rate in m/s)
- `speed::Vector{<:Union{Missing,<:AbstractFloat}}` (velocity in m/s)

By default, all AbstractFloat are set to `Float32`, but can be set to any other
precision in `FlightSet` by the kwarg `Float`.


# Instantiation

    FlightTrack(
      FlightTrack::DataFrame,
      dbID::Union{Int,AbstractString},
      flightID::Union{Missing,AbstractString},
      aircraft::Union{Missing,AbstractString},
      route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
      flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,<:AbstractFloat,<:AbstractFloat}}}},
      useLON::Bool,
      source::String,
      file::AbstractString
    ) -> struct FlightTrack

Construct `FlightTrack` from the `data` `DataFrame` and additonal metainformation
`dbID`, `flightID`, `aircraft` type, `route`, `flex` points in the flight trajectory,
a flag whether longitudes are used as x data in the track (`useLON`), the database
`source`, and `file` names.
For `data`, `time` can be passed `ZonedDateTime`, which will be converted to `UTC`
standard time, or as UTC `DateTime`.

Or construct by directly handing over every field:

    FlightTrack(data::DataFrame, metadata::FlightMetadata)

Checks exist that the order, names, and types of the `data` `DataFrame` are correct.
"""
struct FlightData{T} <: FlightTrack{T}
  data::DataFrame
  metadata::FlightMetadata

  """ Unmodified constructor for `FlightData` with basic checks for correct `data`"""
  function FlightData{T}(data::DataFrame, metadata::FlightMetadata{T}) where T<:AbstractFloat
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
end #struct FlightData


""" Modified constructor with variable checks and some automated calculation of fields """
function FlightTrack{T}(
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
  # T = promote_type(eltype(lon), eltype(lat))
  metadata = FlightMetadata{T}(dbID,flightID,route,aircraft,t,lat,lon,useLON,flex,source,file)

  # Instatiate new FlightTrack
  return FlightData{T}(DataFrame(time=t,lat=lat,lon=lon,alt=alt,heading=heading,climb=climb,speed=speed),metadata)
end #constructor 2 FlightTrack


""" External constructor for emtpy FlightTrack struct """
FlightTrack{T}() where T = FlightData(DataFrame(time = DateTime[],
  lat = T[], lon = T[], alt=T[], heading = Int[], climb = T[],
  speed = T[]), FlightMetadata{T}())


"""
# struct FlightSet

Database for aircraft data of different database types with fields:
- `inventory::Vector{FlightTrack}`
- `archive::Vector{FlightTrack}`
- `onlineData::Vector{FlightTrack}`
- `metadata`: `DBMetatadata` with information about
  - `altmin::Real`: Minimum altitude threshold for which flight data is considered
  - `date`: date range of the dataset
  - `created`: time of creation
  - `loadtime`
  - `remarks` (any additional data or comments)

## inventory
Flight data from csv files.

## archive
Commercial flight data by FlightAware from a csv file.

## onlineData
Online data from the FlightAware website copied to whitespace-separated files.

## metadata
Immutable struct `SetMetadata`.


# Instantiation

Instatiate by giving a String with identifiers of the `DBtype` and an equal number
of `folder` paths as characters in the `DBtype` `String`. Optionally add a minimum
altitude threshold for the data (default = `15000`) and any remarks
(comments or additional data). Define the delimiter in the input files of the
online data with the keyword `odelim`. Use any character or string as delimiter.
By default (`odelim=nothing`), auto-detection is used.

    FlightSet(DBtype::String, folder::Union{String, Vector{String}}...;
      Float::DataType=Float32, altmin::Real=5_000, remarks=nothing,
      odelim::Union{Nothing,Char,String}=nothing)

`DBtype` can be identified with:
- `1` or `i`: VOLPE AEDT inventory
- `2` or `a`: FlightAware archived data (commercially available)
- `3` or `o`: flightaware.com online data

By default, all values are read in as `Float32`, but can be set to any other
precision by the `Float` kwarg; `altmin` set the minimum threshold above which
flight tracking points are considered. Set the delimiter of the input files with
kwarg `odelim`, if delimiters are any character different from whitespace. Any
`remarks` can be attached to `FlightSet`.

Alternatively, instatiate directly with the fields of `FlightSet`, where the correct
database type is checked, and wrong datasets are removed in every field.
"""
struct FlightSet{T} <: PrimarySet{T}
  inventory::Vector{FlightData{T}}
  archive::Vector{FlightData{T}}
  onlineData::Vector{FlightData{T}}
  metadata::SetMetadata{T}

  """
  Unmodified constructor for `FlightSet` with basic checks for correct dataset type
  in each dataset field.
  """
  function FlightSet{T}(inventory::Vector{FlightData{T}}, archive::Vector{FlightData{T}},
    onlineData::Vector{FlightData{T}}, metadata::SetMetadata{T}) where {T}

    # Check for correct dataset type in each vector and for correct floating point precision
    inventory = checkDBtype(inventory, "VOLPE")
    archive = checkDBtype(archive, "FlightAware")
    onlineData = checkDBtype(onlineData, "flightaware.com")

    new{T}(inventory, archive, onlineData, metadata)
  end #constructor 1 FlightSet

  """
  Modified constructor creating the database from an identifer of the
  database type and the respective folder path for that database.
  """
  function FlightSet{T}(DBtype::String, folder::String...; altmin::Real=5000,
    remarks=nothing, odelim::Union{Nothing,Char,String}=nothing) where T

    # Save time of database creation
    tstart = Dates.now()
    # Check DBtype addresses all folder paths
    if length(DBtype) ≠ length(folder)
      throw(ArgumentError("Number of characters in `DBtype` must match length of vararg `folder`"))
    end
    # Find database types
    i1 = [findall(isequal('i'), lowercase(DBtype)); findall(isequal('1'), DBtype)]
    i2 = [findall(isequal('a'), lowercase(DBtype)); findall(isequal('2'), DBtype)]
    i3 = [findall(isequal('o'), lowercase(DBtype)); findall(isequal('3'), DBtype)]

    # Load databases for each type
    # VOLPE AEDT inventory
    ifiles = String[]
    for i in i1
      findfiles!(ifiles, folder[i], ".csv")
    end
    inventory = loadInventory(ifiles...; Float=T, altmin=altmin)
    # FlightAware commercial archive
    ifiles = String[]
    for i in i2
      findfiles!(ifiles, folder[i], ".csv")
    end
    archive = loadArchive(ifiles...; Float=T, altmin=altmin)
    ifiles = String[]
    for i in i3
      findfiles!(ifiles, folder[i], ".tsv", ".txt", ".dat")
    end
    onlineData = loadOnlineData(ifiles...; Float=T, altmin=altmin, delim=odelim)
    tmin, tmax = if isempty([inventory; archive; onlineData])
      tstart, tstart
    else
      minimum([[f.metadata.date.start for f in inventory];
        [f.metadata.date.start for f in archive];
        [f.metadata.date.start for f in onlineData]]),
      maximum([[f.metadata.date.stop for f in inventory];
        [f.metadata.date.stop for f in archive];
        [f.metadata.date.stop for f in onlineData]])
    end

    tend = Dates.now()
    tc = tz.ZonedDateTime(tend, tz.localzone())
    loadtime = Dates.canonicalize(Dates.CompoundPeriod(tend - tstart))

    @info string("FlightSet loaded in ",
      "$(join(loadtime.periods[1:min(2,length(loadtime.periods))], ", ")) to",
      "\n▪ inventory ($(length(inventory)) entries)\n▪ archive ($(length(archive)) entries)\n",
      "▪ onlineData ($(length(onlineData)) entries)\n▪ metadata")

    new{T}(inventory, archive, onlineData,
      SetMetadata{T}(altmin, (start=tmin, stop=tmax), tc, loadtime, remarks))
  end # constructor 2 FlightSet
end #struct FlightSet


""" Default FlightSet/FlightData constructor for Float32 """
FlightSet(DBtype::String, folder::String...; altmin::Real=5000,
  remarks=nothing, odelim::Union{Nothing,Char,String}=nothing) =
  FlightSet{Float32}(DBtype, folder...; altmin=altmin, remarks=remarks, odelim=odelim)


## Define structs related to cloud data

"""
# struct CloudTrack{T<:AbstractFloat}

Store Data related to a single cloud track.

## Fields

- `time::Vector{DateTime}`
- `lat::Vector{T}`
- `lon::Vector{T}`
- `metadata::CloudMetadata`

Default floating point precision is `Float32`.
"""
struct CloudData{T} <: PrimaryTrack{T}
  data::DataFrame
  metadata::CloudMetadata

  """ Unmodified constructor for `CloudTrack` with basic checks for correct `data`"""
  function CloudData{T}(data::DataFrame, metadata::CloudMetadata) where T

    # Column checks and warnings
    standardnames = ["time", "lat", "lon"]
    standardtypes = [Union{DateTime,Vector{DateTime}},
      Vector{T}, Vector{T}]
    bounds = (:lat => (-90, 90), :lon => (-180, 180))
    checkcols!(data, standardnames, standardtypes, bounds,
      "CloudTrack", metadata.ID)
    new{T}(data,metadata)
  end #constructor 1 CloudTrack
end #struct CloudTrack


""" CloudTrack constructor for CloudData """
CloudTrack{T}(data::DataFrame, metadata::CloudMetadata) where T =
  CloudData{T}(data, metadata)

""" CloudTrack default Float32 constructor for CloudData """
CloudTrack(data::DataFrame, metadata::CloudMetadata) =
  CloudData{Float32}(data, metadata)



"""
# struct CloudSet{T} <: PrimarySet{T}

Database for cloud track data with fields:
- `tracks::Vector{CloudTrack}`
- `metadata::SetMetadata`
"""
struct CloudSet{T} <: PrimarySet{T}
  tracks::Vector{CloudTrack{T}}
  metadata::SetMetadata{T}

  """ unmodified constructor for CloudSet """
  CloudSet{T}(tracks::Vector{CloudTrack}, metadata::SetMetadata) where T = new{T}(tracks, metadata)

  """
  Modified constructor creating the database from mat files in the given folder
  or any subfolder using the floating point precision given by `Float`.
  """
  function CloudSet{T}(folders::String...; remarks=nothing) where T
    tstart = Dates.now()
    # Scan folders for HDF4 files
    files = String[]
    for folder in folders
      try findfiles!(files, folder, ".mat")
      catch
        @warn "read error; data skipped" folder
      end
    end

    # Load cloud tracks from mat files into TrackMatcher in Julia format
    tracks = loadCloudTracks(files...; Float=T)
    # Calculate load time
    tend = Dates.now()
    tc = tz.ZonedDateTime(tend, tz.localzone())
    loadtime = Dates.canonicalize(Dates.CompoundPeriod(tend - tstart))
    # For now find min/max times in CloudTracks
    tmin = minimum(t.data.time[1] for t in tracks)
    tmax = maximum(t.data.time[end] for t in tracks)

    # Instantiate CloudDB
    new{T}(tracks, SetMetadata{T}(NaN, (start=tmin, stop=tmax), tc, loadtime, remarks))
  end #modified constructor 2
end #struct CloudDB


""" Default CloudSet construct for Float32 """
CloudSet(folders::String...; remarks=nothing) = CloudSet{Float32}(folders...; remarks=remarks)


## Define structs related to sat data

"""
# struct SatData

Satellite data with fields
- `data::DataFrame`
- `metadata::SatMetadata`

The `DataFrame` of `data` has columns for the satellite time and position and indices
for data files with additional information:
`
- `time::Vector{DateTime}`                        (time stamp of measurement)
- `lat::Vector{<:Union{Missing,<:AbstractFloat}}`   (latitude)
- `lon::Vector{<:Union{Missing,<:AbstractFloat}}`   (longitude)
- `fileindex::Vector{Int}`                        (index for filenames in `SatMetadata`)


# Instantiation

    SatData(folders::String...; Float::DataType=Float32, type::Symbol=:undef,
      remarks=nothing) -> struct SatData

Construct `SatData` from any number of absolute or relative folder paths given as string.
SatData searches for hdf files in all folders recursively and determines the data type
(`CLay` or `CPro`) from the file names of the first 50 files unless the type is specified
with a Symbol `:CLay` or `:CPro` by the keyword argument `type`. Only one type of satellite
data can be stored in `SatData`. By default, all values are read in as `Float32`,
but can be set to any other precision by the `Float` kwarg.
"""
struct SatData
  data::DataFrame
  metadata::SatMetadata

  """ Unmodified constructor for `SatData` with basic checks for correct `data`"""
  function SatData(data::DataFrame, metadata::SatMetadata)
    standardnames = ["time", "lat", "lon", "fileindex"]
    standardtypes = [Vector{DateTime}, Vector{<:AbstractFloat},
      Vector{<:AbstractFloat}, Vector{Int}]
    bounds = (:time => (DateTime(2006), Dates.now()), :lat => (-90,90),
      :lon => (-180,180), :fileindex => (1, Inf))
    checkcols!(data, standardnames, standardtypes, bounds, "CLay")
    new(data, metadata)
  end #constructor 1 SatData

  """
  Modified constructor creating the database from mat files in the given `folders`
  or any subfolder using the floating point precision given by `Float`. The sat data
  `type` is determined from the first 50 files in the database unless directly
  specified `type`. Any `remarks` can be added to the metadata.
  """
  function SatData(
    folders::String...;
    Float::DataType=Float32,
    type::Symbol=:undef,
    remarks=nothing
  )
    tstart = Dates.now()
    # Scan folders for HDF4 files
    files = String[]
    for folder in folders
      try findfiles!(files, folder, ".hdf")
      catch
        @warn "read error; data skipped" folder
      end
    end
    # If type of satellite data is not defined, find it based on first 50 file names (~2 days)
    type = type == :CLay || type == :CPro ? string(type) :
      count(occursin.("CLay", files[1:min(length(files), 50)])) ≥
      count(occursin.("CPro", files[1:min(length(files), 50)])) ? "CLay" : "CPro"
    # Select files of type based on file name
    satfiles = files[occursin.(type, basename.(files))]
    wrongtype = length(files) - length(satfiles)
    wrongtype > 0 &&
      @warn "$wrongtype files with wrong satellite type detected; data skipped"
    # Create empty struct, if no data files were found
    isempty(satfiles) && return SatData(Float)
    # Start MATLAB session
    ms = mat.MSession()
    # Initialise arrays
    utc = Vector{Vector{DateTime}}(undef, length(satfiles))
    lat = Vector{Vector{Float}}(undef, length(satfiles))
    lon = Vector{Vector{Float}}(undef, length(satfiles))
    fileindex = Vector{Vector{Int}}(undef, length(satfiles))
    # Loop over files
    prog = pm.Progress(length(satfiles), "load sat data...")
    for (i, file) in enumerate(satfiles)
      # Find files with cloud layer data
      utc[i], lat[i], lon[i], fileindex[i] = try
        # Extract time
        mat.put_variable(ms, :file, file)
        mat.eval_string(ms, "clear t\ntry\nt = hdfread(file, 'Profile_UTC_Time');\nend")
        t = mat.jarray(mat.get_mvariable(ms, :t))[:,2]
        # Extract lat/lon
        mat.eval_string(ms, "clear longitude\ntry\nlongitude = hdfread(file, 'Longitude');\nend")
        longitude = mat.jarray(mat.get_mvariable(ms, :longitude))[:,2]
        mat.eval_string(ms, "clear latitude\ntry\nlatitude = hdfread(file, 'Latitude');\nend")
        latitude = mat.jarray(mat.get_mvariable(ms, :latitude))[:,2]
        # Save time converted to UTC and lat/lon
        convertUTC.(t), latitude, longitude, [i for index in t]
      catch
        # Skip data on failure and warn
        @warn "read error in CALIPSO granule; data skipped"  granule = splitext(basename(file))[1]
        DateTime[], Float[], Float[], Int[]
      end
      # Monitor progress for progress bar
      pm.next!(prog, showvalues = [(:date,Dates.Date(splitdir(dirname(file))[2], "y_m_d"))])
    end #loop over files
    pm.finish!(prog)

    # Close MATLAB session
    mat.close(ms)
    # Calculate time span of satellite data
    sattime = [utc...;]
    # Find data range
    tmin = minimum(sattime)
    tmax = maximum(sattime)
    # Save computing times
    tend = Dates.now()
    loadtime = Dates.canonicalize(Dates.CompoundPeriod(tend - tstart))

    # Instantiate new struct
    @info string("SatData data loaded in ",
      "$(join(loadtime.periods[1:min(2,length(loadtime.periods))], ", ")) to",
      "\n▪ data ($(length(sattime)) data rows)\n  – time\n  – lat\n  – lon",
      "\n  – fileindex\n▪ metadata")
    satdata = DataFrame(time=sattime, lat=[lat...;],
      lon=[lon...;], fileindex=[fileindex...;])
    new(satdata, SatMetadata(satfiles, (start=tmin, stop=tmax), loadtime, remarks=remarks))
  end #constructor 2 SatData
end #struct SatData


""" External constructor for emtpy SatData struct """
SatData(Float::DataType=Float32) = SatData(DataFrame(time=DateTime[], lat=Float[],
  lon=Float[], fileindex=Int[]), SatMetadata(Dict{Int,String}(), :undef,
  (start=Dates.now(), stop=Dates.now()), Dates.now(),
  Dates.canonicalize(Dates.CompoundPeriod())))
