"""
# Module TrackMatcher

Find intersections between different trajectories. The module is aimed to find
intersections between aircraft and satellite tracks, but can be modified for use
with ship or cloud tracks.

## Exported structs

- `FlightDB` stores flight track data and other relevant aircraft related data
  from 3 different inventories:
  - `inventory`: VOLPE AEDT inventory
  - `archive`: commercially available database by FlightAware
  - `onlineData`: free online data by FlightAware
- `DBMetadata` stores metadata for the flight database (`FLightDB`)
- `FlightData` stores data of a single flight in `FlightDB`
- `FlightMetadata` holds metadata to every flight
- `SatData` stores CALIPSO position and time and a file index for the granule of each data line
- `CLay` CALIPSO cloud layer data
- `CPro` CALIPSO cloud profile data
- `SatMetadata` stores metadata of the CALIPSO data
- `Intersection`: positions of intersections between aircraft and satellite trajectories,
  together with time difference between crossings; and `FlightData` and `CPro`/`CLay` in the
  vicinity of the intersection as well as information about the accuracy of the data
- `XMetadata` stores metadata for the `Intersection` data
"""
module TrackMatcher

## Import Julia packages
import DataFrames; const df = DataFrames
import CSV
import Dates
import TimeZones; const tz = TimeZones
import Distances; const dist = Distances
import MATLAB; const mat = MATLAB
import Statistics; const stats = Statistics
import ProgressMeter; const pm = ProgressMeter
import Logging; const logg = Logging

# Import structs and functions from packages
import PCHIP: Polynomial, pchip, interpolate
import DataFrames.DataFrame
import Dates: DateTime, Date, Time
import TimeZones.ZonedDateTime

# Define Logger with log level
logger = try logg.SimpleLogger(logfile, logg.Debug)
catch; logg.ConsoleLogger(stdout, logg.Debug)
end
logg.global_logger(logger)


## Define own Metadata structs
"""
# struct FlightMetadata{T<:AbstractFloat}

Immutable struct to hold metadata for `FlightData` of the `FlightDB` with fields

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
in `FlightDB` with the kwarg `Float`.

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

`FlightMetadata` is constructed automatically, when `FlightData` is instatiated using
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
struct FlightMetadata{T<:AbstractFloat}
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
  function FlightMetadata(
    dbID::Union{Int,AbstractString},
    flightID::Union{Missing,AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
    aircraft::Union{Missing,AbstractString},
    date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,T}} where T<:AbstractFloat,
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}} where T<:AbstractFloat},
    useLON::Bool,
    source::AbstractString,
    file::AbstractString
  )
    T = typeof(flex.latmin)
    new{T}(dbID, flightID, route, aircraft, date, area, flex, useLON, source, file)
  end #constructor 1 FlightMetadata

  """
  Modified constructor for FlightMetadata with some automated construction of fields
  and variable checks.
  """
  function FlightMetadata(
    dbID::Union{Int,AbstractString},
    flightID::Union{Missing,AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
    aircraft::Union{Missing,AbstractString},
    date::Vector{DateTime},
    lat::Vector{<:Union{Missing,<:AbstractFloat}},
    lon::Vector{<:Union{Missing,<:AbstractFloat}},
    useLON::Bool,
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}} where T<:AbstractFloat,
    source::AbstractString,
    file::AbstractString
  )
    T = eltype(lat)
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
# struct DBMetadata

Immutable struct with additional information of databases:

- `altmin::Real`: Minimum altitude threshold for which flight data is considered
- `date`: NamedTuple with entries `start`/`stop` giving the time range of the database
- `created`: time of creation of database
- `loadtime`: time it took to read data files and load it to the struct
- `remarks`: any additional data or comments that can be attached to the database
"""
struct DBMetadata
  altmin::Real
  date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  created::Union{DateTime,ZonedDateTime}
  loadtime::Dates.CompoundPeriod
  remarks
end #struct DBMetadata


"""
# struct XMetadata

Immutable struct with additional information of intersection data:

- `maxtimediff`: maximum time difference in minutes allowed between
  satellite overpass and aircraft passing at intersection
- `stepwidth`: stepwidth (in meters, partially internally converted to degrees at equator)
  used to interpolate track data
- `Xradius`: radius in meters around an intersection in which further intersections
  will be removed as duplicates due to the interpolation algorithm
- `lidarrange::NamedTuple{(:top,:bottom),Tuple{Real,Real}}`: user defined level thresholds
  for top/bottom heights for which CALIPSO data is considered
- `lidarprofile::NamedTuple`: CALIPSO lidar height levels used in the current dataset
- `sattype::Symbol`: cloud layer or profile data types used for the satellite date
- `satdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}`: date range of
  satellite data
- `altmin::Real`: minimum threshold for which flight data is considered
- `flightdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}`: date range
  of flight data
- `created`: time of creation of database
- `loadtime`: time it took to find intersections and load it to the struct
- `remarks`: any additional data or comments that can be attached to the database

XMetadata can be instantiated using a `Tuple` or `NamedTuple` for the `lidarrange`
"""
struct XMetadata
  maxtimediff::Int
  stepwidth::Real
  Xradius::Real
  lidarrange::NamedTuple{(:top,:bottom),Tuple{Real,Real}}
  lidarprofile::NamedTuple
  sattype::Symbol
  satdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  altmin::Real
  flightdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  created::Union{DateTime,ZonedDateTime}
  loadtime::Dates.CompoundPeriod
  remarks

  function XMetadata(
    maxtimediff::Int,
    stepwidth::Real,
    Xradius::Real,
    lidarrange::NamedTuple{(:top,:bottom),Tuple{Real,Real}},
    lidarprofile::NamedTuple,
    sattype::Symbol,
    satdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    altmin::Real,
    flightdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    created::Union{DateTime,ZonedDateTime},
    loadtime::Dates.CompoundPeriod,
    remarks
  )
    new(maxtimediff, stepwidth, Xradius,lidarrange, lidarprofile,
      sattype, satdates, altmin, flightdates, created, loadtime, remarks)
  end #constructor 1 XMetaData

  function XMetadata(
    maxtimediff::Int,
    stepwidth::Real,
    Xradius::Real,
    lidarrange::Tuple{Real,Real},
    lidarprofile::NamedTuple,
    sattype::Symbol,
    satdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    altmin::Real,
    flightdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    created::Union{DateTime,ZonedDateTime},
    loadtime::Dates.CompoundPeriod,
    remarks=nothing
  )
    new(maxtimediff, stepwidth, Xradius,(top=lidarrange[1], bottom=lidarrange[2]),
      lidarprofile, sattype, satdates, altmin, flightdates, created, loadtime, remarks)
  end #constructor 2 XMetaData
end #struct XMetaData


## Define structs related to flight data
"""
# struct FlightData{T<:AbstractFloat}

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
precision in `FlightDB` by the kwarg `Float`.


# Instantiation

    FlightData(
      flightdata::DataFrame,
      dbID::Union{Int,AbstractString},
      flightID::Union{Missing,AbstractString},
      aircraft::Union{Missing,AbstractString},
      route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
      flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,<:AbstractFloat,<:AbstractFloat}}}},
      useLON::Bool,
      source::String,
      file::AbstractString
    ) -> struct FlightData

Construct `FlightData` from the `data` `DataFrame` and additonal metainformation
`dbID`, `flightID`, `aircraft` type, `route`, `flex` points in the flight trajectory,
a flag whether longitudes are used as x data in the track (`useLON`), the database
`source`, and `file` names.
For `data`, `time` can be passed `ZonedDateTime`, which will be converted to `UTC`
standard time, or as UTC `DateTime`.

Or construct by directly handing over every field:

    FlightData(data::DataFrame, metadata::FlightMetadata)

Checks exist that the order, names, and types of the `data` `DataFrame` are correct.
"""
struct FlightData{T<:AbstractFloat}
  data::DataFrame
  metadata::FlightMetadata

  """ Unmodified constructor for `FlightData` with basic checks for correct `data`"""
  function FlightData(data::DataFrame, metadata::FlightMetadata)

    # Column checks and warnings
    standardnames = ["time", "lat", "lon", "alt", "heading", "climb", "speed"]
    standardtypes = [Union{DateTime,Vector{DateTime}},
      Vector{<:AbstractFloat}, Vector{<:AbstractFloat},
      Vector{<:Union{Missing,<:AbstractFloat}},
      Vector{<:Union{Missing,Int}},
      Vector{<:Union{Missing,<:AbstractFloat}},
      Vector{<:Union{Missing,<:AbstractFloat}}]
    bounds = (:lat => (-90, 90), :lon => (-180, 180), :alt => (0,Inf),
      :heading => (0, 360), :speed => (0, Inf))
    checkcols!(data, standardnames, standardtypes, bounds,
      metadata.source, metadata.dbID)
    T = eltype(data.lon)
    new{T}(data,metadata)
  end #constructor 1 FlightData

  """ Modified constructor with variable checks and some automated calculation of fields """
  function FlightData(
    flightdata::DataFrame,
    dbID::Union{Int,AbstractString},
    flightID::Union{Missing,AbstractString},
    aircraft::Union{Missing,AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}} where T<:AbstractFloat}},
    useLON::Bool,
    source::String,
    file::AbstractString
  )
    # Check dataframe columns of flight data; fill missing columns with missing values
    t = getproperty(flightdata, :time)
    t = t isa Vector{ZonedDateTime} ? [zt.utc_datetime for zt in t] : t
    lat = getproperty(flightdata, :lat)
    lon = getproperty(flightdata, :lon)
    alt = getproperty(flightdata, :alt)
    heading = df.hasproperty(flightdata, :heading) ? df.getproperty(flightdata, :heading) :
      [missing for i in t]
    climb = df.hasproperty(flightdata, :climb) ? df.getproperty(flightdata, :climb) :
      [missing for i in t]
    speed = df.hasproperty(flightdata, :speed) ? df.getproperty(flightdata, :speed) :
      [missing for i in t]
    T = eltype(lon)
    metadata = FlightMetadata(dbID,flightID,route,aircraft,t,lat,lon,useLON,flex,source,file)

    # Instatiate new FlightData
    new{T}(DataFrame(time=t,lat=lat,lon=lon,alt=alt,heading=heading,climb=climb,speed=speed),metadata)
  end #constructor 2 FlightData
end #struct FlightData


"""
# struct FlightDB

Database for aircraft data of different database types with fields:
- `inventory::Vector{FlightData}`
- `archive::Vector{FlightData}`
- `onlineData::Vector{FlightData}`
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
Immutable struct `DBMetadata`.


# Instantiation

Instatiate by giving a String with identifiers of the `DBtype` and an equal number
of `folder` paths as characters in the `DBtype` `String`. Optionally add a minimum
altitude threshold for the data (default = `15000`) and any remarks
(comments or additional data). Define the delimiter in the input files of the
online data with the keyword `odelim`. Use any character or string as delimiter.
By default (`odelim=nothing`), auto-detection is used.

    FlightDB(DBtype::String, folder::Union{String, Vector{String}}...;
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
`remarks` can be attached to `FlightDB`.

Alternatively, instatiate directly with the fields of `FlightDB`, where the correct
database type is checked, and wrong datasets are removed in every field.
"""
struct FlightDB
  inventory::Vector{FlightData}
  archive::Vector{FlightData}
  onlineData::Vector{FlightData}
  metadata::DBMetadata

  """
  Unmodified constructor for `FlightDB` with basic checks for correct dataset type
  in each dataset field.
  """
  function FlightDB(inventory::Vector{FlightData}, archive::Vector{FlightData},
    onlineData::Vector{FlightData}, metadata::DBMetadata)

    inventory = checkDBtype(inventory, "VOLPE AEDT")
    archive = checkDBtype(archive, "FlightAware")
    onlineData = checkDBtype(onlineData, "flightaware.com")

    new(inventory, archive, onlineData, metadata)
  end #constructor 1 FlightDB

  """
  Modified constructor creating the database from an identifer of the
  database type and the respective folder path for that database.
  """
  function FlightDB(DBtype::String, folder::String...; Float::DataType = Float32,
    altmin::Real=5_000, remarks=nothing, odelim::Union{Nothing,Char,String}=nothing)

    # Save time of database creation
    tstart = Dates.now()
    # Check DBtype addresses all folder paths
    if length(DBtype) ≠ length(folder)
      throw(ArgumentError("Number of characters in `DBtype` must match length of vararg `folder`"))
    end
    # Find database types
    i1 = [findall(isequal('i'), DBtype); findall(isequal('1'), DBtype)]
    i2 = [findall(isequal('a'), DBtype); findall(isequal('2'), DBtype)]
    i3 = [findall(isequal('o'), DBtype); findall(isequal('3'), DBtype)]

    # Load databases for each type
    # VOLPE AEDT inventory
    ifiles = String[]
    for i in i1
      findfiles!(ifiles, folder[i], ".csv")
    end
    inventory = loadInventory(ifiles...; Float=Float, altmin=altmin)
    # FlightAware commercial archive
    ifiles = String[]
    for i in i2
      findfiles!(ifiles, folder[i], ".csv")
    end
    archive = loadArchive(ifiles...; Float=Float, altmin=altmin)
    ifiles = String[]
    for i in i3
      findfiles!(ifiles, folder[i], ".tsv", ".txt", ".dat")
    end
    onlineData = loadOnlineData(ifiles...; Float=Float, altmin=altmin, delim=odelim)
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

    @info string("FlightDB data loaded in ",
      "$(join(loadtime.periods[1:min(2,length(loadtime.periods))], ", ")) to",
      "\n▪ inventory ($(length(inventory)) entries)\n▪ archive ($(length(archive)) entries)\n",
      "▪ onlineData ($(length(onlineData)) entries)\n▪ metadata")

    new(inventory, archive, onlineData,
      DBMetadata(altmin, (start=tmin, stop=tmax), tc, loadtime, remarks))
  end # constructor 2 FlightDB
end #struct FlightDB


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

  function SatData(data::DataFrame, metadata::SatMetadata)
    standardnames = ["time", "lat", "lon", "fileindex"]
    standardtypes = [Vector{DateTime}, Vector{<:AbstractFloat},
      Vector{<:AbstractFloat}, Vector{Int}]
    bounds = (:time => (DateTime(2006), Dates.now()), :lat => (-90,90),
      :lon => (-180,180), :fileindex => (1, Inf))
    checkcols!(data, standardnames, standardtypes, bounds, "CLay")
    new(data, metadata)
  end #constructor 1 SatData

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


"""
# struct CLay

CALIOP cloud layer `data` stored in a `DataFrame` with columns:
- `time::Vector{DateTime}` (time index)
- `lat::Vector{AbstractFloat}` (latitude position of current time index)
- `lon::Vector{AbstractFloat}` (longitude position of current time index)
- `layer::Vector{NamedTuple{(:top,:base),Tuple{Vector{<:AbstractFloat},Vector{<:AbstractFloat}}}}`
  (layer top/base heights in meters)
- `feature::Vector{Vector{Symbol}}` (symbol features of a layer)
- `OD::Vector{<:Vector{<:AbstractFloat}}` (layer optical depth)
- `IWP::Vector{<:Vector{<:Union{Missing,<:AbstractFloat}}}` (layer ice water path)
- `Ttop::Vector{<:Vector{<:AbstractFloat}}` (layer top temperature)
- `Htropo::Vector{<:AbstractFloat}` (tropopause height at current time index)
- `night::BitVector` (flag for nights (`true`))
- `averaging::Vector{<:Vector{Int}}` (horizontal averaging in km)

# Instantiation

    CLay(ms::mat.MSession, files::Vector{String},
      lidarrange::Tuple{Real,Real}=(15_000,-Inf), altmin::Real=5000, Float::DataType=Float32)

Construct `CLay` from a list of file names (including directories) and a running
MATLAB session `ms` and save data, if layers are within the bounds
of `lidarrange` and above flight `altmin` threshold. By default, all values are
saved as `Float32`, but can be set to any other precision by the `Float` kwarg.

Or construct `CLay` by directly handing over the `DataFrame` where the names, order,
and types of each columns are checked and attempted to correct:

    CLay(data::DataFrame) -> struct CLay
"""
struct CLay
  data::DataFrame

  """ Unmodified constructor for `CLay` """
  function CLay(data::DataFrame)
    standardnames = ["time", "lat", "lon",
      "layer", "feature", "OD", "IWP", "Ttop", "Htropo", "night", "averaging"]
      standardtypes = [Vector{DateTime}, Vector{<:AbstractFloat}, Vector{<:AbstractFloat},
      Vector{NamedTuple{(:top,:base),Tuple{T,T}}} where T<:Vector{<:AbstractFloat},
      Vector{Vector{Symbol}}, Vector{<:Vector{<:AbstractFloat}},
      Vector{<:Vector{<:Union{Missing,<:AbstractFloat}}}, Vector{<:Vector{<:AbstractFloat}},
      Vector{<:AbstractFloat}, BitVector, Vector{<:Vector{Int}}]
    bounds = (:lat => (-90,90), :lon => (-180,180))
    checkcols!(data, standardnames, standardtypes, bounds, "CLay")
    new(data)
  end #constructor 1 CLay

  """
  Modified constructor of `CLay` reading data from hdf `files` using MATLAB session `ms`
  in the `lidarrange` (top to bottom), if data is above `altmin`.
  """
  function CLay(ms::mat.MSession, files::Vector{String},
    lidarrange::Tuple{Real,Real}=(15_000,-Inf), altmin::Real=5000, Float::DataType=Float32)
    # Return default empty struct if files are empty
    isempty(files) && return CLay(Float)
    # Initialise arrays
    # essential data
    utc = Vector{Vector{DateTime}}(undef, length(files))
    lat = Vector{Vector{Float}}(undef, length(files))
    lon = Vector{Vector{Float}}(undef, length(files))
    # non-essential data
    layers = Vector{Vector{NamedTuple{(:top,:base),Tuple{Vector{Float},Vector{Float}}}}}(undef, length(files))
    feature = Vector{Vector{Vector{Symbol}}}(undef,length(files))
    OD = Vector{Vector{Vector{Float}}}(undef,length(files))
    IWP = Vector{Vector{Vector{Union{Missing,Float}}}}(undef,length(files))
    Ttop = Vector{Vector{Vector{Float}}}(undef,length(files))
    Htropo = Vector{Vector{Float}}(undef, length(files))
    night = Vector{BitVector}(undef, length(files))
    averaging = Vector{Vector{Vector{Int}}}(undef,length(files))
    # Convert mininmum flight altitude to meters
    altmin = ft2km(Float(altmin))
    # Loop over files
    for (i, file) in enumerate(files)
      ## Retrieve cloud layer data; assumes faulty files are filtered by SatData
      # Extract time
      mat.put_variable(ms, :file, file)
      mat.eval_string(ms, "clear t\ntry\nt = hdfread(file, 'Profile_UTC_Time');\nend")
      utc[i] = convertUTC.(mat.jarray(mat.get_mvariable(ms, :t))[:,2])
      # Extract lat/lon
      mat.eval_string(ms, "clear longitude\ntry\nlongitude = hdfread(file, 'Longitude');\nend")
      lon[i] = mat.jarray(mat.get_mvariable(ms, :longitude))[:,2]
      mat.eval_string(ms, "clear latitude\ntry\nlatitude = hdfread(file, 'Latitude');\nend")
      lat[i] = mat.jarray(mat.get_mvariable(ms, :latitude))[:,2]
      # Save time converted to UTC and lat/lon
      # utc[i], lon[i], lat[i] = convertUTC.(t), longitude, latitude

      ## Extract layer top/base, layer features and optical depth from hdf files
      mat.eval_string(ms, "clear basealt\ntry\nbasealt = hdfread(file, 'Layer_Base_Altitude');\nend")
      mat.eval_string(ms, "clear topalt\ntry\ntopalt = hdfread(file, 'Layer_Top_Altitude');\nend")
      basealt = mat.jarray(mat.get_mvariable(ms, :basealt))
      topalt = mat.jarray(mat.get_mvariable(ms, :topalt))
      mat.eval_string(ms, "clear FCF\ntry\nFCF = hdfread(file, 'Feature_Classification_Flags');\nend")
      FCF = mat.jarray(mat.get_mvariable(ms, :FCF))
      mat.eval_string(ms, "clear FOD\ntry\nFOD = hdfread(file, 'Feature_Optical_Depth_532');\nend")
      FOD = mat.jarray(mat.get_mvariable(ms, :FOD))
      mat.eval_string(ms, "clear IWPath\ntry\nIWPath = hdfread(file, 'Ice_Water_Path');\nend")
      IWPath = mat.jarray(mat.get_mvariable(ms, :IWPath))
      mat.eval_string(ms, "clear LTT\ntry\nLTT = hdfread(file, 'Layer_Top_Temperature');\nend")
      LTT = mat.jarray(mat.get_mvariable(ms, :LTT))
      mat.eval_string(ms, "clear Htropo\ntry\nHtropo = hdfread(file, 'Tropopause_Height');\nend")
      Htropo[i] = 1000vec(mat.jarray(mat.get_mvariable(ms, :Htropo)))
      mat.eval_string(ms, "clear daynight\ntry\ndaynight = hdfread(file, 'Day_Night_Flag');\nend")
      night[i] = Bool.(vec(mat.jarray(mat.get_mvariable(ms, :daynight))))
      mat.eval_string(ms, "clear average\ntry\naverage = hdfread(file, 'Horizontal_Averaging');\nend")
      horav = mat.jarray(mat.get_mvariable(ms, :average))
      # Loop over data and convert to TrackMatcher format
      layer = Vector{NamedTuple{(:top,:base),Tuple{Vector{Float},Vector{Float}}}}(undef,length(utc[i]))
      feat = Vector{Vector{Symbol}}(undef,length(utc[i]))
      optdepth = Vector{Vector{Float}}(undef,length(utc[i]))
      icewater = Vector{Vector{Union{Missing,Float}}}(undef,length(utc[i]))
      toptemp = Vector{Vector{Float}}(undef,length(utc[i]))
      average = Vector{Vector{Int}}(undef,length(utc[i]))
      for n = 1:length(utc[i])
        l = findall((basealt[n,:] .> 0) .& (topalt[n,:] .> 0) .& (basealt[n,:] .< lidarrange[1]) .&
          (topalt[n,:] .> lidarrange[2]) .& (topalt[n,:] .> altmin))
        layer[n], feat[n], optdepth[n], toptemp[n], icewater[n], average[n] = if isempty(l)
          (top = Float[], base = Float[]), Symbol[],
          Float[], Float[], Float[], Int[]
        else
          l = findall((basealt[n,:] .> 0) .& (topalt[n,:] .> 0) .& (basealt[n,:] .< lidarrange[1]) .&
            (topalt[n,:] .> lidarrange[2]))
          (top = [1000topalt[n, m] for m in l] , base = [1000basealt[n, m] for m in l]),
          [feature_classification(classification(FCF[n,m])...) for m in l],
          [FOD[n,m] for m in l],
          [LTT[n,m] for m in l],
          [IWPath[n,m] == -9999 ? missing : IWPath[n,m] for m in l],
          [1000horav[n,m] for m in l]
        end
      end # loop over time steps in current file
      layers[i], feature[i], OD[i], IWP[i], Ttop[i], averaging[i] =
        layer, feat, optdepth, icewater, toptemp, average
    end #loop over files

    # Construct and standardise data
    data = DataFrame(time=[utc...;], lat=[lat...;], lon=[lon...;],
      layer=[layers...;], feature=[feature...;], OD=[OD...;],
      IWP=[IWP...;], Ttop=[Ttop...;], Htropo = [Htropo...;],
      night = [night...;], averaging = [averaging...;])
    # Save time, lat/lon arrays in CLay struct
    new(data)
  end #constructor 2 CLay
end #struct CLay


""" External constructor for emtpy CLay struct """
CLay(Float::DataType=Float32) = CLay(DataFrame(time = DateTime[], lat = Float[], lon = Float[],
  layer = NamedTuple{(:top,:base),Tuple{Vector{Float},Vector{Float}}}[],
  feature = Vector{Symbol}[], OD = Vector{Float}[], IWP = Vector{Float}[],
  Ttop = Vector{Float}[], Htropo = Float[], night = BitVector(),
  averaging = Vector{Int}[]))


"""
# struct CPro

CALIOP cloud profile `data` stored in a `DataFrame` with columns:
- `time::Vector{DateTime}` (current time index)
- `lat::Vector{AbstractFloat}` (latitude coordinate for current time index)
- `lon::Vector{AbstractFloat}` (lonitude coordinate for current time index)
- `feature::Vector{<:Vector{<:Union{Missing,Symbol}}}`
  (feature symbols for every height level at current time index)
- `EC532::Vector{<:Vector{<:Union{Missing,AbstractFloat}}}`
  (extinction coefficient at 532nm at every height level in current time index)

# Instantiation

    CPro(ms::mat.MSession, files::Vector{String}, sattime::Vector{DateTime}, lidarprofile::NamedTuple)
      -> struct CPro

Construct `CPro` from a list of file names (including directories) and a running
MATLAB session `ms`. CPro data is only stored in the vicinity of intersections for
the designated `sattime`. Column data is stored height-resolved as defined by the
`lidarprofile`. By default, all values are stored as `Float32`, but can be set to
any other precision by the `Float` kwarg.

Or construct `CPro` by directly handing over the `DataFrame` where the names, order,
and types of each columns are checked and attempted to correct:

    CPro(data::DataFrame) -> struct CPro
"""
struct CPro
  data::DataFrame

  """ unmodified constructor """
  function CPro(data::DataFrame)
    standardnames = ["time", "lat", "lon", "feature", "EC532", "Htropo", "temp",
      "pressure", "rH", "IWC", "deltap", "CADscore", "night"]
    standardtypes = [Vector{DateTime}, Vector{<:AbstractFloat}, Vector{<:AbstractFloat},
      Vector{<:Vector{<:Union{Missing,Symbol}}}, Vector{<:Vector{<:Union{Missing,<:AbstractFloat}}},
      Vector{<:AbstractFloat}, Vector{<:Vector{<:Union{Missing,<:AbstractFloat}}},
      Vector{<:Vector{<:Union{Missing,<:AbstractFloat}}}, Vector{<:Vector{<:Union{Missing,<:AbstractFloat}}},
      Vector{<:Vector{<:Union{Missing,<:AbstractFloat}}}, Vector{<:Vector{<:Union{Missing,<:AbstractFloat}}},
      Vector{<:Vector{<:Union{Missing,Int8}}}, BitVector]
    bounds = (:lat => (-90,90), :lon => (-180,180), :Htropo => (4000,22_000),
      :temp => (-120,60), :pressure => (1,1086), :rH => (0,1.5), :IWC => (0,0.54),
      :deltap => (0,1), :CADscore => (-101,106))
    checkcols!(data, standardnames, standardtypes, bounds, "CPro")
    new(data)
  end #constructor 1 CPro

  """
  Modified constructor of `CPro` reading data from hdf `files` for all given `sattime` indices
  using MATLAB session `ms` and `lidarprofile` data, if data is above `altmin`.
  """
  function CPro(ms::mat.MSession, files::Vector{String}, sattime::Vector{DateTime},
    lidarprofile::NamedTuple, Float::DataType=Float32)
    # Return default empty struct if files are empty
    isempty(files) && return CPro(Float)
    # Initialise arrays
    # essential data
    utc = Vector{Vector{DateTime}}(undef, length(files))
    lat = Vector{Vector{Float}}(undef, length(files))
    lon = Vector{Vector{Float}}(undef, length(files))
    fcf = Vector{Vector{Vector{<:Union{Missing,UInt16}}}}(undef, length(files))
    # non-essential data
    ec532 = Vector{Vector{Vector{<:Union{Missing,Float}}}}(undef, length(files))
    Htropo = Vector{Vector{Float}}(undef, length(files))
    temp = Vector{Vector{Vector{<:Union{Missing,Float}}}}(undef, length(files))
    pres = Vector{Vector{Vector{<:Union{Missing,Float}}}}(undef, length(files))
    rH = Vector{Vector{Vector{<:Union{Missing,Float}}}}(undef, length(files))
    iwc = Vector{Vector{Vector{<:Union{Missing,Float}}}}(undef, length(files))
    deltap = Vector{Vector{Vector{<:Union{Missing,Float}}}}(undef, length(files))
    cad = Vector{Vector{Vector{<:Union{Missing,Int8}}}}(undef, length(files))
    night = Vector{BitVector}(undef, length(files))
    # Loop over files with cloud profile data
    for (i, file) in enumerate(files)
      ## Retrieve cloud profile data; assumes faulty files are filtered by SatData
      # Extract time
      mat.put_variable(ms, :file, file)
      mat.eval_string(ms, "clear t\ntry\nt = hdfread(file, 'Profile_UTC_Time');\nend")
      utc[i] = convertUTC.(mat.jarray(mat.get_mvariable(ms, :t))[:,2])
      # Extract lat/lon
      mat.eval_string(ms, "clear longitude\ntry\nlongitude = hdfread(file, 'Longitude');\nend")
      lon[i] = mat.jarray(mat.get_mvariable(ms, :longitude))[:,2]
      mat.eval_string(ms, "clear latitude\ntry\nlatitude = hdfread(file, 'Latitude');\nend")
      lat[i] = mat.jarray(mat.get_mvariable(ms, :latitude))[:,2]
      fcf[i] = get_lidarcolumn(UInt16, ms, "Atmospheric_Volume_Description", lidarprofile,
        coarse=false)
      # Extract non-essential data
      ec532[i] = get_lidarcolumn(Float, ms, "Extinction_Coefficient_532", lidarprofile,
        missingvalues = -9999)
      mat.eval_string(ms, "clear Htropo\ntry\nHtropo = hdfread(file, 'Tropopause_Height');\nend")
      Htropo[i] = 1000vec(mat.jarray(mat.get_mvariable(ms, :Htropo)))
      temp[i] = get_lidarcolumn(Float, ms, "Temperature", lidarprofile, missingvalues = -9999)
      pres[i] = get_lidarcolumn(Float, ms, "Pressure", lidarprofile, missingvalues = -9999)
      rH[i] = get_lidarcolumn(Float, ms, "Relative_Humidity", lidarprofile, missingvalues = -9999)
      iwc[i] = get_lidarcolumn(Float, ms, "Ice_Water_Content_Profile", lidarprofile,
        missingvalues = -9999)
      deltap[i] = get_lidarcolumn(Float, ms, "Particulate_Depalarization_Ratio_Profile_532",
        lidarprofile, missingvalues = -9999)
      cad[i] = get_lidarcolumn(Int8, ms, "CAD_Score", lidarprofile, coarse=false,
        missingvalues = -127)
      mat.eval_string(ms, "clear daynight\ntry\ndaynight = hdfread(file, 'Day_Night_Flag');\nend")
      night[i] = Bool.(vec(mat.jarray(mat.get_mvariable(ms, :daynight))))
    end #loop over files

    # Rearrange time vector and get time range
    utc = [utc...;]
    idx = [findfirst(utc .== t) for t in sattime]
    # Rearrange FCF vector and convert to symbols
    fcf = [fcf...;]
    avd =  Vector{Vector{Union{Missing,Symbol}}}(undef, length(fcf))
    for i = 1:length(fcf)
      vect = Vector{Union{Missing,Symbol}}(undef, length(fcf[i]))
      for j = 1:length(fcf[i])
        vect[j] = ismissing(fcf[i][j]) ? missing :
          feature_classification(classification(fcf[i][j])...)
      end
      avd[i] = vect
    end
    # Construct and standardise data
    data = DataFrame(time=utc[idx], lat=[lat...;][idx], lon=[lon...;][idx],
      feature=avd[idx], EC532=[ec532...;][idx], Htropo = [Htropo...;][idx],
      temp=[temp...;][idx], pressure = [pres...;][idx], rH = [rH...;][idx],
      IWC = [iwc...;][idx], deltap = [deltap...;][idx],
      CADscore = [cad...;][idx], night = [night...;][idx])
    # Save time, lat/lon arrays, and feature classification flags (FCF) in CPro struct
    new(data)
  end #constructor 2 CPro
end #struct CPro


""" External constructor for emtpy CPro struct """
CPro(Float::DataType=Float32) = CPro(DataFrame(time = DateTime[], lat = Float[], lon = Float[],
  feature = Vector{Symbol}[], EC532 = Vector{Float}[], Htropo = Float[], temp = Vector{Float}[],
  pressure = Vector{Float}[], rH = Vector{Float}[], IWC = Vector{Float}[],
  deltap = Vector{Float}[], CADscore = Vector{Int8}[], night = BitVector()))

## Define structs related to intersection data

"""
# struct Intersection

Intersection-related data with fields
- `data::DataFrame`
- `tracked::DataFrame`
- `accuracy::DataFrame`
- `metadata::XMetadata`


## data

Data related to spatial and temporal coordinates of intersections between satellite
and flight tracks and the meteorological conditions at the intersections.

`DataFrame` columns are:
- `id::Vector{String}`: unique intersection identifier
- `lat::Vector{<:AbstractFloat}`: latitude of intersection
- `lon::Vector{<:AbstractFloat}`: longitude of intersection
- `tdiff::Vector{Dates.CompoundPeriod}`: time difference between flight and satellite overpass
- `tflight::Vector{DateTime}`: time of aircraft at intersection
- `tsat::Vector{DateTime}`: time of satellite at intersection
- `feature::Vector{<:Union{Missing,Symbol}}`: atmospheric conditions at intersection


## tracked

Original track data in the vicinity of the intersection.

`DataFrame` columns are:
- `id: Vector{String}`: unique intersection identifier
- `flight::Vector{FlightData}`: `FlightData` in the vicinity of the intersection
- `CPro::Vector{CPro}`: `CPro` CALIOP profile data in the vicinity of the intersection
- `CLay::Vector{CLay}`: `CLay` CALIOP layer data in the vicinity of the intersection


## accuracy

Measures about the accuracy of the intersection calculations and the quality of the track data.

`DataFrame` columns are:
- `id: Vector{String}`: unique intersection identifier
- `intersection::Vector{<:AbstractFloat}`: accuracy of the intersection calculation in meters
- `flightcoord::Vector{<:AbstractFloat}`: distance of nearest tracked flight data
  to calculated intersection in meters
- `satcoord::Vector{<:AbstractFloat}`: distance of nearest tracked sat data
  to calculated intersection in meters
- `flighttime::Vector{Dates.CompoundPeriod}`: time difference between measurement
  of nearest tracked flight data and calculated time of aircraft at intersection
- `sattime::Vector{Dates.CompoundPeriod}`: time difference between measurement
  of nearest tracked sat data and calculated time of satellite at intersection


# Instantiation

    Intersection(
      flights::FlightDB,
      sat::SatData,
      savesecondsattype::Bool=false;
      maxtimediff::Int=30,
      flightspan::Int=0,
      satspan::Int=15,
      lidarrange::Tuple{Real,Real}=(15_000,-Inf),
      stepwidth::Real=1000,
      Xradius::Real=20_000,
      Float::DataType=Float32,
      remarks=nothing
    ) -> struct Intersection

Construct `Intersection` from the preloaded `FlightDB` and `SatData` with the option
to save the other satellite data type not used in `sat` (either `CLay` or `CPro`),
when `savesecondsattype` is set to `true`. Folder structure and file names must
be identical only with `CLay`/`CPro` interchanged for this option to work.

The following parameters can be set to influence intersection calculations or
data saved to the struct:
- `maxtimediff::Int=30`: maximum time difference allowed between aircraft passage
  and satellite overpass at intersection
- `flightspan::Int=0`: Number of additional data points of original track data
  saved in the vicinity of the intersection and stored in `Intersection.tracked.flight`
- `satspan::Int=15`: Number of additional data points of original track data
  saved in the vicinity of the intersection and stored in `Intersection.tracked.CPro`
  and `Intersection.tracked.CLay`
- `lidarrange::Tuple{Real,Real}=(15,-Inf)`: lidar measurements saved for column heights
  between `(max, min)` (set to `Inf`/`-Inf` to store all values up to top/bottom)
- `stepwidth::Real=1000`: step width of interpolation in flight and sat tracks
  in meters (partially internally converted to degrees at equator)
- `Xradius::Real=20_000`: radius in meters within which multiple finds of an
  intersection are disregarded and only the most accurate is counted
- `Float::DataType=Float32`: Set the precision of floating point numbers
  (single precision by default)
- `remarks=nothing`: any data or remarks attached to the metadata

Or construct `Intersection` by directly handing over the different `DataFrame`s where the names, order,
and types of each columns are checked and attempted to correct, together with the metadata:

    function Intersection(
      data::DataFrame,
      tracked::DataFrame,
      accuracy::DataFrame,
      metadata::XMetadata
    ) -> struct Intersection
"""
struct Intersection
  data::DataFrame
  tracked::DataFrame
  accuracy::DataFrame
  metadata::XMetadata


  """ Unmodified constructor for `Intersection` """
  function Intersection(
    data::DataFrame,
    tracked::DataFrame,
    accuracy::DataFrame,
    metadata::XMetadata
  )
    # Check data
    standardnames = ["id", "lat", "lon", "alt", "tdiff", "tflight", "tsat", "feature"]
    standardtypes = [Vector{String}, Vector{<:AbstractFloat}, Vector{<:AbstractFloat},
      Vector{<:AbstractFloat}, Vector{Dates.CompoundPeriod}, Vector{DateTime},
      Vector{DateTime}, Vector{<:Union{Missing,Symbol}}]
    bounds = (:lat => (-90,90), :lon => (-180,180), :alt => (0, Inf))
    checkcols!(data, standardnames, standardtypes, bounds, "Intersection.data")
    # Check tracked (measured data)
    standardnames = ["id", "flight", "CPro", "CLay"]
    standardtypes = [Vector{String}, Vector{FlightData}, Vector{CPro}, Vector{CLay}]
    bounds = ()
    checkcols!(tracked, standardnames, standardtypes, bounds, "Intersection.tracked",
      essentialcols = [1])
    # Check accuracy
    standardnames = ["id", "intersection", "flightcoord", "satcoord", "flighttime", "sattime"]
    standardtypes = [Vector{String}, Vector{<:AbstractFloat}, Vector{<:AbstractFloat},
      Vector{<:AbstractFloat}, Vector{Dates.CompoundPeriod}, Vector{Dates.CompoundPeriod}]
    bounds = ()
    checkcols!(accuracy, standardnames, standardtypes, bounds, "Intersection.accuracy",
      essentialcols = [1])
    new(data, tracked, accuracy, metadata)
  end #constructor 1 Intersection


  """ Modified constructor with some automated calculations of the intersection data. """
  function Intersection(
    flights::FlightDB,
    sat::SatData,
    savesecondsattype::Bool=false;
    maxtimediff::Int=30,
    flightspan::Int=0,
    satspan::Int=15,
    lidarrange::Tuple{Real,Real}=(15_000,-Inf),
    stepwidth::Real=1000,
    Xradius::Real=20_000,
    epsilon::Real=-1,
    tolerance::Real=-1,
    Float::DataType=Float32,
    remarks=nothing
  )
    # Initialise DataFrames with Intersection data and monitor start time
    tstart = Dates.now()
    Xdata = DataFrame(id=String[], lat=Float[], lon=Float[],
      alt=Float[], tdiff=Dates.CompoundPeriod[], tflight = DateTime[],
      tsat = DateTime[], feature = Union{Missing,Symbol}[])
    track = DataFrame(id=String[], flight=FlightData[], CPro=CPro[], CLay=CLay[])
    accuracy = DataFrame(id=String[], intersection=Float[], flightcoord=Float[],
      satcoord=Float[], flighttime=Dates.CompoundPeriod[], sattime=Dates.CompoundPeriod[])
    # Get lidar altitude levels
    lidarprofile = get_lidarheights(lidarrange, Float)
    # Save stepwidth in degrees at equator using Earth's equatorial circumference to convert
    degsteps  = stepwidth*360/40_075_017
    # Calculate default tolerances
    epsilon == -1 && (epsilon = 2stepwidth)
    tolerance == -1 && (tolerance = stepwidth)
    # New MATLAB session
    ms = mat.MSession()
    # Loop over data from different datasets and interpolate track data and time, throw error on failure
    flightdata = [[getfield(flights, f) for f in fieldnames(FlightDB)[1:end-1]]...;]
    prog = pm.Progress(length(flightdata), "find intersections...")
    for flight in flightdata
      try
        # Find sat tracks in the vicinity of flight tracks, where intersections are possible
        overlap = findoverlap(flight, sat, maxtimediff)
        if isempty(overlap)
          pm.next!(prog, showvalues = [(:hits, length(Xdata.id)),
            (:featured, length(Xdata.id[.!ismissing.(Xdata.feature) .&
            (Xdata.feature .≠ :no_signal) .& (Xdata.feature .≠ :clear)]))])
          continue
        end
        # Interpolate trajectories with PCHIP method
        sattracks = interpolate_satdata(sat, overlap)
        flighttracks = interpolate_flightdata(flight, degsteps)
        # Calculate intersections and store data and metadata in DataFrames
        currdata, currtrack, curraccuracy = find_intersections(ms, flight,
          flighttracks, flights.metadata.altmin, sat, sattracks, maxtimediff,
          Xradius, epsilon, tolerance, lidarprofile, lidarrange,
          flightspan, satspan, savesecondsattype,Float)
        append!(Xdata, currdata); append!(track, currtrack)
        append!(accuracy, curraccuracy)
      catch err
        @debug begin
          @show flight.metadata.dbID
          rethrow(err)
        end
        # Issue warning on failure of interpolating track or time data
        @warn("Track data and/or time could not be interpolated. Data ignored.",
          dataset = flight.metadata.source, flight = flight.metadata.dbID)
      end
      # Monitor progress for progress bar
      pm.next!(prog, showvalues = [(:hits, length(Xdata.id)),
        (:featured, length(Xdata.id[.!ismissing.(Xdata.feature) .&
        (Xdata.feature .≠ :no_signal) .& (Xdata.feature .≠ :clear)]))])
    end #loop over flights
    pm.finish!(prog)
    # Close MATLAB session after looping over all data
    mat.close(ms)
    # Calculate load time
    tend = Dates.now()
    tc = tz.ZonedDateTime(tend, tz.localzone())
    loadtime = Dates.canonicalize(Dates.CompoundPeriod(tend - tstart))
    # Return Intersections after completion
    @info string("Intersection data ($(length(Xdata[!,1])) matches) loaded in ",
      "$(join(loadtime.periods[1:min(2,length(loadtime.periods))], ", ")) to",
      "\n▪ data\n▪ tracked\n▪ accuracy\n▪ metadata")
    new(Xdata, track, accuracy, XMetadata(maxtimediff, stepwidth, Xradius,
      lidarrange, lidarprofile, sat.metadata.type, sat.metadata.date,
      flights.metadata.altmin, flights.metadata.date, tc, loadtime, remarks))
  end #constructor Intersection
end #struct Intersection


## Export structs
export FlightDB, FlightData, SatData, CLay, CPro, Intersection,
       FlightMetadata, SatMetadata, DBMetadata, XMetadata


## Import functions for Julia include files
include("auxiliary.jl")       # helper functions
include("lidar.jl")           # functions related to processing CALIOP lidar data
include("loadFlightData.jl")  # functions related to loading flight databases/datasets
include("match.jl")           # functions related to finding track intersections

end # module TrackMatcher
