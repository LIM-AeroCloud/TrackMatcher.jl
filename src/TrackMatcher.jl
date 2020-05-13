"""
# Module TrackMatcher

To find intersection between different trajectories. The module is aimed to find
intersections between aircraft and satellite tracks, but can be used for ship
or cloud tracks as well.

## Exported structs

- `FlightDB` stores flight track data and other relevant aircraft related data
  from 3 different inventories:
  - `inventory`: VOLPE AEDT inventory
  - `archive`: commercially available database by FlightAware
  - `onlineData`: free online data by FlightAware
- `FlightData` stores data of a single flight in `FlightDB`
- `FlightMetadata` holds metadata to every flight
- `SatDB` stores CALIPSO cloud layer and profile data from the CALIOP satellite
- `CLay` CALIPSO cloud layer data
- `CPro` CALIPSO cloud profile data
- `SatMetadata` with metadata of the CALIPSO cloud profile data
- `DBMetadata` with metadata of the databases for satellite and flight data
- `Intersection`: position of intersection between aircraft and satellite trajectory,
  together with time difference between crossing; and `FlightData` and `SatDB` in the
  vicinity of the intersection as well as information about the accuracy of the data
- `XMetadata` with metadata for the `Intersection` data
"""
module TrackMatcher

## Import Julia packages
import CSV
import DataFrames; const df = DataFrames
import Dates
import TimeZones; const tz = TimeZones
import Geodesy; const geo = Geodesy
import MATLAB; const mat = MATLAB
import Statistics; const stats = Statistics
import ProgressMeter; const pm = ProgressMeter
import Logging; const logg = Logging
# Import structs from packages
import DataFrames.DataFrame
import Dates: DateTime, Date, Time
import TimeZones.ZonedDateTime


# Define Logger with log level
logger = try logg.SimpleLogger(logfile, logg.Info)
catch; logg.ConsoleLogger(stdout, logg.Info)
end
logg.global_logger(logger)


## Define own Metadata structs
"""
# struct FlightMetadata

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

## dbID
Database ID – integer counter for `inventory` and FlightAware `onlineData`,
String with information about `FlightID`, `route`, and scheduled arrival for
FlightAware archived data.

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
a modified constructor and `dbID`, `flightID`, `aircraft` type, `route`, and `file`.
Fields `area` and `date` are calculated from `lat`/`lon`, and `date` vectors.

    function FlightMetadata(
      dbID::Union{Int,AbstractString},
      flightID::Union{Missing,AbstractString},
      route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
      aircraft::Union{Missing,AbstractString},
      date::Vector{DateTime},
      lat::Vector{<:Union{Missing,AbstractFloat}},
      lon::Vector{<:Union{Missing,AbstractFloat}},
      useLON::Bool,
      flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,AbstractFloat,AbstractFloat}}}},
      source::AbstractString,
      file::AbstractString
    ) -> struct FlightMetadata

Or construct `FlightMetadata` by directly handing over every field:

    function FlightMetadata(
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
struct FlightMetadata
  dbID::Union{Int,AbstractString}
  flightID::Union{Missing,AbstractString}
  route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}}
  aircraft::Union{Missing,AbstractString}
  date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,AbstractFloat}}
  flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,AbstractFloat,AbstractFloat}}}}
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
    area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,AbstractFloat}},
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,AbstractFloat,AbstractFloat}}}},
    useLON::Bool,
    source::AbstractString,
    file::AbstractString
  )

    new(dbID, flightID, route, aircraft, date, area, flex, useLON, source, file)
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
    lat::Vector{<:Union{Missing,AbstractFloat}},
    lon::Vector{<:Union{Missing,AbstractFloat}},
    useLON::Bool,
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,AbstractFloat,AbstractFloat}}}},
    source::AbstractString,
    file::AbstractString
  )

    elonmax = isempty(lon[lon.≥0]) ? NaN : maximum(lon[lon.≥0])
    elonmin = isempty(lon[lon.≥0]) ? NaN : minimum(lon[lon.≥0])
    wlonmax = isempty(lon[lon.<0]) ? NaN : maximum(lon[lon.<0])
    wlonmin = isempty(lon[lon.<0]) ? NaN : minimum(lon[lon.<0])
    area = (latmin=minimum(lat), latmax=maximum(lat),
      elonmin=elonmin, elonmax=elonmax, wlonmin=wlonmin, wlonmax=wlonmax)
    new(dbID, flightID, route, aircraft, (start=date[1], stop=date[end]), area,
      flex, useLON, source, file)
  end #constructor 2 FlightMetadata
end #struct FlightMetadata


"""


"""
struct SatMetadata
  files::Dict{Int,String}
  type::Symbol
  date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  created::Union{DateTime,ZonedDateTime}
  loadtime::Dates.CompoundPeriod
  remarks

  function SatMetadata(
    files::Vector{String},
    date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    loadtime::Dates.CompoundPeriod=Dates.canonicalize(Dates.CompoundPeriod());
    remarks=nothing
  )
    # Find type of satellite data based on first 50 files (~2 days)
    type = occursin("CLay", files[1]) ≥
      count(occursin.("CPro", files[1:min(length(files), 50)])) ? :Clay : :CPro
    # Create a new instance of SatMetadata
    new(Dict(enumerate(files)), type, date, tz.now(tz.localzone()), loadtime, remarks)
  end #constructor 2 SatMetadata
end #struct SatMetadata



"""
# struct CProMetadata

Immutable struct with metadata for CALIOP cloud layer profiles storing information
about the height levels of the measurements.
"""
struct CProMetadata
  lidarrange
  lidarlevels

    function CProMetadata(
      lidarrange::NamedTuple{(:top,:bottom), Tuple{Real,Real}},
      lidarlevels::NamedTuple{(:coarse,:fine,:itop,:ibottom,:i30),
        Tuple{T,T,Int,Int,Int}}
    ) where T <: Vector{<:AbstractFloat}
      new(lidarrange, lidarlevels)
    end

    function CProMetadata(
      lidarrange::Tuple{Real,Real},
      lidarlevels::NamedTuple{(:coarse,:fine,:itop,:ibottom,:i30),
        Tuple{T,T,Int,Int,Int}}
    ) where T <: Vector{<:AbstractFloat}
      new((top=lidarrange[1], bottom=lidarrange[2]), lidarlevels)
    end
end #struct SatMetadata


"""
# struct DBMetadata

Immutable struct with additional information of databases:

- `date`: NamedTuple with entries `start`/`stop` giving the time range of the database
- `created`: time of creation of database
- `loadtime`: time it took to read data files and load it to the struct
- `remarks`: any additional data or comments that can be attached to the database
"""
struct DBMetadata
  date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  created::Union{DateTime,ZonedDateTime}
  loadtime::Dates.CompoundPeriod
  remarks
end #struct DBMetadata


"""
# struct XMetadata

Immutable struct with additional information of intersection data:

- `maxtimediff`: maximum time difference allowed between satellite overpass and
  aircraft passing at intersection
- `stepwidth`: stepwidth (in degrees at equator) used to interpolate track data
- `Xradius`: radius in meters around an intersection in which further intersections
  will be removed as duplicates due to the interpolation algorithm
- `preferred`: Symbol, which satellite data type was preferred to find intersections in
  flight and satellite trajectories (either `:CLay` (default) or `:CPro`)
- `created`: time of creation of database
- `loadtime`: time it took to find intersections and load it to the struct
- `remarks`: any additional data or comments that can be attached to the database
"""
struct XMetadata
  maxtimediff::Int
  stepwidth::AbstractFloat
  Xradius::Real
  lidarrange::Tuple{Real,Real}
  lidarprofile::NamedTuple
  created::Union{DateTime,ZonedDateTime}
  loadtime::Dates.CompoundPeriod
  remarks
end #struct XMetaData


## Define structs related to flight data
"""
# struct FlightData

Aircraft data with fields
- `data::DataFrame`
- `metadata::FlightMetadata`

The `DataFrame` of `data` has columns in the following order with the respective types:

- `time::Vector{DateTime}`                        (time stamp of measurement)
- `lat::Vector{<:Union{Missing,AbstractFloat}}`   (latitude)
- `lon::Vector{<:Union{Missing,AbstractFloat}}`   (longitude)
- `alt::Vector{<:Union{Missing,AbstractFloat}}`   (altitude)
- `heading::Vector{<:Union{Missing,Int}}`         (heading/direction of the aircraft)
- `climb::Vector{<:Union{Missing,Int}}`           (climbing (positive)/sinking (negative values) rate)
- `speed::Vector{<:Union{Missing,AbstractFloat}}` (velocity)


# Instantiation

    FlightData(time::Vector{ZonedDateTime}, lat::Vector{<:Union{Missing,AbstractFloat}},
      lon::Vector{<:Union{Missing,AbstractFloat}}, alt::Vector{<:Union{Missing,AbstractFloat}},
      heading::Vector{<:Union{Missing,Int}}, climb::Vector{<:Union{Missing,Int}},
      speed::Vector{<:Union{Missing,AbstractFloat}}, dbID::Union{Int,AbstractString},
      flightID::Union{Missing,AbstractString}, aircraft::Union{Missing,AbstractString},
      route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
      file::AbstractString) -> struct FlightData

Construct `FlightData` from columns for the `data` `DataFrame` and additonal information
`dbID`, `flightID`, `aircraft` type, `route`, and `file` name for `FlightMetadata`.
For `data`, `time` is the only exception, which is given as `ZonedDateTime` and
will be converted to `UTC` standard time.

Or construct by directly handing over every field:

    FlightData(time::Vector{DateTime}, lat::Vector{<:Union{Missing,AbstractFloat}},
      lon::Vector{<:Union{Missing,AbstractFloat}}, alt::Vector{<:Union{Missing,AbstractFloat}},
      heading::Vector{<:Union{Missing,Int}}, climb::Vector{<:Union{Missing,Int}},
      speed::Vector{<:Union{Missing,AbstractFloat}}, metadata::FlightMetadata)

Checks exist that the order, names, and types of the `data` `DataFrame` are correct.
"""
struct FlightData
  data::DataFrame
  metadata::FlightMetadata

  """ Unmodified constructor for `FlightData` with basic checks for correct `data`"""
  function FlightData(data::DataFrame, metadata::FlightMetadata)

    # Column checks and warnings
    standardnames = [:time, :lat, :lon, :alt, :heading, :climb, :speed]
    standardtypes = [Union{DateTime,Vector{DateTime}},
      Union{AbstractFloat,Vector{<:AbstractFloat}},
      Union{AbstractFloat,Vector{<:AbstractFloat}},
      Union{Missing,AbstractFloat,Vector{<:Union{Missing,AbstractFloat}}},
      Union{Missing,Int,Vector{<:Union{Missing,Int}}},
      Union{Missing,Int,Vector{<:Union{Missing,Int}}},
      Union{Missing,AbstractFloat,Vector{<:Union{Missing,AbstractFloat}}}]
    bounds = [(0,Inf), (0, 360), (-Inf, Inf), (0, Inf)]
    data = checkcols(data, standardnames, standardtypes, bounds,
      metadata.source, metadata.dbID)
    new(data,metadata)
  end #constructor 1 FlightData

  """ Modified constructor with variable checks and some automated calculation of fields """
  function FlightData(
    flightdata::DataFrame,
    dbID::Union{Int,AbstractString},
    flightID::Union{Missing,AbstractString},
    aircraft::Union{Missing,AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,AbstractFloat,AbstractFloat}}}},
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
    metadata = FlightMetadata(dbID,flightID,route,aircraft,t,lat,lon,useLON,flex,source,file)

    # Instatiate new FlightData
    new(DataFrame(time=t,lat=lat,lon=lon,alt=alt,heading=heading,climb=climb,speed=speed),metadata)
  end #constructor 2 FlightData
end #struct FlightData


"""
# struct FlightDB

Database for aircraft data of different database types with fields:
- `inventory::Vector{FlightData}`
- `archive::Vector{FlightData}`
- `onlineData::Vector{FlightData}`
- `metadata`: `DBMetatadata` with information about
  - `date`: date range of the dataset
  - `created`: time of creation
  - `loadtime`
  - `remarks` (any additional data or comments)

## inventory
Flight data from csv files.

## archive
Commercial flight data by FlightAware.

## onlineData
Online data from FlightAware website.

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
      altmin::Int=15_000, remarks=nothing, odelim::Union{Nothing,Char,String}=nothing)

`DBtype` can be identified with:
- `1` or `i`: VOLPE AEDT inventory
- `2` or `a`: FlightAware archived data (commercially available)
- `3` or `o`: flightaware.com online data

Or instatiate directly with the fields of `FlightDB`, where the correct database
type is checked, and wrong datasets are removed in every field.
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
  function FlightDB(DBtype::String, folder::Union{String, Vector{String}}...;
    altmin::Int=15_000, remarks=nothing, odelim::Union{Nothing,Char,String}=nothing)

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
    inventory = loadInventory(ifiles, altmin=altmin)
    # FlightAware commercial archive
    ifiles = String[]
    for i in i2
      findfiles!(ifiles, folder[i], ".csv")
    end
    archive = loadArchive(ifiles, altmin=altmin)
    ifiles = String[]
    for i in i3
      if VERSION ≥ v"1.2"
        findfiles!(ifiles, folder[i], ".txt", ".dat")
      else
        findfiles!(ifiles, folder[i], ".txt")
        findfiles!(ifiles, folder[i], ".dat")
      end
    end
    onlineData = loadOnlineData(ifiles, altmin=altmin, delim=odelim)
    tmin = minimum([[f.metadata.date.start for f in inventory];
      [f.metadata.date.start for f in archive];
      [f.metadata.date.start for f in onlineData]])
    tmax = maximum([[f.metadata.date.stop for f in inventory];
      [f.metadata.date.stop for f in archive];
      [f.metadata.date.stop for f in onlineData]])

    tend = Dates.now()
    tc = tz.ZonedDateTime(tend, tz.localzone())
    loadtime = Dates.canonicalize(Dates.CompoundPeriod(tend - tstart))

    @info string("FlightDB data loaded to properties in ",
      "$(join(loadtime.periods[1:min(2,length(loadtime.periods))], ", "))",
      "\n▪ inventory\n▪ archive\n▪ onlineData\n▪ metadata")

    new(inventory, archive, onlineData,
      DBMetadata((start=tmin, stop=tmax), tc, loadtime, remarks))
  end # constructor 2 FlightDB
end #struct FlightDB


## Define structs related to sat data

struct SatData
  data::DataFrame
  metadata::SatMetadata

  function SatData(folders::String...; remarks=nothing)
    tstart = Dates.now()
    # Scan folders for HDF4 files
    files = String[];
    for folder in folders
      try findfiles!(files, folder, ".hdf")
      catch
        @warn "Read error in $folder.\nData skipped"
      end
    end
    # Create empty struct, if no data files were found
    isempty(files) && return new(DataFrame(time=DateTime[], lat=AbstractFloat[],
      lon=AbstractFloat[], fileindex=Int[]), SatMetadata(files, (start=tstart, stop=tstart),
      Dates.canonicalize(Dates.CompoundPeriod()), remarks=remarks))
    # Find type of satellite data based on first 50 files (~2 days)
    type = count(occursin.("CLay", files[1:min(length(files), 50)])) ≥
      count(occursin.("CPro", files[1:min(length(files), 50)])) ? "Clay" : "CPro"
    # Start MATLAB session
    ms = mat.MSession()
    # Initialise arrays
    utc = Vector{Vector{DateTime}}(undef, length(files))
    lat = Vector{Vector{AbstractFloat}}(undef, length(files))
    lon = Vector{Vector{AbstractFloat}}(undef, length(files))
    fileindex = Vector{Vector{Int}}(undef, length(files))
    # Loop over files
    prog = pm.Progress(length(files), "load sat data...")
    for (i, file) in enumerate(files)
      if !occursin(type, basename(file))
        @warn string("Wrong satellite type.\n",
          "Data must be of type $type with indication in the file name. Data skipped.")
        utc[i], lat[i], lon[i], fileindex[i] =
          DateTime[], AbstractFloat[], AbstractFloat[], Int[]
        continue
      end
      # Find files with cloud layer data
      try
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
        utc[i] = convertUTC.(t)
        lon[i] = longitude
        lat[i] = latitude
        fileindex[i] = [i for index in t]
      catch
        # Skip data on failure and warn
        @warn "read error in CALIPSO granule $(splitext(basename(file))[1])\ndata skipped"
        utc[i], lat[i], lon[i], fileindex[i] =
          DateTime[], AbstractFloat[], AbstractFloat[], Int[]
      end
      # Monitor progress for progress bar
      pm.next!(prog, showvalues = [(:date,Dates.Date(splitdir(dirname(file))[2], "y_m_d"))])
    end #loop over files
    pm.finish!(prog)

    # Close MATLAB session
    mat.close(ms)
    # Calculate time span of satellite data
    sattime = [DateTime[]; utc...]
    tmin = minimum(sattime)
    tmax = maximum(sattime)
    tend = Dates.now()
    # Save computing times
    loadtime = Dates.canonicalize(Dates.CompoundPeriod(tend - tstart))

    # Instantiate new struct
    @info string("SatDB data loaded to properties in ",
      "$(join(loadtime.periods[1:min(2,length(loadtime.periods))], ", "))",
      "\n▪ CLay\n▪ CPro\n▪ metadata")
    new(DataFrame(time=sattime, lat=[AbstractFloat[]; lat...],
      lon=[AbstractFloat[]; lon...], fileindex=[Int[]; fileindex...]),
      SatMetadata(files, (start=tmin, stop=tmax), loadtime, remarks=remarks))
  end #constructor 2 SatData
end #struct SatData





"""
# struct CLay

CALIOP cloud layer `data` stored in a `DataFrame` with columns:
- `time::Vector{DateTime}`
- `lat::Vector{AbstractFloat}`
- `lon::Vector{AbstractFloat}`

# Instantiation

    CLay(ms::mat.MSession, files::Vector{String}) -> struct CLay

Construct `CLay` from a list of file names (including directories) and a running
MATLAB session.

Or construct `CLay` by directly handing over the `DataFrame` where the names, order,
and types of each columns are checked and attempted to correct:

    CLay(data::DataFrame) -> struct CLay
"""
struct CLay
  data::DataFrame

  """ Unmodified constructor for `CLay` """
  function CLay(data::DataFrame)
    standardnames = [:time, :lat, :lon]
    standardtypes = [Vector{DateTime}, Vector{<:AbstractFloat}, Vector{<:AbstractFloat}]
    bounds = Tuple{Real,Real}[]
    data = checkcols(data, standardnames, standardtypes, bounds, "CLay", nothing)
    new(data)
  end #constructor 1 CLay

  """
  Modified constructor of `CLay` reading data from hdf files given in `folders...`
  using MATLAB session `ms`.
  """
  function CLay(ms::mat.MSession, files::Vector{String})
    # Initialise arrays
    utc = Vector{DateTime}[]; lon = Vector{AbstractFloat}[]; lat = Vector{AbstractFloat}[]
    # Loop over files
    prog = pm.Progress(length(files), "load CLay data...")
    for file in files
      # Find files with cloud layer data
      try
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
        push!(utc, convertUTC.(t))
        push!(lon, longitude)
        push!(lat, latitude)
      catch
        # Skip data on failure and warn
        @warn string("read error in CALIPSO granule ",
          "$(splitext(basename(file))[1])\ndata skipped")
      end
      # Monitor progress for progress bar
      pm.next!(prog, showvalues = [(:date,Dates.Date(splitdir(dirname(file))[2], "y_m_d"))])
    end #loop over files
    pm.finish!(prog)

    # Save time, lat/lon arrays in CLay struct
    new(DataFrame(time=[DateTime[]; utc...], lat=[AbstractFloat[]; lat...],
      lon=[AbstractFloat[]; lon...]))
  end #constructor 2 CLay
end #struct CLay


""" External constructor for emtpy CLay struct """
CLay() = CLay(DataFrame(time = DateTime[], lat = AbstractFloat[], lon = AbstractFloat[]))


"""
# struct CPro

CALIOP cloud profile `data` stored in a `DataFrame` with columns:
- `time::Vector{DateTime}`
- `lat::Vector{AbstractFloat}`
- `lon::Vector{AbstractFloat}`
- `FCF::Vector{UInt16}`
- `EC532::Vector{<:Union{Missing,Symbol}}`

Additional `SatMetadata` are stored in a immutable struct with fields:
- `lidarrange::NamedTuple{(:top,:bottom), Tuple{Real,Real}}`
- `lidarlevels::NamedTuple{(:coarse,:fine,:itop,:ibottom,:i30),Tuple{Vector{<:AbstractFloat},Vector{<:AbstractFloat},Int,Int,Int}}`

# Instantiation

    CPro(ms::mat.MSession, files::Vector{String}, lidar::NamedTuple, lidarrange::Tuple{Real,Real})
      -> struct CPro

Construct `CPro` from a list of file names (including directories) and a running
MATLAB session. Furthermore, provide a `NamedTuple` `lidar` with information about
the heights of the lidar levels and indices of important height levels in the
vectorized data.

Or construct `CPro` by directly handing over the `DataFrame` where the names, order,
and types of each columns are checked and attempted to correct:

    CPro(data::DataFrame, metadata::SatMetadata) -> struct CPro
"""
struct CPro
  data::DataFrame
  metadata::CProMetadata

  """ unmodified constructor """
  function CPro(data::DataFrame, metadata::CProMetadata)
    standardnames = [:time, :lat, :lon, :FCF, :EC532]
    standardtypes = [Vector{DateTime}, Vector{<:AbstractFloat}, Vector{<:AbstractFloat},
      Vector{<:Vector{<:Union{Missing,UInt16}}}, Vector{<:Vector{<:Union{Missing,AbstractFloat}}}]
    bounds = Tuple{Real,Real}[]
    data = checkcols(data, standardnames, standardtypes, bounds, "CPro", nothing)
    new(data, metadata)
  end #constructor 1

  function CPro(ms::mat.MSession, files::Vector{String}, sattime::Vector{DateTime},
    lidar::NamedTuple, lidarrange::Tuple{Real,Real})
    # Initialise arrays
    # essential data
    utc = Vector{Vector{DateTime}}(undef, length(files))
    lat = Vector{Vector{AbstractFloat}}(undef, length(files))
    lon = Vector{Vector{AbstractFloat}}(undef, length(files))    # non-essential data
    avd = Vector{Vector{Vector{<:Union{Missing,UInt16}}}}(undef, length(files))
    ec532 = Vector{Vector{Vector{<:Union{Missing,Float32}}}}(undef, length(files))
    # Loop over files
    prog = pm.Progress(length(files), "load CPro data...")
    # Loop over files with cloud profile data
    for (i, file) in enumerate(files)
      # Retrieve essential data
      try
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
        utc[i] = convertUTC.(t)
        lon[i] = longitude
        lat[i] = latitude
      catch err
        @debug rethrow(err)
        # Warn on failure
        @warn string("read error in CALIPSO granule ",
          "$(splitext(basename(file))[1])\ndata skipped")
        # Monitor progress and skip to next file
        pm.next!(prog, showvalues = [(:date,Dates.Date(splitdir(dirname(file))[2], "y_m_d"))])
        continue
      end
      # Retrieve non-essential data
      # Extract feature classification flags
      avd[i] = try get_lidarcolumn(avd, ms, i, "Atmospheric_Volume_Description", lidar)
      catch
        @debug rethrow(err)
        @warn "missing Atmosphericc Volume Description for granule in $(splitext(basename(file))[1])"
        t = mat.jarray(mat.get_mvariable(ms, :t))[:,2]
        [[missing for i = 1:length(lidar.fine)] for i in t]
      end
      ec532[i] = try get_lidarcolumn(ec532, ms, i, "Extinction_Coefficient_532", lidar,
        true, missingvalues = -9999)
      catch err
        rethrow(err)
        @warn "missing Extinction Coefficient 532nm for granule in $(splitext(basename(file))[1])"
        t = mat.jarray(mat.get_mvariable(ms, :t))[:,2]
        [[missing for i = 1:length(lidar.coarse)] for i in t]
      end
      # Monitor progress for progress bar
      pm.next!(prog, showvalues = [(:date,Dates.Date(splitdir(dirname(file))[2], "y_m_d"))])
    end #loop over files
    pm.finish!(prog)

    utc = [DateTime[]; utc...]
    idx = [findfirst(utc .== t) for t in sattime]
    # Save time, lat/lon arrays, and feature classification flags (FCF) in CPro struct
    new(DataFrame(time=utc[idx], lat=[AbstractFloat[]; lat...][idx], lon=[AbstractFloat[]; lon...][idx],
      FCF=[Vector{<:Union{Missing,UInt16}}[]; avd...][idx],
      EC532=[Vector{<:Union{Missing,AbstractFloat}}[]; ec532...][idx]),
      CProMetadata(lidarrange, lidar))
  end
end #struct CPro


""" External constructor for emtpy CPro struct """
CPro() = CPro(DataFrame(time = DateTime[], lat = AbstractFloat[], lon = AbstractFloat[],
	FCF = UInt16[], EC532 = AbstractFloat[]), CProMetadata((-Inf, Inf), get_lidarheights((-Inf, Inf))))


## Define structs related to intersection data

struct Intersection
  data::DataFrame
  tracked::DataFrame
  accuracy::DataFrame
  metadata::XMetadata

  #=
  """ Unmodified constructor for `Intersection` """
  function Intersection(lat::AbstractFloat, lon::AbstractFloat, tdiff::Dates.CompoundPeriod,
    accuracy::AbstractFloat, cirrus::Bool, sat::SatDB, flight::FlightData)
    new(lat, lon, tdiff, accuracy, cirrus, sat, flight)
  end #constructor 1 Intersection
  =#

  """ Modified constructor with some automated calculations of the intersection data. """
  function Intersection(flights::FlightDB, sat::SatData, secondary::Bool=:false;
    maxtimediff::Int=30, flightspan::Int=0, satspan::Int=15,
    lidarrange::Tuple{Real,Real}=(15,-Inf), savesecondsattype::Bool=true,
    stepwidth::AbstractFloat=0.01, Xradius::Real=5000, remarks=nothing)
    # Initialise DataFrames with Intersection data and monitor start time
    tstart = Dates.now()
    Xdata = DataFrame(id=String[], lat=AbstractFloat[], lon=AbstractFloat[],
      tdiff=Dates.CompoundPeriod[], tflight = DateTime[], tsat = DateTime[],
      feature = Union{Missing,Symbol}[])
    track = DataFrame(id=String[], flight=FlightData[], CPro=CPro[], CLay=CLay[])
    accuracy = DataFrame(id=String[], intersection=AbstractFloat[], flightcoord=AbstractFloat[],
      satcoord=AbstractFloat[], flighttime=Dates.CompoundPeriod[], sattime=Dates.CompoundPeriod[])
    # Get lidar altitude levels
    lidar = get_lidarheights(lidarrange)
    # New MATLAB session
    ms = mat.MSession()
    # Loop over data and interpolate track data and time, throw error on failure
    @pm.showprogress 1 "find intersections..." for flight in
      [flights.inventory; flights.archive; flights.onlineData]
      try
        # Find sat tracks in the vicinity of flight tracks, where intersections are possible
        overlap = findoverlap(flight, sat, maxtimediff)
        isempty(overlap) && continue
        # Interpolate trajectories using MATLAB's pchip routine
        sattracks = interpolate_satdata(ms, sat, overlap, flight.metadata)
        flighttracks = interpolate_flightdata(ms, flight, stepwidth)
        # Calculate intersections and store data and metadata in DataFrames
        currdata, currtrack, curraccuracy = find_intersections(ms, flight, flighttracks,
          sat, sattracks, maxtimediff, stepwidth, Xradius, lidar, lidarrange,
          flightspan, satspan, savesecondsattype)
        append!(Xdata, currdata); append!(track, currtrack)
        append!(accuracy, curraccuracy)
      catch e
        @debug begin
          @show flight.metadata.dbID
          rethrow(e)
        end
        # Issue warning on failure of interpolating track or time data
        @warn string("Track data and/or time could not be interpolated for flight ",
          "$(flight.metadata.dbID) of $(flight.metadata.source) dataset. Data ignored.")
      end
    end #loop over flights
    # Close MATLAB session after looping over all data
    mat.close(ms)
    # Calculate load time
    tend = Dates.now()
    tc = tz.ZonedDateTime(tend, tz.localzone())
    loadtime = Dates.canonicalize(Dates.CompoundPeriod(tend - tstart))
    # Return Intersections after completion
    @info string("Intersection data loaded to properties in ",
      "$(join(loadtime.periods[1:min(2,length(loadtime.periods))], ", "))",
      "\n▪ data\n▪ tracked\n▪ accuracy\n▪ metadata")
    new(Xdata, track, accuracy,
      XMetadata(maxtimediff,stepwidth,Xradius,lidarrange,lidar,tc,loadtime,remarks))
  end #constructor Intersection
end #struct Intersection


## Export structs
export FlightDB, FlightData, SatDB, CLay, CPro, Intersection,
       FlightMetadata, SatMetadata, DBMetadata, XMetadata


## Import functions for Julia include files
include("auxiliary.jl")       # helper functions
include("lidar.jl")           # functions related to processing CALIOP data
include("loadFlightData.jl")  # functions related to loading flight databases/datasets
include("match.jl")           # functions related to finding track intersections

end # module TrackMatcher
