"""
# Module TrackMatcher

To find intersection between different trajectories. The module is aimed to find
intersections between aircraft and satellite tracks, but can be used for flight
or cloud tracks as well.

## Public structs

- `FlightDB` stores flight track data and other relevant aircraft related data
  from 3 different inventories:
  - `inventory`: VOLPE AEDT inventory
  - `archive`: commercially available database by FlightAware
  - `onlineData`: free online data by FlightAware
- `FlightData` stores `FlightDB` data of a single flight
- `FlightMetadata` holds metadata to every flight
- `SatDB` stores CALIPSO cloud layer and profile data from the CALIOP satellite
- `CLay` CALIPSO cloud layer data
- `CPro` CALIPSO cloud profile data
- `Intersection`: position of intersection between aircraft and satellite trajectory,
  together with time difference between crossing and `FlightData` and `SatDB` in the
  vicinity of the intersection


## Public functions

- `intersection` finds intersections in the trajectories of aircrafts and satellites
  stored in `FlightDB` and `SatDB` and returns a vector of `Intersection` instances
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
- `area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),Tuple{Float64,Float64,Float64,Float64,Float64,Float64}}`
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

## file
String holding the absolute folder path and file name.


# Instantiation

`FlightMetadata` is constructed automatically, when `FlightData` is instatiated using
a modified constructor and `dbID`, `flightID`, `aircraft` type, `route`, and `file`.
Fields `area` and `date` are calculated from `lat`/`lon`, and `date` vectors.

    FlightMetadata(dbID::Union{Int,AbstractString},
      flightID::Union{Missing,AbstractString}, aircraft::Union{Missing,AbstractString},
      route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
      lat::Vector{<:Union{Missing,Float64}}, lon::Vector{<:Union{Missing,Float64}},
      date::Vector{DateTime}, file::AbstractString) -> struct FlightMetadata

Or construct `FlightMetadata` by directly handing over every field:

    FlightMetadata(dbID::Union{Int,AbstractString}, flightID::Union{Missing,AbstractString},
      route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
      aircraft::Union{Missing,AbstractString}, date::Vector{DateTime},
      lat::Vector{<:Union{Missing,Float64}}, lon::Vector{<:Union{Missing,Float64}},
      useLON::Bool,
      flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,Float64,Float64}}}},
      file::AbstractString)
"""
struct FlightMetadata
  dbID::Union{Int,AbstractString}
  flightID::Union{Missing,AbstractString}
  route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}}
  aircraft::Union{Missing,AbstractString}
  date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,Float64}}
  flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,Float64,Float64}}}}
  useLON::Bool
  source::AbstractString
  file::AbstractString

  """ Unmodified constructor for `FlightMetadata` """
  function FlightMetadata(dbID::Union{Int,AbstractString}, flightID::Union{Missing,AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
    aircraft::Union{Missing,AbstractString}, date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,Float64}},
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,Float64,Float64}}}},
    useLON::Bool, source::AbstractString, file::AbstractString)

    new(dbID, flightID, route, aircraft, date, area, flex, useLON, source, file)
  end #constructor 1 FlightMetadata

  """
  Modified constructor for FlightMetadata with some automated construction of fields
  and variable checks.
  """
  function FlightMetadata(dbID::Union{Int,AbstractString}, flightID::Union{Missing,AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
    aircraft::Union{Missing,AbstractString}, date::Vector{DateTime},
    lat::Vector{<:Union{Missing,Float64}}, lon::Vector{<:Union{Missing,Float64}},
    useLON::Bool,
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,Float64,Float64}}}},
    source::AbstractString, file::AbstractString)

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

- `cirrus`: flag for cirrus clouds at flight altitude level
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
  stepwidth::Float64
  Xradius::Real
  preferred::Symbol
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

- `time::Vector{DateTime}`                  (time stamp of measurement)
- `lat::Vector{<:Union{Missing,Float64}}`   (latitude)
- `lon::Vector{<:Union{Missing,Float64}}`   (longitude)
- `alt::Vector{<:Union{Missing,Float64}}`   (altitude)
- `heading::Vector{<:Union{Missing,Int}}`   (heading/direction of the aircraft)
- `climb::Vector{<:Union{Missing,Int}}`     (climbing (positive)/sinking (negative values) rate)
- `speed::Vector{<:Union{Missing,Float64}}` (velocity)


# Instantiation

    FlightData(time::Vector{ZonedDateTime}, lat::Vector{<:Union{Missing,Float64}},
      lon::Vector{<:Union{Missing,Float64}}, alt::Vector{<:Union{Missing,Float64}},
      heading::Vector{<:Union{Missing,Int}}, climb::Vector{<:Union{Missing,Int}},
      speed::Vector{<:Union{Missing,Float64}}, dbID::Union{Int,AbstractString},
      flightID::Union{Missing,AbstractString}, aircraft::Union{Missing,AbstractString},
      route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
      file::AbstractString) -> struct FlightData

Construct `FlightData` from columns for the `data` `DataFrame` and additonal information
`dbID`, `flightID`, `aircraft` type, `route`, and `file` name for `FlightMetadata`.
For `data`, `time` is the only exception, which is given as `ZonedDateTime` and
will be converted to `UTC` standard time.

Or construct by directly handing over every field:

    FlightData(time::Vector{DateTime}, lat::Vector{<:Union{Missing,Float64}},
      lon::Vector{<:Union{Missing,Float64}}, alt::Vector{<:Union{Missing,Float64}},
      heading::Vector{<:Union{Missing,Int}}, climb::Vector{<:Union{Missing,Int}},
      speed::Vector{<:Union{Missing,Float64}}, metadata::FlightMetadata)

Checks exist that the order, names, and types of the `data` `DataFrame` are correct.
"""
struct FlightData
  data::DataFrame
  metadata::FlightMetadata

  """ Unmodified constructor for `FlightData` with basic checks for correct `data`"""
  function FlightData(data::DataFrame, metadata::FlightMetadata)

    # Column checks and warnings
    standardnames = [:time, :lat, :lon, :alt, :heading, :climb, :speed]
    standardtypes = [Union{DateTime,Vector{DateTime}}, Union{Float64,Vector{Float64}},
      Union{Float64,Vector{Float64}}, Union{Missing,Float64,Vector{<:Union{Missing,Float64}}},
      Union{Missing,Int,Vector{<:Union{Missing,Int}}},
      Union{Missing,Int,Vector{<:Union{Missing,Int}}},
      Union{Missing,Float64,Vector{<:Union{Missing,Float64}}}]
    bounds = [(0,Inf), (0, 360), (-Inf, Inf), (0, Inf)]
    data = checkcols(data, standardnames, standardtypes, bounds,
      metadata.source, metadata.dbID)
    new(data,metadata)
  end #constructor 1 FlightData

  """ Modified constructor with variable checks and some automated calculation of fields """
  function FlightData(flightdata::DataFrame, dbID::Union{Int,AbstractString},
    flightID::Union{Missing,AbstractString}, aircraft::Union{Missing,AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,Float64,Float64}}}},
    useLON::Bool, source::String, file::AbstractString)

    # Check dataframe columns of flight data; fill missing columns with missing values
    t = getproperty(flightdata, :time)
    t = t isa Vector{ZonedDateTime} ? [zt.utc_datetime for zt in t] : t
    lat = getproperty(flightdata, :lat)
    lon = getproperty(flightdata, :lon)
    alt = getproperty(flightdata, :alt)
    heading = hasproperty(flightdata, :heading) ? getproperty(flightdata, :heading) :
      [missing for i in t]
    climb = hasproperty(flightdata, :climb) ? getproperty(flightdata, :climb) :
      [missing for i in t]
    speed = hasproperty(flightdata, :speed) ? getproperty(flightdata, :speed) :
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
"""
# struct CLay

CALIOP cloud layer `data` stored in a `DataFrame` with columns:
- `time::Vector{DateTime}`
- `lat::Vector{Float64}`
- `lon::Vector{Float64}`

# Instantiation

    CLay(ms::mat.MSession, files::String...) -> struct CLay

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
    standardtypes = [Vector{DateTime}, Vector{Float64}, Vector{Float64}]
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
    utc = DateTime[]; lon = Float64[]; lat = Float64[]
    # Loop over files
    prog = pm.Progress(length(files), "load CLay data...")
    for file in files
      # Find files with cloud layer data
      utc, lon, lat = try
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
        [utc; convertUTC.(t)],
        [lon; longitude],
        [lat; latitude]
      catch
        # Skip data on failure and warn
        @warn string("read error in CALIPSO granule ",
          "$(splitext(basename(file))[1])\ndata skipped")
        utc, lon, lat
      end
      # Monitor progress for progress bar
      pm.next!(prog, showvalues = [(:date,Dates.Date(splitdir(dirname(file))[2], "y_m_d"))])
    end #loop over files
    pm.finish!(prog)

    # Save time, lat/lon arrays in CLay struct
    new(DataFrame(time=utc, lat=lat, lon=lon))
  end #constructor 2 CLay
end #struct CLay


"""
# struct CPro

CALIOP cloud profile `data` stored in a `DataFrame` with columns:
- `time::Vector{DateTime}`
- `lat::Vector{Float64}`
- `lon::Vector{Float64}`

# Instantiation

    CPro(ms::mat.MSession, files::String...) -> struct CPro

Construct `CPro` from a list of file names (including directories) and a running
MATLAB session.

Or construct `CPro` by directly handing over the `DataFrame` where the names, order,
and types of each columns are checked and attempted to correct:

    CPro(data::DataFrame) -> struct CPro
"""
struct CPro
  data::DataFrame

  """ Unmodified constructor for `CPro` """
  function CPro(data::DataFrame)
    standardnames = [:time, :lat, :lon, :FCF, :EC532]
    standardtypes = [Vector{DateTime}, Vector{Float64}, Vector{Float64},
      Vector{Vector{<:Union{Missing,UInt16}}}, Vector{Vector{<:Union{Missing,Float32}}}]
    bounds = Tuple{Real,Real}[]
    data = checkcols(data, standardnames, standardtypes, bounds, "CPro", nothing)
    new(data)
  end #constructor 1 CPro

  """
  Modified constructor of `CPro` reading data from hdf files given in `folders...`
  using MATLAB session `ms`.
  """
  function CPro(ms::mat.MSession, files::Vector{String}, lidar::NamedTuple)
    # Initialise arrays
    # essential data
    utc = DateTime[]; lon = Float64[]; lat = Float64[]
    # non-essential data
    avd = Vector{<:Union{Missing,UInt16}}[]; ec532 = Vector{<:Union{Missing,Float32}}[]
    # Loop over files
    prog = pm.Progress(length(files), "load CPro data...")
    # Loop over files with cloud profile data
    for file in files
      # Retrieve essential data
      utc, lon, lat = try
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
        [utc; convertUTC.(t)],
        [lon; longitude],
        [lat; latitude]
      catch
        # Skip data on failure and warn
        @warn string("read error in CALIPSO granule ",
          "$(splitext(basename(file))[1])\ndata skipped")
        utc, lon, lat
      end
      # Retrieve non-essential data
      # Extract feature classification flags
      avd = try append_lidardata!(avd, ms, "Atmospheric_Volume_Description", lidar)
      catch
        [avd; [missing for i = 1:length(lidar.refined)]]
      end
      ec532 = try append_lidardata!(ec532, ms, "Extinction_Coefficient_532", lidar,
        true, missingvalues = -9999)
        catch
          [ec532; [missing for i = 1:length(lidar.coarse)]]
        end
      # Monitor progress for progress bar
      pm.next!(prog, showvalues = [(:date,Dates.Date(splitdir(dirname(file))[2], "y_m_d"))])
    end #loop over files
    pm.finish!(prog)

    # Save time, lat/lon arrays, and feature classification flags (FCF) in CPro struct
    new(DataFrame(time=utc, lat=lat, lon=lon, FCF=avd, EC532=ec532))
  end #constructor 2 CPro
end #struct CPro


"""
# struct SatDB

Immutable struct with fields

- `CLay::CLay`
- `CPro::CPro`
- `metadata::DBMetadata`


## CLay and CPro

CALIOP satellite data currently holding time as `DateTime` in `UTC`
and position (`lat`/`lon`) of cloud layer and profile data in a `DataFrame`.

## metadata

Immutable struct of type `DBMetadata` with fields:
- `created::Union{ZonedDateTime,DateTime}`: time of creation
- `loadtime::CompoundPeriod`: time it to to read data files and load data to struct
- `remarks`: Any additional data or comments

# Instantiation

    SatDB(folder::String...; remarks=nothing) -> struct SatDB

Construct a CALIOP satellite database from HDF4 files (CALIOP version 4.x)
in `folder` or any subfolder (several folders can be given as vararg).
Attach comments or any data with keyword argument `remarks`.

Or construct by directly handing over struct fields (remarks are an optional
argument defaulting to `nothing`):

    SatDB(CLay::CLay, CPro::CPro, created::Union{DateTime,ZonedDateTime}, remarks=nothing)
"""
struct SatDB
  CLay::CLay
  CPro::CPro
  metadata::DBMetadata

  """ Unmodified constructor for `SatDB` """
  function SatDB(clay::CLay, cpro::CPro, metadata::DBMetadata)
    new(clay, cpro, metadata)
  end #constructor 1 SatDb

  """
  Automated constructor scanning for `HDF4` in `folders`; any data or comments
  can be attached in the field remarks.
  """
  function SatDB(folders::String...; remarks=nothing)
    tstart = Dates.now()
    # Scan folders for HDF4 files
    files = String[];
    for folder in folders
      findfiles!(files, folder, ".hdf")
    end
    # Start MATLAB session
    ms = mat.MSession()
    # Get lidar altitude levels
    lidar = get_lidarheights()
    # Load CALIPSO cloud layer and profile data based on file names
    clay = CLay(ms, files[occursin.("CLay", basename.(files))])
    cpro = CPro(ms, files[occursin.("CPro", basename.(files))], lidar)
    # Close MATLAB session
    mat.close(ms)
    # Calculate time span of satellite data
    tmin = minimum([clay.data.time; cpro.data.time])
    tmax = maximum([clay.data.time; cpro.data.time])
    tend = Dates.now()
    # Save computing times
    tc = tz.ZonedDateTime(tend, tz.localzone())
    loadtime = Dates.canonicalize(Dates.CompoundPeriod(tend - tstart))

    # Instantiate new struct
    @info string("SatDB data loaded to properties in ",
      "$(join(loadtime.periods[1:min(2,length(loadtime.periods))], ", "))",
      "\n▪ CLay\n▪ CPro\n▪ metadata")
    new(clay, cpro, DBMetadata((start=tmin, stop=tmax), tc, loadtime, remarks))
  end #constructor 2 SatDB
end #struct SatDB


## Define structs related to intersection data
"""
# struct Intersection

Immutable struct with fields

- `data::DataFrame`
- `tracked::DataFrame`
- `accuracy::DataFrame`
- `metadata::XMetadata`

## data

DataFrame `data` holds the interpolated spatial and temporal coordinates of all
calcalated intersections in the current dataset in columns:

- `id::Vector{String}`
- `lat::Vector{Float64}`
- `lon::Vector{Float64}`
- `tdiff::Vector{Dates.CompoundPeriod}`
- `tflight::Vector{DateTime}`
- `tsat::Vector{DateTime}`


## tracked

DataFrame `tracked` holds the actual measured flight and satellite data closest
to the intersection with additional ±`flightspan` and ±`satspan` datapoints in
columns:

- `id::Vector{String}` (same as in `data`)
- `flight::Vector{FlightData}`
- `sat::Vector{SatDB}`


## accuracy

DataFrame `accuracy` list the accuracy achived during the calculation of the
intersections in columns:

- `id::Vector{String}`
- `intersection::Vector{Float64}`: distance between calculated intersection using
  either sat or flight data
- `flightcoord`:: distance between calculated intersection and closest flight measurement
- `satcoord`:: distance between calculated intersection and closest satellite measurement
- `flighttime`:: difference between calculated time of aircraft at intersection and
  time of the closest flight measurement
- `sattime`:: difference between calculated time of satellite overpass at intersection and
  time of the closest satellite measurement


## metadata

`XMetadata` with flags used during instatiation and time of creation and computation time
for reproducable results.


# Instatiation

    Intersection(flights::FlightDB, sat::SatDB, sattype::Symbol=:CLay;
      maxtimediff::Int=30, flightspan::Int=0, satspan::Int=15, stepwidth::Float64=0.01,
      Xradius::Real=5000, remarks=nothing)

Instantiate with the `flights` and `sat` database and the `sattype` of the sat data
preferred for the calculations (`CLay` by default). The following kwargs
(with default values) can be used to define parameters of the calculation:

- `maxtimediff::Int=30`: maximum time difference in minutes allowed between satellite overpass and
  aircraft passing at intersection
- `flightspan::Int=0`: ± additional measurement points to the measurement closest to the intersection
  saved in `tracked.flight`
- `satspan::Int=15`: ± additional measurement points to the measurement closest to the intersection
  saved in `tracked.sat`
- `stepwidth`: stepwidth in degrees (lat/lon at equator) used for the interpolation of track data
- `Xradius`: radius in meters around an intersection in which further intersections
  will be removed as duplicates due to the interpolation algorithm
- remarks: any additional data or comments attached to the metadata of the struct
"""
struct Intersection
  data::DataFrame
  tracked::DataFrame
  accuracy::DataFrame
  metadata::XMetadata

  #=
  """ Unmodified constructor for `Intersection` """
  function Intersection(lat::Float64, lon::Float64, tdiff::Dates.CompoundPeriod,
    accuracy::Float64, cirrus::Bool, sat::SatDB, flight::FlightData)
    new(lat, lon, tdiff, accuracy, cirrus, sat, flight)
  end #constructor 1 Intersection
  =#

  """ Modified constructor with some automated calculations of the intersection data. """
  function Intersection(flights::FlightDB, sat::SatDB, sattype::Symbol=:CLay;
    maxtimediff::Int=30, flightspan::Int=0, satspan::Int=15, stepwidth::Float64=0.01,
    Xradius::Real=5000, remarks=nothing)
    # Initialise DataFrames with Intersection data and monitor start time
    tstart = Dates.now()
    idata = DataFrame(id=String[], lat=Float64[], lon=Float64[],
      tdiff=Dates.CompoundPeriod[], tflight = DateTime[], tsat = DateTime[],
      feature = Symbol[])
    track = DataFrame(id=String[], flight=FlightData[], sat=SatDB[])
    accuracy = DataFrame(id=String[], intersection=Float64[], flightcoord=Float64[],
      satcoord=Float64[], flighttime=Dates.CompoundPeriod[], sattime=Dates.CompoundPeriod[])
    # New MATLAB session
    ms = mat.MSession()
    # Loop over data and interpolate track data and time, throw error on failure
    @pm.showprogress 1 "find intersections..." for flight in
      [flights.inventory; flights.archive; flights.onlineData]
      try
        # Find sat tracks in the vicinity of flight tracks, where intersections are possible
        overlap = findoverlap(flight, sat, sattype, maxtimediff)
        isempty(overlap.ranges) && continue
        # Interpolate trajectories using MATLAB's pchip routine
        sattracks = interpolate_satdata(ms, sat, overlap, flight.metadata)
        flighttracks = interpolate_flightdata(ms, flight, stepwidth)
        # Calculate intersections and store data and metadata in DataFrames
        currdata, currtrack, curraccuracy = find_intersections(flight, flighttracks,
          sat, overlap.type, sattracks, maxtimediff, stepwidth, Xradius, flightspan, satspan)
        append!(idata, currdata); append!(track, currtrack)
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
    new(idata, track, accuracy,
      XMetadata(maxtimediff,stepwidth,Xradius,sattype,tc,loadtime,remarks))
  end #constructor Intersection
end #struct Intersection


## Export structs
export FlightDB, FlightData, SatDB, CLay, CPro, Intersection,
       FlightMetadata, DBMetadata, XMetadata


## Import functions for Julia include files
include("auxiliary.jl")       # helper functions
include("lidar.jl")           # functions related to processing CALIOP data
include("loadFlightData.jl")  # functions related to loading flight databases/datasets
include("match.jl")           # functions related to finding track intersections

end # module TrackMatcher
