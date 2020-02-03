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
- `MetaData` holds metadata to every flight
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

# Track changes during development
# using Revise

# Import Julia packages
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
import Dates.DateTime, Dates.Date, Dates.Time
import TimeZones.ZonedDateTime


# Define Logger with log level
logger = logg.ConsoleLogger(stdout, logg.Info)
logg.global_logger(logger)


### Define own structs
"""
# struct MetaData

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

`MetaData` is constructed automatically, when `FlightData` is instatiated using
a modified constructor and `dbID`, `flightID`, `aircraft` type, `route`, and `file`.
Fields `area` and `date` are calculated from `lat`/`lon`, and `date` vectors.

    MetaData(dbID::Union{Int,AbstractString},
      flightID::Union{Missing,AbstractString}, aircraft::Union{Missing,AbstractString},
      route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
      lat::Vector{<:Union{Missing,Float64}}, lon::Vector{<:Union{Missing,Float64}},
      date::Vector{DateTime}, file::AbstractString) -> struct MetaData

Or construct `MetaData` by directly handing over every field:

    MetaData(dbID::Union{Int,AbstractString}, flightID::Union{Missing,AbstractString},
      route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
      aircraft::Union{Missing,AbstractString}, date::Vector{DateTime},
      lat::Vector{<:Union{Missing,Float64}}, lon::Vector{<:Union{Missing,Float64}},
      useLON::Bool,
      flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,Float64,Float64}}}},
      file::AbstractString)
"""
struct MetaData
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

  """ Unmodified constructor for `Metadata` """
  function MetaData(dbID::Union{Int,AbstractString}, flightID::Union{Missing,AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
    aircraft::Union{Missing,AbstractString}, date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,Float64}},
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,Float64,Float64}}}},
    useLON::Bool, source::AbstractString, file::AbstractString)

    new(dbID, flightID, route, aircraft, date, area, flex, useLON, source, file)
  end #constructor 1 MetaData


  """
  Modified constructor for MetaData with some automated construction of fields
  and variable checks.
  """
  function MetaData(dbID::Union{Int,AbstractString}, flightID::Union{Missing,AbstractString},
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
  end #constructor 2 MetaData
end #struct MetaData


"""
# struct FlightData

Aircraft data with fields
- `data::DataFrame`
- `metadata::MetaData`

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
`dbID`, `flightID`, `aircraft` type, `route`, and `file` name for `MetaData`.
For `data`, `time` is the only exception, which is given as `ZonedDateTime` and
will be converted to `UTC` standard time.

Or construct by directly handing over every field:

    FlightData(time::Vector{DateTime}, lat::Vector{<:Union{Missing,Float64}},
      lon::Vector{<:Union{Missing,Float64}}, alt::Vector{<:Union{Missing,Float64}},
      heading::Vector{<:Union{Missing,Int}}, climb::Vector{<:Union{Missing,Int}},
      speed::Vector{<:Union{Missing,Float64}}, metadata::MetaData)

Checks exist that the order, names, and types of the `data` `DataFrame` are correct.
"""
struct FlightData
  data::DataFrame
  metadata::MetaData

  """ Unmodified constructor for `FlightData` with basic checks for correct `data`"""
  function FlightData(data::DataFrame, metadata::MetaData)

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
  function FlightData(time::Vector{ZonedDateTime}, lat::Vector{<:Union{Missing,Float64}},
    lon::Vector{<:Union{Missing,Float64}}, alt::Vector{<:Union{Missing,Float64}},
    heading::Vector{<:Union{Missing,Int}}, climb::Vector{<:Union{Missing,Int}},
    speed::Vector{<:Union{Missing,Float64}}, dbID::Union{Int,AbstractString},
    flightID::Union{Missing,AbstractString}, aircraft::Union{Missing,AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,Float64,Float64}}}},
    useLON::Bool, source::String, file::AbstractString)

    t = [t.utc_datetime for t in time]
    lat = checklength(lat, t)
    lon = checklength(lon, t)
    alt = checklength(alt, t)
    heading = checklength(heading, t)
    climb = checklength(climb, t)
    speed = checklength(speed, t)
    metadata = MetaData(dbID,flightID,route,aircraft,t,lat,lon,useLON,flex,source,file)

    new(DataFrame(time=t,lat=lat,lon=lon,alt=alt,heading=heading,climb=climb,speed=speed),metadata)
  end #constructor 2 FlightData
end #struct FlightData


"""
# struct FlightDB

Database for aircraft data of different database types with fields:
- `inventory::Vector{FlightData}`
- `archive::Vector{FlightData}`
- `onlineData::Vector{FlightData}`
- `created::Union{DateTime,ZonedDateTime}`
- `remarks`

## inventory
Flight data from csv files.

## archive
Commercial flight data by FlightAware.

## onlineData
Online data from FlightAware website.

## created
Time of creation as `DateTime` or `ZonedDateTime` (default).

## remarks
Any data that can be attached to `FlightData` with keyword argument `remarks`.


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
  created::Union{DateTime,ZonedDateTime}
  remarks

  """
  Unmodified constructor for `FlightDB` with basic checks for correct dataset type
  in each dataset field and presets for field created (now) and remarks (nothing).
  """
  function FlightDB(inventory::Vector{FlightData},
    archive::Vector{FlightData}, onlineData::Vector{FlightData},
    created::Union{DateTime,ZonedDateTime}=tz.now(tz.localzone()),
    remarks=nothing)

    inventory = checkDBtype(inventory, "VOLPE AEDT")
    archive = checkDBtype(archive, "FlightAware")
    onlineData = checkDBtype(onlineData, "flightaware.com")

    new(inventory, archive, onlineData, created, remarks)
  end #constructor 1 FlightDB

  """
  Modified constructor creating the database from an identifer of the
  database type and the respective folder path for that database.
  """
  function FlightDB(DBtype::String, folder::Union{String, Vector{String}}...;
    altmin::Int=15_000, remarks=nothing, odelim::Union{Nothing,Char,String}=nothing)

    # Check DBtype addresses all folder paths
    if length(DBtype) ≠ length(folder)
      throw(ArgumentError("Number of characters in `DBtype` must match length of vararg `folder`"))
    end
    # Save time of database creation
    tc = tz.now(tz.localzone())
    # Find database types
    i1 = [findall(isequal('i'), DBtype); findall(isequal('1'), DBtype)]
    i2 = [findall(isequal('a'), DBtype); findall(isequal('2'), DBtype)]
    i3 = [findall(isequal('o'), DBtype); findall(isequal('3'), DBtype)]

    # Load databases for each type
    # VOLPE AEDT inventory
    ifiles = String[]
    for i in i1
      ifiles = findFiles(ifiles, folder[i], ".csv")
    end
    inventory = loadInventory(ifiles, altmin=altmin)
    # FlightAware commercial archive
    ifiles = String[]
    for i in i2
      ifiles = findFiles(ifiles, folder[i], ".csv")
    end
    archive = loadArchive(ifiles, altmin=altmin)
    ifiles = String[]
    for i in i3
      ifiles = findFiles(ifiles, folder[i], ".txt", ".dat")
    end
    onlineData = loadOnlineData(ifiles, altmin=altmin, delim=odelim)

    @info "done loading data to properties\n- inventory\n- archive\n- onlineData"

    new(inventory, archive, onlineData, tc, remarks)
  end # constructor 2 FlightDB
end #struct FlightDB


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
  function CLay(ms::mat.MSession, folders::String...)
    # Scan folders for HDF4 files
    files = String[];
    for folder in folders
      files = findFiles(files, folder, ".hdf")
    end
    # Initialise arrays
    utc = DateTime[]; lon = Float64[]; lat = Float64[]
    # Loop over files
    @pm.showprogress 1 "load CLay data..." for file in files
      # Find files with cloud layer data
      if occursin("CLay", basename(file))
        # Extract time and convert to UTC
        t = mat.mxcall(ms, :hdfread,1,file,"Profile_UTC_Time")[:,2]
        utc = [utc; convertUTC.(t)]
        # Extract lat/lon
        lon = [lon; mat.mxcall(ms, :hdfread,1,file, "Longitude")[:,2]]
        lat = [lat; mat.mxcall(ms, :hdfread,1,file, "Latitude")[:,2]]
      end
    end

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
    standardnames = [:time, :lat, :lon]
    standardtypes = [Vector{DateTime}, Vector{Float64}, Vector{Float64}]
    bounds = Tuple{Real,Real}[]
    data = checkcols(data, standardnames, standardtypes, bounds, "CPro", nothing)
    new(data)
  end #constructor 1 CPro

  """
  Modified constructor of `CPro` reading data from hdf files given in `folders...`
  using MATLAB session `ms`.
  """
  function CPro(ms::mat.MSession, folders::String...)
    # Scan folders for HDF4 files
    files = String[];
    for folder in folders
      files = findFiles(files, folder, ".hdf")
    end
    # Initialise arrays
    utc = DateTime[]; lon = Float64[]; lat = Float64[]
    # Loop over files
    @pm.showprogress 1 "load CPro data..." for file in files
      # Find files with cloud profile data
      if occursin("CPro", basename(file))
        # Extract time and convert to UTC
        t = mat.mxcall(ms, :hdfread,1,file,"Profile_UTC_Time")[:,2]
        utc = [utc; convertUTC.(t)]
        # Extract lat/lon
        lon = [lon; mat.mxcall(ms, :hdfread,1,file, "Longitude")[:,2]]
        lat = [lat; mat.mxcall(ms, :hdfread,1,file, "Latitude")[:,2]]
      end
    end

    # Save time, lat/lon arrays in CLay struct
    new(DataFrame(time=utc, lat=lat, lon=lon))
  end #constructor 2 CPro
end #struct CPro


"""
# struct SatDB

Immutable struct with fields

- `CLay::CLay`
- `CPro::CPro`
- `created::Union{DateTime,ZonedDateTime}`
- `remarks`

## CLay and CPro

CALIOP satellite data currently holding time as `DateTime` in `UTC`
and position (`lat`/`lon`) of cloud layer and profile data in a `DataFrame`.

## created

Time of creation of satellite database as `DateTime` or with timezone information
as `ZonedDateTime` (default).

## remarks
Any data can be attached to the satellite data with the keyword `remarks`.

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
  created::Union{DateTime,ZonedDateTime}
  remarks

  """ Unmodified constructor for `SatDB` """
  function SatDB(CLay::CLay, CPro::CPro,
    created::Union{DateTime,ZonedDateTime}=tz.now(tz.localzone()), remarks=nothing)
    new(CLay, CPro, created, remarks)
  end #constructor 1 SatDb

  """
  Automated constructor scanning for `HDF4` in `folders`; any data or comments
  can be attached in the field remarks.
  """
  function SatDB(folders::String...; remarks=nothing)
    ms = mat.MSession()
    cl = CLay(ms, folders...)
    cp = CPro(ms, folders...)
    mat.close(ms)
    tc = tz.now(tz.localzone())

    @info "done loading data to properties\n- CLay\n- CPro"
    new(cl, cp, tc, remarks)
  end #constructor 2 SatDB
end #struct SatDB


"""
# struct Intersection

Immutable struct with fields

- `lat::Float64`
- `lon::Float64`
- `tdiff::Dates.CompoundPeriod`
- `accuracy::Float64`
- `cirrus::Bool`
- `sat::SatDB`
- `flight::FlightData`

## lat/lon

Position of intersection in degrees.


## tdiff

Time difference between aircraft and satellite overpass at intersection.
Positive time differences mean satellite overpass before flight overpass,
negative times mean flight reaches intersection before satellite.


## accuracy

Accuracy during of the interpolation in meters.


## cirrus

Flag for cirrus clouds at flight level at the intersection.


## sat

Satellite cloud layer and profile data (as available) stored in `SatDB` for ±15
time steps.


## flight

`FlightData` for the uninterpolated timepoint closest to the overpass at the intersection.


# Instatiation

Instantiate with position (`lat`/`lon`) of Intersection, the time difference `tdiff`
of the flight and satellite overpass, the `accuracy` of the interpolation, a flag
for `cirrus` clouds at flight level at the intersection, `SatDB` with any available
`CLay` and `CPro` data in the vicinity of the intersection (default ± 15 time steps),
and `FlightData` with the closest measured point to the intersection.

    Intersection(flight::FlightData, sat::SatDB, sattype::Symbol,
      tflight::DateTime, tsat::DateTime, lat::Float64, lon::Float64,
      flightspan::Int, satspan::Int, accuracy::Float64)

Or construct directly from the fields in `Intersection`:

    Intersection(lat::Float64, lon::Float64, tdiff::Dates.CompoundPeriod,
      accuracy::Float64, cirrus::Bool, sat::SatDB, flight::FlightData)
"""
struct Intersection
  # tflight::DateTime
  # tsat::DateTime
  lat::Float64
  lon::Float64
  tdiff::Dates.CompoundPeriod
  accuracy::Float64
  cirrus::Bool
  sat::SatDB
  flight::FlightData

  """ Unmodified constructor for `Intersection` """
  function Intersection(lat::Float64, lon::Float64, tdiff::Dates.CompoundPeriod,
    accuracy::Float64, cirrus::Bool, sat::SatDB, flight::FlightData)
    new(lat, lon, tdiff, accuracy, cirrus, sat, flight)
  end #constructor 1 Intersection

  """ Modified constructor with some automated calculations of the intersection data. """
  function Intersection(flight::FlightData, sat::SatDB, sattype::Symbol,
    tflight::DateTime, tsat::DateTime, lat::Float64, lon::Float64,
    flightspan::Int, satspan::Int, accuracy::Float64)

    # Calculate time difference between flight and satellite overpass at intersection
    tdiff = Dates.canonicalize(Dates.CompoundPeriod(tflight - tsat))
    # Find the index (DataFrame row) of the intersection in the flight data
    tf = argmin(abs.(flight.data.time .- tflight))

    # Construct FlightData at Intersection
    flightdata =
      FlightData(extract_timespan(flight.data, tf, flightspan), flight.metadata)

    # Get satellite data used to find the intersection and find DataFrame row of intersection
    satprim = getfield(sat, sattype).data
    tsp = argmin(abs.(satprim.time .- tsat))
    # Retrieve DataFrame at Intersection ± 15 time steps
    primdata = extract_timespan(satprim, tsp, satspan)

    # Switch to other satellite data (CLay/CPro) and check whether data is available
    # at intersection and retrieve ±15 time steps as well
    sattype = swap_sattype(sattype)
    secdata = try satsec = getfield(sat, sattype).data
      tss = argmin(abs.(satsec.time .- tsat))
      extract_timespan(satsec, tss, satspan)
    catch
      # Return an empty DataFrame, if no data is available
      DataFrame(time=DateTime[], lat=Float64[], lon=Float64[])
    end
    # Save satellite data in SatDB
    satdb = sattype == :CPro ? SatDB(CLay(primdata), CPro(secdata), sat.created, sat.remarks) :
      SatDB(CLay(secdata), CPro(primdata), sat.created, sat.remarks)

    # Instatiate new Intersection
    new(lat, lon, tdiff, accuracy, false, satdb, flightdata)
  end
end #constructor 2 Intersection

# Needed for julia 1.0.x?:
# SatDB(CLay::CLay, CPro::CPro, created::Union{DateTime,ZonedDateTime}) = SatDB(CLay, CPro, created, nothing)
# SatDB(CLay::CLay, CPro::CPro) = SatDB(CLay, CPro, tz.now(tz.localtime()), nothing)


export FlightDB,
       FlightData,
       MetaData,
       CLay,
       CPro,
       SatDB,
       Intersection,
       intersection


include("auxiliary.jl")       # helper functions
include("loadFlightData.jl")  # functions related to loading flight databases/datasets
include("match.jl")           # functions related to finding track intersections

end # module TrackMatcher
