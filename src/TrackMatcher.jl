module TrackMatcher

# Track changes during development
# using Revise

# Import Julia packages
import CSV
import DataFrames; const df = DataFrames
import TimeZones; const tz = TimeZones
import Dates
import MATLAB; const mat = MATLAB
import Statistics; const stats = Statistics
import ProgressMeter; const pm = ProgressMeter
import Logging; const logg = Logging
# Import structs from packages
import DataFrames.DataFrame
import Dates.DateTime, Dates.Date, Dates.Time
import TimeZones.ZonedDateTime


# Define Logger with log level
logger = logg.ConsoleLogger(stdout, logg.Debug)
logg.global_logger(logger)


### Define own structs

struct PCHIP
  interpolate::Vector{<:Function}
  interval::Vector{<:NamedTuple{(:xmin,:xmax,:ymin,:ymax),<:NTuple{4,<:Float64}}}
  useLON::Bool
  session::mat.MSession

  function PCHIP(x::Vector{<:Float64}, y::Vector{<:Float64},
    flex::Vector{<:UnitRange}, id::Int64, useLON::Bool, ms::mat.MSession)
    itp = Function[]; itv = NamedTuple{(:xmin,:xmax,:ymin,:ymax),<:NTuple{4,<:Float64}}[]
    for (i, f) in enumerate(flex)
      p = "p$(id)_$i"
      mat.put_variable(ms, :x, x[f])
      mat.put_variable(ms, :y, y[f])
      mat.eval_string(ms, "$p = pchip(x, y);");
      pp = mat.get_mvariable(ms, Symbol("$p"))
      push!(itp, Minterpolate(ms, pp))
      push!(itv, (xmin=x[f[1]], xmax=x[f[end]], ymin=minimum(y[f]), ymax=maximum(y[f])))
    end

    new(itp, itv, useLON, ms)
  end
end



"""
# struct MetaData

Immutable struct to hold metadata for `FlightData` of the `FlightDB` with fields

- `dbID::Union{Int,AbstractString}`
- `flightID::Union{Missing,AbstractString}`
- `aircraft::Union{Missing,AbstractString}`
- `route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}}`
- `area::NamedTuple{(:latmin,:latmax,:plonmin,:plonmax,:nlonmin,:nlonmax),Tuple{Float64,Float64,Float64,Float64,Float64,Float64}}`
- `date::NamedTuple{(:start,:stop),Tuple{ZonedDateTime,ZonedDateTime}}`
- `file::AbstractString`

## dbID
Database ID – integer counter for `inventory` and FlightAware `onlineData`,
String with information about `FlightID`, `route`, and scheduled arrival.

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
- `plonmin`
- `plonmax`
- `nlonmin`
- `nlonmax`

## date
`NamedTuple` with fields `start` and `stop` for start and end time of the current
flight.

## file
String holding the absolute folder path and file name.


# Instantiation

    MetaData(dbID::Union{Int,AbstractString},
      flightID::Union{Missing,AbstractString}, aircraft::Union{Missing,AbstractString},
      route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
      lat::Vector{<:Union{Missing,Float64}}, lon::Vector{<:Union{Missing,Float64}},
      date::Vector{ZonedDateTime}, file::AbstractString) -> struct MetaData

Construct `MetaData` from `dbID`, `flightID`, `aircraft` type, `route`, and `file`.
Fields `area` and `date` are calculated from `lat`/`lon`, and `date` vectors.
"""
struct MetaData
  dbID::Union{Int,AbstractString}
  flightID::Union{Missing,AbstractString}
  route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}}
  aircraft::Union{Missing,AbstractString}
  date::NamedTuple{(:start,:stop),Tuple{ZonedDateTime,ZonedDateTime}}
  area::NamedTuple{(:latmin,:latmax,:plonmin,:plonmax,:nlonmin,:nlonmax),NTuple{6,Float64}}
  flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,Float64,Float64}}}}
  useLON::Bool
  file::AbstractString

  function MetaData(dbID::Union{Int,AbstractString}, flightID::Union{Missing,AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
    aircraft::Union{Missing,AbstractString}, date::Vector{ZonedDateTime},
    lat::Vector{<:Union{Missing,Float64}}, lon::Vector{<:Union{Missing,Float64}},
    useLON::Bool,
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,Float64,Float64}}}},
    file::AbstractString)

    plonmax = isempty(lon[lon.≥0]) ? NaN : maximum(lon[lon.≥0])
    plonmin = isempty(lon[lon.≥0]) ? NaN : minimum(lon[lon.≥0])
    nlonmax = isempty(lon[lon.<0]) ? NaN : maximum(lon[lon.<0])
    nlonmin = isempty(lon[lon.<0]) ? NaN : minimum(lon[lon.<0])
    area = (latmin=minimum(lat), latmax=maximum(lat),
      plonmin=plonmin, plonmax=plonmax, nlonmin=nlonmin, nlonmax=nlonmax)
    new(dbID, flightID, route, aircraft, (start=date[1], stop=date[end]), area,
      flex, useLON, file)
  end #constructor MetaData
end #struct MetaData


"""
# struct FlightData

Aircraft data with fields
- `time::Vector{ZonedDateTime}`
- `lat::Vector{<:Union{Missing,Float64}}`
- `lon::Vector{<:Union{Missing,Float64}}`
- `alt::Vector{<:Union{Missing,Float64}}`
- `heading::Vector{<:Union{Missing,Int}}`
- `climb::Vector{<:Union{Missing,Int}}`
- `speed::Vector{<:Union{Missing,Float64}}`
- `metadata::MetaData`

## time
Vector of `ZonedDateTime`

## lat/lon
Vectors of `Float64` with ranges -90°...90° and -180°...180°.

## alt
Vector of `Float64` with altitude in feet.

## heading
Vector of `Int` with course heading in degrees.

## climb
Vector of `Int` with climbing (positive) / sinking (negative) rate in feet (0 = level).

## speed
Vector of `Float64` in knots.


# Instantiation

    FlightData(time::Vector{ZonedDateTime}, lat::Vector{<:Union{Missing,Float64}},
      lon::Vector{<:Union{Missing,Float64}}, alt::Vector{<:Union{Missing,Float64}},
      heading::Vector{<:Union{Missing,Int}}, climb::Vector{<:Union{Missing,Int}},
      speed::Vector{<:Union{Missing,Float64}}, dbID::Union{Int,AbstractString},
      flightID::Union{Missing,AbstractString}, aircraft::Union{Missing,AbstractString},
      route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
      file::AbstractString) -> struct FlightData

Construct `FlightData` from fields and additonal information `dbID`, `flightID`,
`aircraft` type, `route`, and `file` name for `MetaData`.
"""
struct FlightData
  time::Vector{ZonedDateTime}
  lat::Vector{<:Union{Missing,Float64}}
  lon::Vector{<:Union{Missing,Float64}}
  alt::Vector{<:Union{Missing,Float64}}
  heading::Vector{<:Union{Missing,Int}}
  climb::Vector{<:Union{Missing,Int}}
  speed::Vector{<:Union{Missing,Float64}}
  metadata::MetaData

  function FlightData(time::Vector{ZonedDateTime}, lat::Vector{<:Union{Missing,Float64}},
    lon::Vector{<:Union{Missing,Float64}}, alt::Vector{<:Union{Missing,Float64}},
    heading::Vector{<:Union{Missing,Int}}, climb::Vector{<:Union{Missing,Int}},
    speed::Vector{<:Union{Missing,Float64}}, dbID::Union{Int,AbstractString},
    flightID::Union{Missing,AbstractString}, aircraft::Union{Missing,AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,Float64,Float64}}}},
    useLON::Bool, file::AbstractString)

    lat = checklength(lat, time)
    lon = checklength(lon, time)
    alt = checklength(alt, time)
    heading = checklength(heading, time)
    climb = checklength(climb, time)
    speed = checklength(speed, time)
    metadata = MetaData(dbID,flightID,route,aircraft,time,lat,lon,useLON,flex,file)

    new(time,lat,lon,alt,heading,climb,speed,metadata)
  end #constructor FlightData
end #struct FlightData


"""
# struct FlightDB

Database for aircraft data of different database types with fields:
- `inventory::Vector{FlightData}`
- `archive::Vector{FlightData}`
- `onlineData::Vector{FlightData}`
- `created::Union{DateTime,tz.ZonedDateTime}`
- `remarks`

## inventory
Flight data from csv files.

## archive
Commercial flight data by FlightAware.

## onlineData
Online data from FlightAware website.

## created
Time of creation as `DateTime` (or `ZonedDateTime`).

## remarks
Any data that can be attached to `FlightData` with keyword argument `remarks`.


# Instantiation

Use function `loadFlightDB` for an easy instatiation of `FlightDB`.
"""
struct FlightDB
  inventory::Vector{FlightData}
  archive::Vector{FlightData}
  onlineData::Vector{FlightData}
  created::Union{DateTime,tz.ZonedDateTime}
  remarks
end #struct FlightDB


"""
# struct CLay

CALIOP cloud layer data with fields:
- `time::Vector{ZonedDateTime}`
- `lat::Vector{Float64}`
- `lon::Vector{Float64}`

# Instantiation

    CLay(files::String...) -> struct CLay

Construct `CLay` from a list of file names (including directories).
"""
struct CLay
  time::Vector{ZonedDateTime}
  lat::Vector{Float64}
  lon::Vector{Float64}

  function CLay(folders::String...)
    # Scan folders for HDF4 files
    files = String[];
    for folder in folders
      files = findFiles(files, folder, ".hdf")
    end
    # Initialise arrays
    utc = ZonedDateTime[]; lon = []; lat = []
    # Loop over files
    @pm.showprogress 1 "load CLay data..." for file in files
      # Find files with cloud layer data
      if occursin("CLay", basename(file))
        # Extract time and convert to UTC
        t = mat.mxcall(:hdfread,1,file,"Profile_UTC_Time")[:,2]
        utc = [utc; convertUTC.(t)]
        # Extract lat/lon
        lon = [lon; mat.mxcall(:hdfread,1,file, "Longitude")[:,2]]
        lat = [lat; mat.mxcall(:hdfread,1,file, "Latitude")[:,2]]
      end
    end

    # Save time, lat/lon arrays in CLay struct
    new(utc, lat, lon)
  end #constructor CLay
end #struct CLay


"""
# struct CPro

CALIOP cloud profile data with fields:
- `time::Vector{ZonedDateTime}`
- `lat::Vector{Float64}`
- `lon::Vector{Float64}`

# Instantiation

    CPro(files::String...) -> struct CPro

Construct `CPro` from a list of file names (including directories).
"""
struct CPro
  time::Vector{ZonedDateTime}
  lat::Vector{Float64}
  lon::Vector{Float64}

  function CPro(folders::String...)
    # Scan folders for HDF4 files
    files = String[];
    for folder in folders
      files = findFiles(files, folder, ".hdf")
    end
    # Initialise arrays
    utc = ZonedDateTime[]; lon = []; lat = []
    # Loop over files
    @pm.showprogress 1 "load CPro data..." for file in files
      # Find files with cloud profile data
      if occursin("CPro", basename(file))
        # Extract time and convert to UTC
        t = mat.mxcall(:hdfread,1,file,"Profile_UTC_Time")[:,2]
        utc = [utc; convertUTC.(t)]
        # Extract lat/lon
        lon = [lon; mat.mxcall(:hdfread,1,file, "Longitude")[:,2]]
        lat = [lat; mat.mxcall(:hdfread,1,file, "Latitude")[:,2]]
      end
    end

    # Save time, lat/lon arrays in CLay struct
    new(utc, lat, lon)
  end #constructor CPro
end #struct CPro


"""
# struct SatDB

Immutable struct with fields

- `CLay::CLay`
- `CPro::CPro`
- `created::DateTime`
- `remarks`

## CLay and CPro

CALIOP satellite data currently holding time as `ZonedDateTime`
and position (`lat`/`lon`) of cloud layer and profile data.

## created

Time of creation of satellite database as `DateTime`.

## remarks
Any data can be attached to the satellite data with the keyword `remarks`.

# Instantiation

    SatDB(folder::String...; remarks=nothing) -> struct SatDB

Construct a CALIOP satellite database from HDF4 files (CALIOP version 4.x)
in `folder` or any subfolder (several folders can be given as vararg).
Attach comments or any data with keyword argument `remarks`.
"""
struct SatDB
  CLay::CLay
  CPro::CPro
  created::DateTime
  remarks

  function SatDB(folders::String...; remarks=nothing)
    cl = CLay(folders...)
    cp = CPro(folders...)
    tc = Dates.now()

    new(cl, cp, tc, remarks)
  end #constructor SatDB
end #struct SatDB

export loadFlightDB,
       FlightDB,
       FlightData,
       MetaData,
       CLay,
       CPro,
       SatDB


include("auxiliary.jl")
include("loadFlightData.jl")

end # module TrackMatcher
