module TrackMatcher

# Track changes during development
# using Revise

# Import Julia packages
import CSV
import DataFrames; const df = DataFrames
import TimeZones; const tz = TimeZones
import Dates
import MATLAB; const mat = MATLAB
import PyCall; const py = PyCall
# import PyPlot; const plt = PyPlot
import ProgressMeter; const pm = ProgressMeter
# Import structs from packages
import DataFrames.DataFrame
import Dates.DateTime, Dates.Date, Dates.Time
import TimeZones.ZonedDateTime
# Import interpolations from scipy for PCHIP interpolations
# Conda.add("pyplot")
const itp = py.PyNULL()
# const plt = py.PyNULL()
function __init__()
  copy!(itp, py.pyimport_conda("scipy.interpolate", "scipy"))
  # copy!(plt, py.pyimport("matplotlib.pyplot"))
end


### Define own structs

struct MetaData
  dbID::Union{Int,AbstractString}
  flightID::Union{Missing,AbstractString}
  aircraft::Union{Missing,AbstractString}
  route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}}
  area::NamedTuple{(:latmin,:latmax,:plonmin,:plonmax,:nlonmin,:nlonmax),
        Tuple{AbstractFloat,AbstractFloat,AbstractFloat,AbstractFloat,AbstractFloat,AbstractFloat}}
  date::NamedTuple{(:start,:stop),Tuple{ZonedDateTime,ZonedDateTime}}
  file::AbstractString

  function MetaData(dbID::Union{Int,AbstractString},
    flightID::Union{Missing,AbstractString}, aircraft::Union{Missing,AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
    lat::Vector{<:Union{Missing,AbstractFloat}}, lon::Vector{<:Union{Missing,AbstractFloat}},
    date::Vector{ZonedDateTime}, file::AbstractString)
    plonmax = isempty(lon[lon.≥0]) ? NaN : maximum(lon[lon.≥0])
    plonmin = isempty(lon[lon.≥0]) ? NaN : minimum(lon[lon.≥0])
    nlonmax = isempty(lon[lon.<0]) ? NaN : maximum(lon[lon.<0])
    nlonmin = isempty(lon[lon.<0]) ? NaN : minimum(lon[lon.<0])
    area = (latmin=minimum(lat), latmax=maximum(lat),
      plonmin=plonmin, plonmax=plonmax, nlonmin=nlonmin, nlonmax=nlonmax)
    new(dbID, flightID, aircraft, route, area, (start=date[1], stop=date[end]), file)
  end #constructor MetaData
end #struct MetaData

struct FlightData
  time::Vector{ZonedDateTime}
  lat::Vector{<:Union{Missing,AbstractFloat}}
  lon::Vector{<:Union{Missing,AbstractFloat}}
  alt::Vector{<:Union{Missing,AbstractFloat}}
  heading::Vector{<:Union{Missing,Int}}
  climb::Vector{<:Union{Missing,Int}}
  speed::Vector{<:Union{Missing,AbstractFloat}}
  metadata::MetaData

  function FlightData(time::Vector{ZonedDateTime}, lat::Vector{<:Union{Missing,AbstractFloat}},
    lon::Vector{<:Union{Missing,AbstractFloat}}, alt::Vector{<:Union{Missing,AbstractFloat}},
    heading::Vector{<:Union{Missing,Int}}, climb::Vector{<:Union{Missing,Int}},
    speed::Vector{<:Union{Missing,AbstractFloat}}, dbID::Union{Int,AbstractString},
    flightID::Union{Missing,AbstractString}, aircraft::Union{Missing,AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
    file::AbstractString)

    lat = checklength(lat, time)
    lon = checklength(lon, time)
    alt = checklength(alt, time)
    heading = checklength(heading, time)
    climb = checklength(climb, time)
    speed = checklength(speed, time)
    metadata = MetaData(dbID,flightID,aircraft,route,lat,lon,time,file)

    new(time,lat,lon,alt,heading,climb,speed,metadata)
  end #constructor FlightData
end #struct FlightData

struct FlightDB
  inventory::Vector{FlightData}
  archive::Vector{FlightData}
  onlineData::Vector{FlightData}
  created::Union{DateTime,tz.ZonedDateTime}
  remarks
end #struct FlightDB

struct CLay
  time::Vector{ZonedDateTime}
  lat::Vector{AbstractFloat}
  lon::Vector{AbstractFloat}

  function CLay(files::String...)
    utc = ZonedDateTime[]; lon = []; lat = []
    for file in files
      if occursin("CLay", basename(file))
        t = mat.mxcall(:hdfread,1,file,"Profile_UTC_Time")[:,2]
        utc = [utc; convertUTC.(t)]
        lon = [lon; mat.mxcall(:hdfread,1,file, "Longitude")[:,2]]
        lat = [lat; mat.mxcall(:hdfread,1,file, "Latitude")[:,2]]
      end
    end

    new(utc, lat, lon)
  end #constructor CLay
end #struct CLay

struct CPro
  time::Vector{ZonedDateTime}
  lat::Vector{AbstractFloat}
  lon::Vector{AbstractFloat}

  function CPro(files::String...)
    utc = ZonedDateTime[]; lon = []; lat = []
    for file in files
      if occursin("CPro", basename(file))
        t = mat.mxcall(:hdfread,1,file,"Profile_UTC_Time")[:,2]
        utc = [utc; convertUTC.(t)]
        lon = [lon; mat.mxcall(:hdfread,1,file, "Longitude")[:,2]]
        lat = [lat; mat.mxcall(:hdfread,1,file, "Latitude")[:,2]]
      end
    end

    new(utc, lat, lon)
  end #constructor CPro
end #struct CPro

struct SatDB
  CLay::CLay
  CPro::CPro
  created::DateTime
  remarks

  function SatDB(folder::String...; remarks=nothing)
    cl = CLay(folder...)
    cp = CPro(folder...)
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

# Use constructor for FlightData and MetaData:
# - check vector length
# - in MetaData: construct mins/maxs of area and date
# - in area distinguish between pos/neg mins/maxs
# - check format of route
# - construct MetaData within FlightData

# Outsource loop to save FlightData to separate routine loadFlightData

include("auxiliary.jl")
include("loadFlightData.jl")


"""
    loadFlightDB(DBtype, folder...)

documentation
"""
function loadFlightDB(DBtype::String, folder::Union{String, Vector{String}}...; remarks=nothing)
  # Save time of database creation
  tc = Dates.now()
  # Find database types
  if occursin('i', DBtype)  i1 = findfirst(isequal('i'), DBtype)
  else  i1 = findfirst(isequal('1'), DBtype);  end
  if occursin('a', DBtype)  i2 = findfirst(isequal('a'), DBtype)
  else  i2 = findfirst(isequal('2'), DBtype);  end
  if occursin('o', DBtype)  i3 = findfirst(isequal('o'), DBtype)
  else  i3 = findfirst(isequal('3'), DBtype);  end

  # Load databases for each type
  if !isnothing(i1)
    ifiles = String[]
    ifiles = findFiles(ifiles, folder[i1], ".csv")
    inventory = try loadInventory(ifiles)
    catch
      @warn "Flight inventory couldn't be loaded."
      FlightData[]
    end
  else inventory = FlightData[];  end
  if !isnothing(i2)
    ifiles = String[]
    ifiles = findFiles(ifiles, folder[i2], ".csv")
    archive = try loadArchive(ifiles)
    catch
      @warn "FlightAware archive couldn't be loaded."
      FlightData[]
    end
  else archive = FlightData[];  end
  if !isnothing(i3)
    ifiles = String[]
    ifiles = findFiles(ifiles, folder[i3], ".txt", ".dat")
    onlineData = try loadOnlineData(ifiles)
    catch
      @warn "FlightAware online data couldn't be loaded."
      FlightData[]
    end
  else onlineData = FlightData[];  end

  println("\ndone loading data to properties\n- inventory\n- archive\n- onlineData\n", "")

  return FlightDB(inventory, archive, onlineData, tc, remarks)
end # function loadFlightDB

end # module TrackMatcher
