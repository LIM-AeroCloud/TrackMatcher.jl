module TrackMatcher

# Track changes during development
# using Revise

# Import Julia packages
import CSV
import DataFrames; const df = DataFrames
import TimeZones; const tz = TimeZones
import Dates
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
  end
end

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
  end
end #FlightData

# Define own structs
struct FlightDB
  inventory::Vector{FlightData}
  archive::Vector{FlightData}
  onlineData::Vector{FlightData}
  created::Union{DateTime,tz.ZonedDateTime}
  remarks
end

export loadDB,
       loadInventory,
       loadArchive,
       loadOnlineData,
       FlightDB,
       FlightData,
       MetaData

# Use constructor for FlightData and MetaData:
# - check vector length
# - in MetaData: construct mins/maxs of area and date
# - in area distinguish between pos/neg mins/maxs
# - check format of route
# - construct MetaData within FlightData

# Outsource loop to save FlightData to separate routine loadFlightData


function checklength(vect, ref)
  len = length(ref) - length(vect)
  if len > 0
    @warn string("$(len) entries missing in vector compared to reference. ",
      "Missing entries are filled with `missing` at the end of the vector.")
    vect = [vect; [missing for i = 1:len]]
  elseif len < 0
    @warn string("$(-len) additional entries found in vector compared to reference. ",
      "Additional entries at the end of the vector are ignored.")
    vect = vect[1:length(ref)]
  end

  return vect
end

include("inventory.jl")
include("archive.jl")
include("onlineData.jl")

"""
    loadDB(DBtype, folder...)

documentation
"""
function loadDB(DBtype::String, folder::Union{String, Vector{String}}...; remarks=[])
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
    ifiles = findcsv(ifiles, folder[i1])
    inventory = try loadInventory(ifiles)
    catch
      @warn "Flight inventory couldn't be loaded."
      FlightData[]
    end
  else inventory = FlightData[];  end
  if !isnothing(i2)
    ifiles = String[]
    ifiles = findcsv(ifiles, folder[i2])
    archive = try loadArchive(ifiles)
    catch
      @warn "FlightAware archive couldn't be loaded."
      FlightData[]
    end
  else archive = FlightData[];  end
  if !isnothing(i3)
    ifiles = String[]
    ifiles = findtextfiles(ifiles, folder[i3])
    onlineData = try loadOnlineData(ifiles)
    catch
      @warn "FlightAware online data couldn't be loaded."
      FlightData[]
    end
  else onlineData = FlightData[];  end

  println("\ndone loading data to properties\n- inventory\n- archive\n- onlineData\n", "")

  return FlightDB(inventory, archive, onlineData, tc, remarks)
end # function loadDB

end # module TrackMatcher
