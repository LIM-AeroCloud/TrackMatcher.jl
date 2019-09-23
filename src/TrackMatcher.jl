module TrackMatcher

# Track changes during development
using Revise

# Import Julia packages
import CSVFiles; const csv = CSVFiles
import DataFrames; const df = DataFrames
import TimeZones; const tz = TimeZones
import Dates
import Dierckx; const spl = Dierckx
import PyCall; const py = PyCall
# import Conda
import PyPlot; const plt = PyPlot
# Import structs from packages
import DataFrames.DataFrame
import Dates.DateTime
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
  dbID::Union{Int32,String}
  flightID::String
  route::String
  area::Vector{Float32}
  date::Vector{Union{DateTime,tz.ZonedDateTime}}
  filename::String
end

struct flightData
  time::Vector{ZonedDateTime}
  lat::Vector{Union{Missing,Float32}}
  lon::Vector{Union{Missing,Float32}}
  alt::Vector{Union{Missing,Float32}}
  heading::Vector{Union{Missing,Int32}}
  climb::Vector{Union{Missing,Int32}}
  speed::Vector{Union{Missing,Float32}}
  # metadata::MetaData
end #flightData

# Define own structs
struct flightDB
  inventory::Vector{flightData}
  archive::Vector{flightData}
  onlineData::Vector{flightData}
  created::Union{DateTime,tz.ZonedDateTime}
  remarks
end

# Use constructor for flightData and MetaData:
# - check vector length
# - in MetaData: construct mins/maxs of area and date
# - in area distinguish between pos/neg mins/maxs
# - check format of route
# - construct MetaData within flightData

include("inventory.jl")


"""
    loadDB(DBtype, folder...)

documentation
"""
function loadDB(DBtype::String, folder::Union{String, Vector{String}}...; remarks=[])
  # Save time of database creation
  tc = DateTime(Dates.now())
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
    inventory = loadInventory(ifiles)
  else inventory = flightData[];  end
  if !isnothing(i2)
    ifiles = String[]
    ifiles = findtextfiles(ifiles, folder[i2])
    archive = loadArchive(ifiles)
  else archive = flightData[];  end
  if !isnothing(i3)
    ifiles = String[]
    ifiles = findcsv(ifiles, folder[i3])
    onlineData = loadOnlineData(ifiles)
  else onlineData = flightData[];  end

  println("\ndone loading data to properties\n- inventory\n- archive\n- onlineData\n", "")

  return flightDB(inventory, archive, onlineData, tc, remarks)
end # function loadDB

end # module TrackMatcher
