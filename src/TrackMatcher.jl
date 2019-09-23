module TrackMatcher

# Track changes during development
using Revise

# Import Julia packages
import CSVFiles; const csv = CSVFiles
import DataFrames; const df = DataFrames
import TimeZones; const tz = TimeZones
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


pchip = itp.PchipInterpolator([1,2,3,4,5,6,7], [0,0,0,3,6,5.7,6])
spline = spl.Spline1D([1,2,3,4,5,6,7], [0,0,0,3,6,5.7,6], bc="extrapolate")
plt.clf()
plt.scatter([1,2,3,4,5,6,7], [0,0,0,3,6,5.7,6], marker="x", color="k", label="data")
plt.plot(0:0.01:8,[pchip(i) for i in 0:0.01:8], label="pchip")
plt.plot(0:0.01:8,spline(0:0.01:8), label="spline")
plt.legend(); plt.grid(ls=":")
plt.gcf()

# Define own structs
struct flightDB
  inventory::Vector{flightData}
  archive::Vector{flightData}
  onlineData::Vector{flightData}
  created::Union{DateTime,tz.ZonedDateTime}
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

# Use constructor for flightData and MetaData:
# - check vector length
# - in MetaData: construct mins/maxs of area and date
# - in area distinguish between pos/neg mins/maxs
# - check format of route
# - construct MetaData within flightData


struct MetaData
  dbID::Union{Int32,String}
  flightID::String
  route::String
  area::Vector{Float32}
  date::Vector{Union{DateTime,tz.ZonedDateTime}}
end

# Load data
flights = DataFrame(csv.load("data/flightinventory/1_1_2012_SEGMENT.csv", skiplines_begin=3,
  header_exists=false, colnames=[:id, :seg, :lat, :lon, :alt, :month, :day, :year,
  :hour, :min, :sec, :emiss, :temp, :press, :rH, :speed, :segtime, :segdist, :thrust,
  :weight, :fuel, :co, :hc, :nox, :pmnv, :pmso, :pmfo, :co2, :h2o, :sox],
  colparsers=Dict(:id => Int32, :lat => Float32, :lon =>Float32, :alt => Float32, :speed => Float32)))
df.deleterows!(flights, length(flights.id))
flights.time = ZonedDateTime.(flights.year, flights.month, flights.day,
  flights.hour, flights.min, flights.sec, tz.tz"UTC")

# Initialise
global iStart = 1
inventory = flightData[];
# Loop over all data points
for i = 1:length(flights.time)-1
  if flights.id[i] ≠ flights.id[i+1]
    # When flight ID changes save data as flightData and set start index to next dataset
    push!(inventory, flightData(flights.time[iStart:i], flights.lat[iStart:i],
      flights.lon[iStart:i], flights.alt[iStart:i], [missing for j=iStart:i],
      [missing for j=iStart:i], flights.speed[iStart:i]))
    global iStart = i+1
  end
end
# Save last flight of the dataset
i = length(flights.time)
push!(inventory, flightData(flights.time[iStart:i], flights.lat[iStart:i],
  flights.lon[iStart:i], flights.alt[iStart:i], [missing for j=iStart:i],
  [missing for j=iStart:i], flights.speed[iStart:i]))

end # module
