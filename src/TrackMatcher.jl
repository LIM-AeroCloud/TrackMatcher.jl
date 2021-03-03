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
- `SetMetadata` stores metadata for the primary database (`FLightDB` or `CloudDB`)
- `FlightTrack` stores data of a single flight in `FlightDB`
- `FlightMetadata` holds metadata to every flight
- `SatData` stores CALIPSO position and time and a file index for the granule of each data line
- `CLay` CALIPSO cloud layer data
- `CPro` CALIPSO cloud profile data
- `SatMetadata` stores metadata of the CALIPSO data
- `Intersection`: positions of intersections between aircraft and satellite trajectories,
  together with time difference between crossings; and `FlightTrack` and `CPro`/`CLay` in the
  vicinity of the intersection as well as information about the accuracy of the data
- `XMetadata` stores metadata for the `Intersection` data
"""
module TrackMatcher

## Import Julia packages
import DataFrames; const df = DataFrames
import DataStructures; const ds = DataStructures
import CSV
import Dates
import TimeZones; const tz = TimeZones
import Distances; const dist = Distances
import MATLAB; const mat = MATLAB
import IntervalRootFinding; root = IntervalRootFinding
import IntervalArithmetic...
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


## Define time zones for FlightAware online data
zonedict = ds.DefaultDict{String,tz.TimeZone}(tz.localzone())
zonedict["_CET_"] = tz.tz"+0100"
zonedict["_CEST_"] = tz.tz"+0200"


## Define type tree of abstract types
abstract type DataSet{T<:AbstractFloat} end
abstract type MeasuredSet{T} <: DataSet{T} end
abstract type ComputedSet{T} <: DataSet{T} end
abstract type PrimarySet{T} <: MeasuredSet{T} end
abstract type ObservationSet{T} <: MeasuredSet{T} end
abstract type SecondaryTrack{T} <: MeasuredSet{T} end
abstract type PrimaryTrack{T} <: PrimarySet{T} end
abstract type FlightTrack{T} <: PrimaryTrack{T} end
abstract type CloudTrack{T} <: PrimaryTrack{T} end
abstract type SatTrack{T} <: SecondaryTrack{T} end
# abstract type Intersection{T} <: ComputedSet{T} end


## Export types and constructors
export DataSet, MeasuredSet, ComputedSet, PrimarySet, ObservationSet,
       PrimaryTrack, SecondaryTrack, FlightTrack, CloudTrack, FlightData, #CloudData,
       SecondaryTrack, SatData, CLay, CPro, Intersection, #APro, ALay, XData,
       FlightMetadata, SatMetadata, SetMetadata, XMetadata


## Import functions for Julia include files
include("inputtypes.jl")      # structs of concrete types at the end of the type tree
include("outputtypes.jl")     # structs of abstract types with constructors for concrete types
include("datachecks.jl")      # helper functions for data checks
include("dataprocessing.jl")  # helper functions for data processing
include("conversions.jl")     # helper functions for time/unit conversions
include("lidar.jl")           # functions related to processing CALIOP lidar data
include("flightdata.jl")      # functions related to loading flight databases/datasets
include("clouddata.jl")       # functions related to loading cloud track databases/datasets
include("match.jl")           # functions related to finding track intersections

end # module TrackMatcher
