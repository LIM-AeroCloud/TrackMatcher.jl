"""
# Module TrackMatcher

Find intersections between different sets of trajectories. The module is aimed to find
intersections between aircraft and satellite tracks with initial amendments for cloud
tracks and room for further applications in geosciences.

# Type tree

```
DataSet──Data
├─MeasuredSet──MeasuredData
│ ├─PrimarySet
│ │ ├─PrimaryMetadata
│ │ ├─FlightSet
│ │ ├─CloudSet
│ │ └─PrimaryTrack
│ │   ├─FlightTrack──FlightData
│ │   │ └─FlightMetadata
│ │   └─CloudTrack──CloudData
│ │     └─CloudMetadata
│ ├─SecondarySet
│ │ ├─SecondaryMetadata
│ │ ├─SatSet
│ │ └─SecondaryTrack
│ │   └─SatTrack──SatData
│ │     └─SatMetadata
│ └─ObservationSet
│   ├─CPro
│   ├─CLay
│   ├─APro
│   └─ALay
└─ComputedSet
      └─Intersection──XData
        └─XMetadata
```

## Exported structs

### DataSet

- top level abstract type that combines all kinds of TrackMatcher data
- alias constructor for `Data` to load data of 2 different sets (primary and secondary data)
  and calculate intersections in 1 step


### Data

- concrete type to store all types of measured and calculated TrackMatcher data
- constructor to load all primary and secondary data and calculate intersections
  between all matching primary/secondary datasets


### MeasuredSet

- abstract type for all data set with measured track data and/or observations


### PrimarySet

- abstract type for primary data with individual trajectories of a common type
  (currently `FlightSet` or `CloudSet`)
- alias constructors for `FlightSet` and `CloudSet`


### PrimaryMetadata

- Metadata of all primary data sets


### FlightSet

- concrete type to store flight track data and other relevant aircraft related data
  from 3 different inventories:
  - `inventory`: VOLPE AEDT inventory
  - `archive`: commercially available database by FlightAware
  - `onlineData`: free online data by FlightAware


### CloudSet

- concrete type to store individual cloud tracks and related observations


### PrimaryTrack

- abstract type for structs storing individual (primary) track data


### FlightTrack

- abstract type for individual flight tracks
- alias constructor for `FlightData` to load individual flight track data


### FlightData

- concrete type to store data of a single flight within `FlightSet`


### FlightMetadata

- concrete type with metadata to individual `FlightData`


### CloudTrack

- abstract type for individual cloud tracks
- alias constructor for `CloudData` to load individual cloud track data


### CloudData

- concrete type to store data of a single cloud track within `CloudSet`


### CloudMetadata

- concrete type with metadata to individual `CloudData`

### SecondarySet

- abstract type for all secondary (satellite) data

### SatSet

- concrete type to store satellite track data

### SecondaryMetadata

- concrete type with metadata to `SatSet`

### SecondaryTrack

- abstract type for structs storing secondary (satellite) track data


### SatTrack

- abstract type for satellite track data
- alias constructor for `SatData` to load satellite track data in a continuous format


### SatData

- concrete type to store basic satellite track data as continuous track


### SatMetadata

- concrete type with metadata to `SatSet`


### ObservationSet

- abstract type for struct holding data of satellite observations


### CPro

- concrete type holding CALIOP cloud profile data


### CLay

- concrete type holding CALIOP cloud layer data


### APro

- concrete type holding CALIOP aerosol profile data


### ALay

- concrete type holding CALIOP aerosol layer data


### ComputedSet

- abstract type for data of computed intersections


### Intersection

- abstract type for data of calculated intersections and related observations for
  specified sets of primary and secondary trajectories
- alias constructor for `XData`


### XData

- concrete type holding calculated intersections and related observations for
  specified sets of primary and secondary trajectories


### XMetadata

- concrete type for intersection metadata
"""
module TrackMatcher

## Import Julia packages
import DataFrames as df
import DataStructures as ds
import StructArrays as sa
import HDF5 as h5
import CSV
import Dates
import TimeZones as tz
import Distances as dist
import MAT
import IntervalRootFinding as root
import Statistics as stats
import ProgressMeter as pm
import Logging as logg
import IntervalArithmetic as intar

# Import structs and functions from packages
import PCHIP: Polynomial, pchip, interpolate
import DataFrames: DataFrame, AbstractDataFrame
import StructArrays: StructArray
import DataStructures: OrderedDict, SortedSet
import Dates: AbstractDateTime, DateTime, Date, Time
import TimeZones: ZonedDateTime
import IntervalArithmetic.Symbols: (..)

# Define Logger with log level
logger = logg.ConsoleLogger(stdout, logg.Info)
logg.global_logger(logger)


# Disable progress bars in test runs or when explicitly requested.
progress_enabled() = haskey(ENV, "TRACKMATCHER_PROGRESS") ? ENV["TRACKMATCHER_PROGRESS"] == "true" : true


## Define time zones for FlightAware online data
zonedict = Dict{String,tz.TimeZone}(
    "CET" => tz.tz"+0100",
    "CEST" => tz.tz"+0200"
)


## Export types and constructors
export DataSet, Data, MeasuredSet, MeasuredData, ComputedSet, PrimarySet, SecondarySet,
       ObservationSet, FlightSet, CloudSet, SatSet, PrimaryTrack, SecondaryTrack,
       FlightTrack, CloudTrack, SatTrack, FlightData, CloudData, SatData,
       CLay, CPro, Intersection, XData, #APro, ALay,
       FlightMetadata, CloudMetadata, PrimaryMetadata, SecondaryMetadata, XMetadata,
       AerosolType, CloudType, SkyCondition, invalid, clear, cloud, aerosol, stratospheric,
       surface, subsurface, no_signal, undetermined, marine, dust, polluted, remote, polluted_dust,
       smoke, other, low_trasparent, low_opaque, transition_sc, cu, ac, as, ci, cb,
       unknown, clay, cpro, apro, alay


## Source files
include("types.jl")           # type tree and enums
include("flighttypes.jl")     # concrete types/constructors for flight data/datasets
include("flightdata.jl")      # functions related to loading flight databases/datasets
include("cloudtypes.jl")      # concrete types/constructors for cloud data/datasets
include("clouddata.jl")       # functions related to loading cloud track databases/datasets
include("sattypes.jl")        # concrete types/constructors for secondary sat track data and observations
include("datachecks.jl")      # helper functions for data checks
include("conversions.jl")     # helper functions for time/unit conversions
include("outputtypes.jl")     # concrete types/constructors for intersections and combined datasets
include("observations.jl")    # concrete types/constructors for satellite observations
include("match.jl")           # functions related to finding track intersections
include("dataprocessing.jl")  # helper functions for data processing
include("lidar.jl")           # functions related to processing CALIOP lidar data
include("base_overloads.jl")  # overloads for Base equality/approximation functions for TrackMatcher types

end # module TrackMatcher
