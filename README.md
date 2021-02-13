TrackMatcher
============

Overview
--------

`TrackMatcher` is a Julia package to find intersections between two sets of trajectories.
The current version can match aircraft- or cloud-track (primary) data with CALIPSO satellite 
(secondary) ground tracks and store relevant data in the vicinity of the intersection.


Installation
------------

`TrackMatcher` is an unregistered Julia package, but can be installed with the
package manager. `TrackMatcher` uses the `PCHIP.jl` package developed within the
`TrackMatcher` framework. It is a pure Julia implementation of the Peacewise Cubic
Hermite Interpolating Polynomial. `PCHIP` has to be installed prior to `TrackMatcher`
as Julia currently only acknowledges registered dependencies.

```julia
julia> ]
pkg> add https://github.com/LIM-AeroCloud/PCHIP.jl.git
pkg> add https://github.com/LIM-AeroCloud/TrackMatcher.jl.git
```


Usage
-----

In essence, 3 `TrackMatcher` structs are needed to load essential primary and secondary data, 
and find intersections in the stored track data. A full overview is given in the 
[WIKI](https://github.com/LIM-AeroCloud/TrackMatcher.jl/wiki) and this README is 
only meant as a quick reminder of the most important functions.


Loading primary data
--------------------

`FlightData` of individual flights are loaded into a `FlightDB` with vectors of 
`FlightData` for different database types, currently holding

1. VOLPE AEDT inventory (`i` or `1`)
2. FlightAware archived data (`a` or `2`)
3. flightaware.com online data (`o` or `3`)

A convenience constructor for `FlightDB` exists needing only the database types 
listed in a string with the letters or numbers as indicated in the list above and 
the directories of the main database folders. Those folders are searched recursively 
for the respective data files. More than one folder path can be listed for all the database types.
The order in the list is free, but the order of folders must correspond to the order
of dataset identifiers in `DBtype`:

```julia
FlightDB(DBtype::String, folder::Union{String, Vector{String}}...; kwargs)
```

A similar convenience constructor exists for `CloudTrack`s. As only one database type
exists, only the directories are needed:

```julia
CloudDB(folders::String...; kwargs)
```

### kwargs

#### FlightDB and CloudDB
- `Float::DataType=Float32`: Set the default precision of floating point numbers for flight data
- `remarks=nothing`: any data or comments that can be attached to  `DBMetadata`

#### FlightDB only
- `altmin::Int=5000`: minimum altitude threshold for which to consider flight data
- `odelim::Union{Nothing,Char,String}=nothing`: specify the column delimiter in the text files of the online data


Loading CALIOP data from the CALIPSO satellite
----------------------------------------------

CALIPSO positions and overpass times together with a file index of the corresponding
granule hdf file are stored in the `data` field of `SatData`. Only one of the `type`s
cloud profile (`CPro`) or cloud layer (`CLay`) data can be used to construct `SatData`.
The `metadata` holds a `Dict` with the `fileindex` pointing to a file name (including
the absolute folder path). 
__File names/directories must not be changed in order for _TrackMatcher_ to find intersections correctly.__
__It is currently only possible to find intersections on the same system, where the data was loaded.__  
Further information in the `metadata` include the `type` of the satellite data,
the `date` range of the data, the time the database was `created`, the `loadtime`,
and any `remarks` as additional data or comments.

`SatData` can be instatiated, by giving any number of folder strings and any remarks
using the keyword `remarks`. The `folders` are scanned recursively for any hdf file
and the `type` of the satellite data is determined by keywords `CLay` or `CPro` in
the folder/file names. If both types exist in the `folders`, the data type is determined
by the majority in the first 50 file names. Alternatively, the sat data type can
be forced with the keyword `type` set to a `Symbol` `:CPro` or `:CLay`.

```julia
SatData(folders::String...; kwargs...)
```

### kwargs
- `Float::DataType=Float32`: Set the default precision of floating point numbers for satellite data
- `type::Symbol=:undef`: Set the satellite data type to layer (`CLay`) or profile (`CPro`) data
- `remarks=nothing`: any data or comments that can be attached to the metadata of `SatData`

---
> :information_source: **NOTE**
>
> `SatData` is designed to use CALIPSO data provided by the [AERIS/ICARE Data and Services Centre](http://www.icare.univ-lille1.fr/). 
> For the best performance, you should use the same file/folder format as used by ICARE. 
> In particular, Cloud layer files must include the keyword `CLay` in the file name
> and cloud profile data files the keyword `CPro`.
---


Finding intersections in the trajectories of the flight and satellite data
--------------------------------------------------------------------------

Intersections and corresponding accuracies and flight/satellite data in the vicinity of the intersection are stored in the `Intersection` struct.

A convenience constructor exists for automatic calculation of the intersections 
from the `FlightDB` and `SatData` with parameters controlling these calculations. 
Additionally, it can be specified by the keyword `savesecondsattype` whether the 
corresponding satellite data type of the `CLay` or `CPro` data stored in `SatData`
should be saved as well. 
__For this feature to work, folder and file names of `Clay`/`CPro` data must be identical__
__except for the keywords `CLay`/`CPro` swapped.__

Find intersections by instatiating the `Intersection` struct with:

```julia
Intersection(flights::FlightDB, sat::SatData, savesecondsattype::Bool=false; kwargs...)
Intersection(cloud::CloudDB, sat::SatData, savesecondsattype::Bool=false; kwargs...)
```

### kwargs

- `maxtimediff::Int=30`: maximum delay at intersection between aircraft/satellite overpass
- `primspan::Int=0`: number of primary (flight or cloud) data points saved before and after the closest measurement to the intersection
- `secspan::Int=15`: number of secondary (satellite) data points saved before and after the closest measurement to the intersection
- `lidarrange::Tuple{Real,Real}=(15_000,-Inf)`: top/bottom bounds of the lidar column data, between which
  data is stored; use `(Inf, -Inf)` to store the whole column
- `stepwidth::Float64=1000`: stepwidth in degrees (at the equator) used for the 
  interpolation of flight and satellite tracks
- `Xradius::Real=20_000`: Radius in meters, in which multiple intersection finds are
  assumed to correspond to the same intersection and only the intersection with the
  minimum delay between flight and sat overpass is saved
- `expdist::Real=Inf`: maximum threshold for the distance of a found intersection to the nearest measured track point
- `Float::DataType=Float32`: default precision of floating point numbers for intersection data
- `remarks=nothing`: any data or comments that can be attached to the metadata of `Intersection`
