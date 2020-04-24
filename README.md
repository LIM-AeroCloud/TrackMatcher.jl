TrackMatcher
============

Overview
--------

`TrackMatcher` is a Julia package to find intersections between airplane and CALIPSO satellite flight tracks and store relevant data in the vicinity of the intersection.


Installation
------------

`TrackMatcher` is an unregistered Julia package, but is easily installed with

```julia
julia> ]
pkg> add https://github.com/pb866/TrackMatcher.jl.git
pkg> instantiate
```


Usage
-----

In essence, 3 `TrackMatcher` structs are needed to load essential flight and satellite data, and find intersection in the stored track data. A full overview is given in the [WIKI](https://github.com/pb866/TrackMatcher.jl/wiki) and this README is only meant as a quick reminder of the most important functions.


Loading flight data
-------------------

`FlightData` of individual flights are loaded into a `FlightDB` with vectors of `FlightData` of different database types, currently holding

1. VOLPE AEDT inventory (`i` or `1`)
2. FlightAware archived data (`a` or `2`)
3. flightaware.com online data (`o` or `3`)

A convenience constructor for `FlightDB` exists needing only the database type listed in a string with the letters or numbers as indicated in the list above and the directories of the main database folders. Those folders are searched recursively for the respective data files. More than one folder path can be listed for all the database types.

```julia
FlightDB(DBtype::String, folder::Union{String, Vector{String}}...; kwargs)
```

### kwargs
- `altmin::Int=15_000`: minimum altitude threshold for which to consider flight data
- `remarks=nothing`: any data or comments that can be attached to the metadata of `FlightDB`
- `odelim::Union{Nothing,Char,String}=nothing`: specify the column delimiter in the text files of the online data


Loading CALIOP data from the CALIPSO satellite
----------------------------------------------

Satellite cloud layer data (`CLay`) or cloud profile data (`CPro`) can be loaded to the `SatDB` with a convenience constructor giving one or more folder paths with CALIOP data of version `4.x`. The original file names should be kept, as the constructor scan for `CLay` and `CPro` in the file names to automatically assign the data to the correct database.
Any comments or data can be attached with the keyword argument `remarks` to the metadata of `SatDB`.

```julia
SatDB(folders::String...; remarks=nothing)
```


Finding intersections in the trajectories of the flight and satellite data
--------------------------------------------------------------------------

Intersections and corresponding accuracies and flight/satellite data in the vicinity of the intersection are stored in the `Intersection` data type.

A convenience constructor exists, automatically calculating the intersections from the `FlightDB` and `SatDB` data with parameters that control these calculations. Additionally, it can be specified which `sattype` (`CLay` (default) or `CPro`) should be preferred for finding the intersections.


```julia
Intersection(flights::FlightDB, sat::SatDB, sattype::Symbol=:CLay; kwargs)
```

### kwargs

- `maxtimediff::Int=30`: maximum delay at intersection between aircraft/satellite
- `flightspan::Int=0`: number of flight data points saved before and after the closest measurement to the intersection
- `satspan::Int=15`: number of satellite data points saved before and after the closest measurement to the intersection
- `stepwidth::Float64=0.01`: stepwidth in degrees used for the interpolation of flight and satellite tracks
- `remarks=nothing`: any data or comments that can be attached to the metadata of `Intersection`

<!-- - `Xradius::Real=5000` -->
