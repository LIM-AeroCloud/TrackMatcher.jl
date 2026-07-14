# Release notes

## [UNRELEASED]

### Added

- Setup tests ([#55])
- Improve data checks, error handling and logging for loading webdata ([#67])
- Add `HDF5` and `StructArrays` as dependency ([#51], [#58])
- Add method `checklimits` to check array fields of `SatData` are within expected limits ([#51])
- Setup code coverage for _TrackMatcher_ and add badges for the CI status and code coverage
  to README ([#64])
- Improved error handling for lidar heights retrieval ([#75])
  - exact calculation of last fine-grained level height (~29.94m) instead of 30m approximation
  - Set `i.f.h30m` to `1` instead of `0` for a chosen height range below 8.3km to allow the calculation
    of all fine-grained levels instead of throwing an error
  - Add error handling to `get_lidarcolumn` to throw a descriptive error for wrong user input

### Changed

- Update compatibility of dependencies
- Rename `convertFloats!` and `convertUTC` to `convert_floats!` and `convert_utc`, respectively,
  to be in line with Julia conventions
- Refactor `convert_floats!` without type piracy
- Rename fields and variables `useLON` to `use_lon`
- Refactor satellite types and constructors ([#51], [#58])
  - Remove dataframe layer in `SatData` and store `time`, `lat`, `lon` directly in `StructArray`
    fields
  - Remove dataframe layer in observation data `CLay` and `CPro` and store data in separate fields
  - Rename field `remarks` in `SecondaryMetadata` to `attachments` to signal more flexibility
  - Refactor internal `findfiles!` function using `scandir` in connection with `sat_datafiles!`
  - Add additional data checks on instantiation
- Rename `checkbounds!` to `checklimits!`/`checklimits` to avoid confusion with `checkbounds`
  from Base ([#51])
- Optimise load times for flight data by preferring column-wise operations over row-wise operations ([#58])
  and making use of DataFrames functions for grouping and transforming data ([#58])
- Refactor primary flight and cloud data and datasets to use StructArrays and optimize load times ([#58])
  - Rename `dbID` in `FlightData` to `id` and `flightID` to `flight_num` for more consistency
    in the naming scheme and memory usage ([#58])
- Refine the `XData` struct and the calculation of intersection points ([#51], [#58])
- Refine file handling and save root paths and file paths/names separately ([#58])
- Use look-up dictionaries in the set metadata and reduced Strings or `UInt16` indices in the
  track metadata to reduce struct sizes ([#58])
- Use `Enum` instead of `Symbol` to describe the atmospheric conditions derived from the
  Feature Classification Flag (FCF) [#58]
- Disable validity checks during construction of observation data until less restrictive.
  Future re-implementations should only disregard single outliers not the whole column of data
  with outliers. ([#58])
- Simplify retrieval of lidar column data ([#75])
  - refactor lidar.jl
  - use more generous boundaries including the first value outside each threshold for more exact
    calculations at the edges (technically **breaking**)
  - revised return value of `get_lidarheights`
- Use flight number for error handling in `CLay` method of `atmosphericinfo` ([#75])

### Removed

- Remove MATLAB dependency ([#51])
- Ignore Manifest.toml, this should be auto-generated on each system
- Remove constructors for `Float16`, `Float32`, and `Float64` taking `missing` as input to
  avoid type piracy

### Fixed

- Fix constructor for empty `CloudSet` ([#71])
- Ensure empty `Cloudset` is returned, if no cloud data is found in the given path(s) ([#71])
- Fixed errors in UTC time conversion by rounding the converted seconds of the day from the 
  fraction of the day to 3 digits
- Fix and issue that did not load flight webdata with the given time zone but always used the
  local time zone ([#67])
- `earthradius` now always return `Float64` instead of the input precision to ensure correct
  results and no overflow for `Float16`
- Fixed compilation warnings by removing duplicate constructor for empty `FlightSet`

## [v0.5.4]

### Added

- Add issue templates for more guided issue creation

### Changed

- Updated input format of cloud data to use `centrLatLon` instead of `centrLonLat`, use 
  `"filtered_trajectory"` as default struct name instead of `"cloud"`

### Fixed

- Fixed UTC time conversion for dates pre-2010

## [v0.5.3] - 2022-02-20

### Added

- Add raw data and scripts for result plots in paper in folder paperplots

### Changed

- Update WIKI
- Minor update to CSV.read without multithreading

## [v0.5.2] – 2021-11-28

### Fixed

- Fix unmodified constructor of XData

## [v0.5.1] - 2021-10-30

### Changed

- Update to most recent dependencies
- Update constructor for combined TrackMatcher processes (MeasuredSet and DataSet) especially to new
  FlightSet fields

### Fixed

- Bug fixes and performance improvements

## [v0.5.0] - 2021-09-09

### Added

- Introduce type tree
- New overall constructors MeasuredData to load all primary and secondary data at once and Data
  to load all input data and calculate intersections in one go

### Changed

- Restructure satellite data and store data as individual granules
- Revised naming scheme of variables and struct fields

### Fixed

- Performance improvements and fixes

## [v0.4.0] - 2021-02-13

### Added

- Add cloud tracking data; define primary datasets as aircraft or cloud data and secondary
  datasets as satellite data

### Changed

- Use a function to obtain the distance between the primary and secondary trajectory
- Use _IntervalRootFindling.jl_ to obtain roots in the distance function for the location of
  intersections

## [v0.3.0] - 2020-10-26

### Added

- Add experimental CPro data
- Allow a default precision of floating point numbers and remove `AbstractFloat`s
- Allow different _FlightAware_ archive versions

### Changed

- Replace MATLAB PCHIP version by a native Julia PCHIP package
- Use `haversine` function to calculate great circle distance between coordinate pairs rather
  than `distance` function from _Geodesy_ package
- Refine time interpolation using linear interpolation between closest measured track points
- Revise `tolerance` parameters in Intersection
- Require at least Julia 1.5

## [v0.2.0] - 2020-07-18

### Added

- Introduction of `SatData` with `time`/`lat`/`lon`, and a `fileindex` pointing to additional data
- Include `CLay`/`CPro` with additional measurements directly in `Intersection.tracked`
- Store atmospheric conditions in the vicinity of the intersection (at flight level) as `feature`

### Changed

- Complete package overhaul
- Updated track-matching algorithm
- Revised kwargs

### Removed

- `SatDB` struct

## [v0.1.2] - 2020-02-24

### Added

- Julia 1.0.x compatible code
- MATLAB 2016 compatible code
- New README

## [v0.1.1] - 2020-02-06

### Removed

- Comment out _Revise_ package needed during development, but not part of the
  `Project.toml`/`Manifest.toml`

### Fixed

- Correct type of `precision` being written as `prec` in function `find_intersections`

## [v0.1.0] - 2020-02-06

### Added

- Preliminary version that finds intersections between satellite and aircraft flight trajectories.
  Filtering for cirrus clouds at flight level not yet implemented.
