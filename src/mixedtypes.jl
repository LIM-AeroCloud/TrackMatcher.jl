## Define structs related to flight data

"""
# struct FlightTrack{T<:AbstractFloat}

Aircraft data with fields
- `data::DataFrame`
- `metadata::FlightMetadata`

The `DataFrame` of `data` has columns in the following order with the respective types:

- `time::Vector{DateTime}`                        (time stamp of measurement)
- `lat::Vector{<:Union{Missing,<:AbstractFloat}}`   (latitude in deg)
- `lon::Vector{<:Union{Missing,<:AbstractFloat}}`   (longitude in deg)
- `alt::Vector{<:Union{Missing,<:AbstractFloat}}`   (altitude in meters)
- `heading::Vector{<:Union{Missing,Int}}`         (heading/direction of the aircraft in deg)
- `climb::Vector{<:Union{Missing,Int}}`           (climbing (positive)/sinking (negative values) rate in m/s)
- `speed::Vector{<:Union{Missing,<:AbstractFloat}}` (velocity in m/s)

By default, all AbstractFloat are set to `Float32`, but can be set to any other
precision in `FlightDB` by the kwarg `Float`.


# Instantiation

    FlightTrack(
      FlightTrack::DataFrame,
      dbID::Union{Int,AbstractString},
      flightID::Union{Missing,AbstractString},
      aircraft::Union{Missing,AbstractString},
      route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
      flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,<:AbstractFloat,<:AbstractFloat}}}},
      useLON::Bool,
      source::String,
      file::AbstractString
    ) -> struct FlightTrack

Construct `FlightTrack` from the `data` `DataFrame` and additonal metainformation
`dbID`, `flightID`, `aircraft` type, `route`, `flex` points in the flight trajectory,
a flag whether longitudes are used as x data in the track (`useLON`), the database
`source`, and `file` names.
For `data`, `time` can be passed `ZonedDateTime`, which will be converted to `UTC`
standard time, or as UTC `DateTime`.

Or construct by directly handing over every field:

    FlightTrack(data::DataFrame, metadata::FlightMetadata)

Checks exist that the order, names, and types of the `data` `DataFrame` are correct.
"""
struct FlightData{T} <: FlightTrack{T}
  data::DataFrame
  metadata::FlightMetadata

  """ Unmodified constructor for `FlightData` with basic checks for correct `data`"""
  function FlightData{T}(data::DataFrame, metadata::FlightMetadata{T}) where T<:AbstractFloat
    # Column checks and warnings
    standardnames = ["time", "lat", "lon", "alt", "heading", "climb", "speed"]
    standardtypes = [Union{DateTime,Vector{DateTime}},
      Vector{T}, Vector{T},
      Vector{<:Union{Missing,<:T}},
      Vector{<:Union{Missing,Int}},
      Vector{<:Union{Missing,<:T}},
      Vector{<:Union{Missing,<:T}}]
    bounds = (:lat => (-90, 90), :lon => (-180, 180), :alt => (0,Inf),
      :heading => (0, 360), :speed => (0, Inf))
    checkcols!(data, standardnames, standardtypes, bounds,
      metadata.source, metadata.dbID)
    new{T}(data,metadata)
  end #constructor 1 FlightData
end #struct FlightData


""" Modified constructor with variable checks and some automated calculation of fields """
function FlightTrack{T}(
  track::DataFrame,
  dbID::Union{Int,AbstractString},
  flightID::Union{Missing,AbstractString},
  aircraft::Union{Missing,AbstractString},
  route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{<:AbstractString,<:AbstractString}}},
  flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}},
  useLON::Bool,
  source::String,
  file::AbstractString
) where T
  # Check dataframe columns of flight data; fill missing columns with missing values
  t = getproperty(track, :time)
  t = t isa Vector{ZonedDateTime} ? [zt.utc_datetime for zt in t] : t
  lat = getproperty(track, :lat)
  lon = getproperty(track, :lon)
  alt = getproperty(track, :alt)
  heading = df.hasproperty(track, :heading) ? df.getproperty(track, :heading) :
    [missing for i in t]
  climb = df.hasproperty(track, :climb) ? df.getproperty(track, :climb) :
    [missing for i in t]
  speed = df.hasproperty(track, :speed) ? df.getproperty(track, :speed) :
    [missing for i in t]
  # T = promote_type(eltype(lon), eltype(lat))
  metadata = FlightMetadata{T}(dbID,flightID,route,aircraft,t,lat,lon,useLON,flex,source,file)

  # Instatiate new FlightTrack
  return FlightData{T}(DataFrame(time=t,lat=lat,lon=lon,alt=alt,heading=heading,climb=climb,speed=speed),metadata)
end #constructor 2 FlightTrack


""" External constructor for emtpy FlightTrack struct """
FlightTrack{T}() where T = FlightData(DataFrame(time = DateTime[],
  lat = T[], lon = T[], alt=T[], heading = Int[], climb = T[],
  speed = T[]), FlightMetadata{T}())


"""
# struct FlightDB

Database for aircraft data of different database types with fields:
- `inventory::Vector{FlightTrack}`
- `archive::Vector{FlightTrack}`
- `onlineData::Vector{FlightTrack}`
- `metadata`: `DBMetatadata` with information about
  - `altmin::Real`: Minimum altitude threshold for which flight data is considered
  - `date`: date range of the dataset
  - `created`: time of creation
  - `loadtime`
  - `remarks` (any additional data or comments)

## inventory
Flight data from csv files.

## archive
Commercial flight data by FlightAware from a csv file.

## onlineData
Online data from the FlightAware website copied to whitespace-separated files.

## metadata
Immutable struct `DBMetadata`.


# Instantiation

Instatiate by giving a String with identifiers of the `DBtype` and an equal number
of `folder` paths as characters in the `DBtype` `String`. Optionally add a minimum
altitude threshold for the data (default = `15000`) and any remarks
(comments or additional data). Define the delimiter in the input files of the
online data with the keyword `odelim`. Use any character or string as delimiter.
By default (`odelim=nothing`), auto-detection is used.

    FlightDB(DBtype::String, folder::Union{String, Vector{String}}...;
      Float::DataType=Float32, altmin::Real=5_000, remarks=nothing,
      odelim::Union{Nothing,Char,String}=nothing)

`DBtype` can be identified with:
- `1` or `i`: VOLPE AEDT inventory
- `2` or `a`: FlightAware archived data (commercially available)
- `3` or `o`: flightaware.com online data

By default, all values are read in as `Float32`, but can be set to any other
precision by the `Float` kwarg; `altmin` set the minimum threshold above which
flight tracking points are considered. Set the delimiter of the input files with
kwarg `odelim`, if delimiters are any character different from whitespace. Any
`remarks` can be attached to `FlightDB`.

Alternatively, instatiate directly with the fields of `FlightDB`, where the correct
database type is checked, and wrong datasets are removed in every field.
"""
struct FlightDB{T} <: FlightSet{T}
  inventory::Vector{FlightData{T}}
  archive::Vector{FlightData{T}}
  onlineData::Vector{FlightData{T}}
  metadata::DBMetadata{T}

  """
  Unmodified constructor for `FlightDB` with basic checks for correct dataset type
  in each dataset field.
  """
  function FlightDB{T}(inventory::Vector{FlightData{T}}, archive::Vector{FlightData{T}},
    onlineData::Vector{FlightData{T}}, metadata::DBMetadata{T}) where {T}

    # Check for correct dataset type in each vector and for correct floating point precision
    inventory = checkDBtype(inventory, "VOLPE")
    archive = checkDBtype(archive, "FlightAware")
    onlineData = checkDBtype(onlineData, "flightaware.com")

    new{T}(inventory, archive, onlineData, metadata)
  end #constructor 1 FlightDB
end #struct FlightDB


"""
Modified constructor creating the database from an identifer of the
database type and the respective folder path for that database.
"""
function FlightSet{T}(DBtype::String, folder::String...; altmin::Real=5000,
  remarks=nothing, odelim::Union{Nothing,Char,String}=nothing) where T

  # Save time of database creation
  tstart = Dates.now()
  # Check DBtype addresses all folder paths
  if length(DBtype) ≠ length(folder)
    throw(ArgumentError("Number of characters in `DBtype` must match length of vararg `folder`"))
  end
  # Find database types
  i1 = [findall(isequal('i'), lowercase(DBtype)); findall(isequal('1'), DBtype)]
  i2 = [findall(isequal('a'), lowercase(DBtype)); findall(isequal('2'), DBtype)]
  i3 = [findall(isequal('o'), lowercase(DBtype)); findall(isequal('3'), DBtype)]

  # Load databases for each type
  # VOLPE AEDT inventory
  ifiles = String[]
  for i in i1
    findfiles!(ifiles, folder[i], ".csv")
  end
  inventory = loadInventory(ifiles...; Float=T, altmin=altmin)
  # FlightAware commercial archive
  ifiles = String[]
  for i in i2
    findfiles!(ifiles, folder[i], ".csv")
  end
  archive = loadArchive(ifiles...; Float=T, altmin=altmin)
  ifiles = String[]
  for i in i3
    findfiles!(ifiles, folder[i], ".tsv", ".txt", ".dat")
  end
  onlineData = loadOnlineData(ifiles...; Float=T, altmin=altmin, delim=odelim)
  tmin, tmax = if isempty([inventory; archive; onlineData])
    tstart, tstart
  else
    minimum([[f.metadata.date.start for f in inventory];
      [f.metadata.date.start for f in archive];
      [f.metadata.date.start for f in onlineData]]),
    maximum([[f.metadata.date.stop for f in inventory];
      [f.metadata.date.stop for f in archive];
      [f.metadata.date.stop for f in onlineData]])
  end

  tend = Dates.now()
  tc = tz.ZonedDateTime(tend, tz.localzone())
  loadtime = Dates.canonicalize(Dates.CompoundPeriod(tend - tstart))

  @info string("FlightDB data loaded in ",
    "$(join(loadtime.periods[1:min(2,length(loadtime.periods))], ", ")) to",
    "\n▪ inventory ($(length(inventory)) entries)\n▪ archive ($(length(archive)) entries)\n",
    "▪ onlineData ($(length(onlineData)) entries)\n▪ metadata")

  FlightDB{T}(inventory, archive, onlineData,
    DBMetadata{T}(altmin, (start=tmin, stop=tmax), tc, loadtime, remarks))
end # constructor 2 FlightDB


""" Default FlightSet/FlightData constructor for Float32 """
FlightSet(DBtype::String, folder::String...; altmin::Real=5000,
  remarks=nothing, odelim::Union{Nothing,Char,String}=nothing) =
  FlightSet{Float32}(DBtype, folder...; altmin=altmin, remarks=remarks, odelim=odelim)


## Define structs related to cloud data

"""
# struct CloudTrack{T<:AbstractFloat}

Store Data related to a single cloud track.

## Fields

- `time::Vector{DateTime}`
- `lat::Vector{T}`
- `lon::Vector{T}`
- `metadata::CloudMetadata`

Default floating point precision is `Float32`.
"""
struct CloudTrack{T<:AbstractFloat}
  data::DataFrame
  metadata::CloudMetadata

  """ Unmodified constructor for `CloudTrack` with basic checks for correct `data`"""
  function CloudTrack(data::DataFrame, metadata::CloudMetadata)

    # Column checks and warnings
    standardnames = ["time", "lat", "lon"]
    standardtypes = [Union{DateTime,Vector{DateTime}},
      Vector{<:AbstractFloat}, Vector{<:AbstractFloat}]
    bounds = (:lat => (-90, 90), :lon => (-180, 180))
    checkcols!(data, standardnames, standardtypes, bounds,
      "CloudTrack", metadata.ID)
    T = eltype(data.lon)
    new{T}(data,metadata)
  end #constructor 1 CloudTrack
end #struct CloudTrack


"""
# struct CloudDB

Database for cloud track data with fields:
- `tracks::Vector{CloudTrack}`
- `metadata::DBMetadata`
"""
struct CloudDB
  tracks::Vector{CloudTrack}
  metadata::DBMetadata

  """ unmodified constructor for CloudDB """
  CloudDB(tracks::Vector{CloudTrack}, metadata::DBMetadata) = new(tracks, metadata)

  """
  Modified constructor creating the database from mat files in the given folder
  or any subfolder using the floating point precision given by `Float`.
  """
  function CloudDB(folders::String...; Float::DataType=Float32, remarks=nothing)
    tstart = Dates.now()
    # Scan folders for HDF4 files
    files = String[]
    for folder in folders
      try findfiles!(files, folder, ".mat")
      catch
        @warn "read error; data skipped" folder
      end
    end

    # Load cloud tracks from mat files into TrackMatcher in Julia format
    tracks = loadCloudTracks(files...; Float=Float)
    # Calculate load time
    tend = Dates.now()
    tc = tz.ZonedDateTime(tend, tz.localzone())
    loadtime = Dates.canonicalize(Dates.CompoundPeriod(tend - tstart))
    # For now find min/max times in CloudTracks
    tmin = minimum(t.data.time[1] for t in tracks)
    tmax = maximum(t.data.time[end] for t in tracks)

    # Instantiate CloudDB
    new(tracks, DBMetadata(NaN, (start=tmin, stop=tmax), tc, loadtime, remarks))
  end #modified constructor 2
end #struct CloudDB

## Define structs related to sat data

"""
# struct SatData

Satellite data with fields
- `data::DataFrame`
- `metadata::SatMetadata`

The `DataFrame` of `data` has columns for the satellite time and position and indices
for data files with additional information:
`
- `time::Vector{DateTime}`                        (time stamp of measurement)
- `lat::Vector{<:Union{Missing,<:AbstractFloat}}`   (latitude)
- `lon::Vector{<:Union{Missing,<:AbstractFloat}}`   (longitude)
- `fileindex::Vector{Int}`                        (index for filenames in `SatMetadata`)


# Instantiation

    SatData(folders::String...; Float::DataType=Float32, type::Symbol=:undef,
      remarks=nothing) -> struct SatData

Construct `SatData` from any number of absolute or relative folder paths given as string.
SatData searches for hdf files in all folders recursively and determines the data type
(`CLay` or `CPro`) from the file names of the first 50 files unless the type is specified
with a Symbol `:CLay` or `:CPro` by the keyword argument `type`. Only one type of satellite
data can be stored in `SatData`. By default, all values are read in as `Float32`,
but can be set to any other precision by the `Float` kwarg.
"""
struct SatData
  data::DataFrame
  metadata::SatMetadata

  """ Unmodified constructor for `SatData` with basic checks for correct `data`"""
  function SatData(data::DataFrame, metadata::SatMetadata)
    standardnames = ["time", "lat", "lon", "fileindex"]
    standardtypes = [Vector{DateTime}, Vector{<:AbstractFloat},
      Vector{<:AbstractFloat}, Vector{Int}]
    bounds = (:time => (DateTime(2006), Dates.now()), :lat => (-90,90),
      :lon => (-180,180), :fileindex => (1, Inf))
    checkcols!(data, standardnames, standardtypes, bounds, "CLay")
    new(data, metadata)
  end #constructor 1 SatData

  """
  Modified constructor creating the database from mat files in the given `folders`
  or any subfolder using the floating point precision given by `Float`. The sat data
  `type` is determined from the first 50 files in the database unless directly
  specified `type`. Any `remarks` can be added to the metadata.
  """
  function SatData(
    folders::String...;
    Float::DataType=Float32,
    type::Symbol=:undef,
    remarks=nothing
  )
    tstart = Dates.now()
    # Scan folders for HDF4 files
    files = String[]
    for folder in folders
      try findfiles!(files, folder, ".hdf")
      catch
        @warn "read error; data skipped" folder
      end
    end
    # If type of satellite data is not defined, find it based on first 50 file names (~2 days)
    type = type == :CLay || type == :CPro ? string(type) :
      count(occursin.("CLay", files[1:min(length(files), 50)])) ≥
      count(occursin.("CPro", files[1:min(length(files), 50)])) ? "CLay" : "CPro"
    # Select files of type based on file name
    satfiles = files[occursin.(type, basename.(files))]
    wrongtype = length(files) - length(satfiles)
    wrongtype > 0 &&
      @warn "$wrongtype files with wrong satellite type detected; data skipped"
    # Create empty struct, if no data files were found
    isempty(satfiles) && return SatData(Float)
    # Start MATLAB session
    ms = mat.MSession()
    # Initialise arrays
    utc = Vector{Vector{DateTime}}(undef, length(satfiles))
    lat = Vector{Vector{Float}}(undef, length(satfiles))
    lon = Vector{Vector{Float}}(undef, length(satfiles))
    fileindex = Vector{Vector{Int}}(undef, length(satfiles))
    # Loop over files
    prog = pm.Progress(length(satfiles), "load sat data...")
    for (i, file) in enumerate(satfiles)
      # Find files with cloud layer data
      utc[i], lat[i], lon[i], fileindex[i] = try
        # Extract time
        mat.put_variable(ms, :file, file)
        mat.eval_string(ms, "clear t\ntry\nt = hdfread(file, 'Profile_UTC_Time');\nend")
        t = mat.jarray(mat.get_mvariable(ms, :t))[:,2]
        # Extract lat/lon
        mat.eval_string(ms, "clear longitude\ntry\nlongitude = hdfread(file, 'Longitude');\nend")
        longitude = mat.jarray(mat.get_mvariable(ms, :longitude))[:,2]
        mat.eval_string(ms, "clear latitude\ntry\nlatitude = hdfread(file, 'Latitude');\nend")
        latitude = mat.jarray(mat.get_mvariable(ms, :latitude))[:,2]
        # Save time converted to UTC and lat/lon
        convertUTC.(t), latitude, longitude, [i for index in t]
      catch
        # Skip data on failure and warn
        @warn "read error in CALIPSO granule; data skipped"  granule = splitext(basename(file))[1]
        DateTime[], Float[], Float[], Int[]
      end
      # Monitor progress for progress bar
      pm.next!(prog, showvalues = [(:date,Dates.Date(splitdir(dirname(file))[2], "y_m_d"))])
    end #loop over files
    pm.finish!(prog)

    # Close MATLAB session
    mat.close(ms)
    # Calculate time span of satellite data
    sattime = [utc...;]
    # Find data range
    tmin = minimum(sattime)
    tmax = maximum(sattime)
    # Save computing times
    tend = Dates.now()
    loadtime = Dates.canonicalize(Dates.CompoundPeriod(tend - tstart))

    # Instantiate new struct
    @info string("SatData data loaded in ",
      "$(join(loadtime.periods[1:min(2,length(loadtime.periods))], ", ")) to",
      "\n▪ data ($(length(sattime)) data rows)\n  – time\n  – lat\n  – lon",
      "\n  – fileindex\n▪ metadata")
    satdata = DataFrame(time=sattime, lat=[lat...;],
      lon=[lon...;], fileindex=[fileindex...;])
    new(satdata, SatMetadata(satfiles, (start=tmin, stop=tmax), loadtime, remarks=remarks))
  end #constructor 2 SatData
end #struct SatData


""" External constructor for emtpy SatData struct """
SatData(Float::DataType=Float32) = SatData(DataFrame(time=DateTime[], lat=Float[],
  lon=Float[], fileindex=Int[]), SatMetadata(Dict{Int,String}(), :undef,
  (start=Dates.now(), stop=Dates.now()), Dates.now(),
  Dates.canonicalize(Dates.CompoundPeriod())))


## Define structs related to intersection data

"""
# struct Intersection

Intersection-related data with fields
- `data::DataFrame`
- `tracked::DataFrame`
- `accuracy::DataFrame`
- `metadata::XMetadata`


## data

Data related to spatial and temporal coordinates of intersections between satellite
and flight tracks and the meteorological conditions at the intersections.

`DataFrame` columns are:
- `id::Vector{String}`: unique intersection identifier
- `lat::Vector{<:AbstractFloat}`: latitude of intersection
- `lon::Vector{<:AbstractFloat}`: longitude of intersection
- `tdiff::Vector{Dates.CompoundPeriod}`: time difference between flight and satellite overpass
- `tflight::Vector{DateTime}`: time of aircraft at intersection
- `tsat::Vector{DateTime}`: time of satellite at intersection
- `feature::Vector{<:Union{Missing,Symbol}}`: atmospheric conditions at intersection


## tracked

Original track data in the vicinity of the intersection.

`DataFrame` columns are:
- `id: Vector{String}`: unique intersection identifier
- `flight::Vector{FlightTrack}`: `FlightTrack` in the vicinity of the intersection
- `CPro::Vector{CPro}`: `CPro` CALIOP profile data in the vicinity of the intersection
- `CLay::Vector{CLay}`: `CLay` CALIOP layer data in the vicinity of the intersection


## accuracy

Measures about the accuracy of the intersection calculations and the quality of the track data.

`DataFrame` columns are:
- `id: Vector{String}`: unique intersection identifier
- `intersection::Vector{<:AbstractFloat}`: accuracy of the intersection calculation in meters
- `flightcoord::Vector{<:AbstractFloat}`: distance of nearest tracked flight data
  to calculated intersection in meters
- `satcoord::Vector{<:AbstractFloat}`: distance of nearest tracked sat data
  to calculated intersection in meters
- `flighttime::Vector{Dates.CompoundPeriod}`: time difference between measurement
  of nearest tracked flight data and calculated time of aircraft at intersection
- `sattime::Vector{Dates.CompoundPeriod}`: time difference between measurement
  of nearest tracked sat data and calculated time of satellite at intersection


# Instantiation

    Intersection(
      flights::FlightDB,
      sat::SatData,
      savesecondsattype::Bool=false;
      maxtimediff::Int=30,
      primspan::Int=0,
      secspan::Int=15,
      lidarrange::Tuple{Real,Real}=(15_000,-Inf),
      stepwidth::Real=1000,
      Xradius::Real=20_000,
      expdist::Real=Inf,
      Float::DataType=Float32,
      remarks=nothing
    ) -> struct Intersection

Construct `Intersection` from the preloaded `FlightDB` and `SatData` with the option
to save the other satellite data type not used in `sat` (either `CLay` or `CPro`),
when `savesecondsattype` is set to `true`. Folder structure and file names must
be identical only with `CLay`/`CPro` interchanged for this option to work.

The following parameters can be set to influence intersection calculations or
data saved to the struct:
- `maxtimediff::Int=30`: maximum time difference allowed between aircraft passage
  and satellite overpass at intersection
- `primspan::Int=0`: Number of additional data points of original track data
  saved in the vicinity of the intersection and stored in `Intersection.tracked.flight`
- `secspan::Int=15`: Number of additional data points of original track data
  saved in the vicinity of the intersection and stored in `Intersection.tracked.CPro`
  and `Intersection.tracked.CLay`
- `lidarrange::Tuple{Real,Real}=(15,-Inf)`: lidar measurements saved for column heights
  between `(max, min)` (set to `Inf`/`-Inf` to store all values up to top/bottom)
- `stepwidth::Real=1000`: step width of interpolation in flight and sat tracks
  in meters (partially internally converted to degrees at equator)
- `expdist::Real=Inf`: threshold for maximum distant to nearest measured track point;
  if intersection is above threshold, it will be ignored
- `Xradius::Real=20_000`: radius in meters within which multiple finds of an
  intersection are disregarded and only the most accurate is counted
- `Float::DataType=Float32`: Set the precision of floating point numbers
  (single precision by default)
- `remarks=nothing`: any data or remarks attached to the metadata

Or construct `Intersection` by directly handing over the different `DataFrame`s where the names, order,
and types of each columns are checked and attempted to correct, together with the metadata:

    function Intersection(
      data::DataFrame,
      tracked::DataFrame,
      accuracy::DataFrame,
      metadata::XMetadata
    ) -> struct Intersection
"""
struct Intersection
  data::DataFrame
  tracked::DataFrame
  accuracy::DataFrame
  metadata::XMetadata


  """ Unmodified constructor for `Intersection` """
  function Intersection(
    data::DataFrame,
    tracked::DataFrame,
    accuracy::DataFrame,
    metadata::XMetadata
  )
    # Check data
    standardnames = ["id", "lat", "lon", "alt", "tdiff", "tflight", "tsat", "feature"]
    standardtypes = [Vector{String}, Vector{<:AbstractFloat}, Vector{<:AbstractFloat},
      Vector{<:AbstractFloat}, Vector{Dates.CompoundPeriod}, Vector{DateTime},
      Vector{DateTime}, Vector{<:Union{Missing,Symbol}}]
    bounds = (:lat => (-90,90), :lon => (-180,180), :alt => (0, Inf))
    checkcols!(data, standardnames, standardtypes, bounds, "Intersection.data")
    # Check tracked (measured data)
    standardnames = ["id", "flight", "CPro", "CLay"]
    standardtypes = [Vector{String}, Vector{FlightTrack}, Vector{CPro}, Vector{CLay}]
    bounds = ()
    checkcols!(tracked, standardnames, standardtypes, bounds, "Intersection.tracked",
      essentialcols = [1])
    # Check accuracy
    standardnames = ["id", "intersection", "flightcoord", "satcoord", "flighttime", "sattime"]
    standardtypes = [Vector{String}, Vector{<:AbstractFloat}, Vector{<:AbstractFloat},
      Vector{<:AbstractFloat}, Vector{Dates.CompoundPeriod}, Vector{Dates.CompoundPeriod}]
    bounds = ()
    checkcols!(accuracy, standardnames, standardtypes, bounds, "Intersection.accuracy",
      essentialcols = [1])
    new(data, tracked, accuracy, metadata)
  end #constructor 1 Intersection


  """ Modified constructor with some automated calculations of the flight intersection data. """
  function Intersection(
    flights::FlightDB,
    sat::SatData,
    savesecondsattype::Bool=false;
    maxtimediff::Int=30,
    primspan::Int=0,
    secspan::Int=15,
    lidarrange::Tuple{Real,Real}=(15_000,-Inf),
    stepwidth::Real=0.01,
    Xradius::Real=20_000,
    expdist::Real=Inf,
    Float::DataType=Float32,
    remarks=nothing
  )
    # Combine all datasets and find intersections
    track = [[getfield(flights, f) for f in fieldnames(FlightDB)[1:end-1]]...;]
    intersection(track, flights.metadata, sat, savesecondsattype, maxtimediff,
      primspan, secspan, lidarrange, stepwidth, Xradius, expdist, Float, remarks)
  end #constructor 2 Intersection


  """ Modified constructor with some automated calculations of the cloud intersection data. """
  function Intersection(
    cloud::CloudDB,
    sat::SatData,
    savesecondsattype::Bool=false;
    maxtimediff::Int=30,
    primspan::Int=0,
    secspan::Int=15,
    lidarrange::Tuple{Real,Real}=(15_000,-Inf),
    stepwidth::Real=0.01,
    Xradius::Real=20_000,
    expdist::Real=Inf,
    Float::DataType=Float32,
    remarks=nothing
  )
    # Combine all datasets and find intersections
    intersection(cloud.tracks, cloud.metadata, sat, savesecondsattype, maxtimediff,
      primspan, secspan, lidarrange, stepwidth, Xradius, expdist, Float, remarks)
  end #constructor 3 Intersection
end #struct Intersection
