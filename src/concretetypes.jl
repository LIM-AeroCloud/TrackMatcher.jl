## Define own Metadata structs
"""
# struct FlightMetadata{T}

Immutable struct to hold metadata for `FlightTrack` of the `FlightDB` with fields

- `dbID::Union{Int,AbstractString}`
- `flightID::Union{Missing,AbstractString}`
- `route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}}`
- `aircraft::Union{Missing,AbstractString}`
- `date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}`
- `area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,AbstractFloat}}`
- `flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,AbstractFloat,AbstractFloat}}}}`
- `useLON::Bool`
- `source::AbstractString`
- `file::AbstractString`

By default, all `AbstractFloat` are set to Float32 but can be set to any other precision
in `FlightDB` with the kwarg `Float`.

## dbID
Database ID – integer counter for `inventory`,
String with information about `FlightID`, `route`, and/or scheduled arrival for
FlightAware data.

## FlightID and aircraft
Strings with aircraft identification and type.

## route
`NamedTuple` with fields for `orig`in and `dest`ination holding the airport codes.

## area
`NamedTuple` with fields for latitude and Longitude range. For the longitude range,
it is distinguished between positive and negative ranges to avoid problems with
flights passing the date line.

Fields:
- `latmin`
- `latmax`
- `elonmin`
- `elonmax`
- `wlonmin`
- `wlonmax`

## date
`NamedTuple` with fields `start` and `stop` for start and end time of the current
flight.

## flex
`Tuple` of `NamedTuple`s with entries
- `range` (`UnitRange`): flight segment between inflection points of the current flight track
- `min` (`AbstractFloat`): minimum x value in the flight segment
- `max` (`AbstractFloat`): maximum x value in the flight segment

## useLON
Flag (`Bool`) whether to use longitude as x data for track interpolation.

## source
String describing the database source of the current flight:
- `"VOLPE AEDT"`
- `"FlightAware"`
- `"flightaware.com"`

## file
String holding the absolute folder path and file name.


# Instantiation

`FlightMetadata` is constructed automatically, when `FlightTrack` is instatiated using
a modified constructor and `dbID`, `flightID`, `aircraft` type, `route`, `useLON`,
`flex`, `source`, and `file`.
Fields `area` and `date` are calculated from `lat`/`lon`, and `date` vectors.

    FlightMetadata(
      dbID::Union{Int,AbstractString},
      flightID::Union{Missing,AbstractString},
      route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
      aircraft::Union{Missing,AbstractString},
      date::Vector{DateTime},
      lat::Vector{<:Union{Missing,<:AbstractFloat}},
      lon::Vector{<:Union{Missing,<:AbstractFloat}},
      useLON::Bool,
      flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,AbstractFloat,AbstractFloat}}}},
      source::AbstractString,
      file::AbstractString
    ) -> struct FlightMetadata

Or construct `FlightMetadata` by directly handing over every field:

    FlightMetadata(
      dbID::Union{Int,AbstractString},
      flightID::Union{Missing,AbstractString},
      route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
      aircraft::Union{Missing,AbstractString},
      date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
      area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,AbstractFloat}},
      flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange,AbstractFloat,AbstractFloat}}}},
      useLON::Bool,
      source::AbstractString,
      file::AbstractString
    ) -> struct FlightMetadata
"""
struct FlightMetadata{T} <: FlightTrack{T}
  dbID::Union{Int,AbstractString}
  flightID::Union{Missing,AbstractString}
  route::Union{Missing,NamedTuple{(:orig,:dest),Tuple{AbstractString,AbstractString}}}
  aircraft::Union{Missing,AbstractString}
  date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,T}}
  flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}}
  useLON::Bool
  source::AbstractString
  file::AbstractString

  """ Unmodified constructor for `FlightMetadata` """
  function FlightMetadata{T}(
    dbID::Union{Int,AbstractString},
    flightID::Union{Missing,AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
    aircraft::Union{Missing,AbstractString},
    date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,T}},
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}},
    useLON::Bool,
    source::AbstractString,
    file::AbstractString
  ) where T
    # T = typeof(area.latmin)
    new{T}(dbID, flightID, route, aircraft, date, area, flex, useLON, source, file)
  end #constructor 1 FlightMetadata

  """
  Modified constructor for FlightMetadata with some automated construction of fields
  and variable checks.
  """
  function FlightMetadata{T}(
    dbID::Union{Int,AbstractString},
    flightID::Union{Missing,AbstractString},
    route::Union{Missing,NamedTuple{(:orig,:dest),<:Tuple{AbstractString,AbstractString}}},
    aircraft::Union{Missing,AbstractString},
    date::Vector{DateTime},
    lat::Vector{<:Union{Missing,<:AbstractFloat}},
    lon::Vector{<:Union{Missing,<:AbstractFloat}},
    useLON::Bool,
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}},
    source::AbstractString,
    file::AbstractString
  ) where T
    # T = promote_type(eltype(lat), eltype(lon))
    elonmax = isempty(lon[lon.≥0]) ? T(NaN) : maximum(lon[lon.≥0])
    elonmin = isempty(lon[lon.≥0]) ? T(NaN) : minimum(lon[lon.≥0])
    wlonmax = isempty(lon[lon.<0]) ? T(NaN) : maximum(lon[lon.<0])
    wlonmin = isempty(lon[lon.<0]) ? T(NaN) : minimum(lon[lon.<0])
    area = (latmin=minimum(lat), latmax=maximum(lat),
      elonmin=elonmin, elonmax=elonmax, wlonmin=wlonmin, wlonmax=wlonmax)
    new{T}(dbID, flightID, route, aircraft, (start=date[1], stop=date[end]), area,
      flex, useLON, source, file)
  end #constructor 2 FlightMetadata
end #struct FlightMetadata

""" External constructor for emtpy FlightMetadata struct """
FlightMetadata{T}() where T = FlightMetadata("", missing, missing, missing,
  (start=Dates.now(), stop=Dates.now()), (latmin=T(NaN), latmax=T(NaN),
  elonmin=T(NaN), elonmax=T(NaN), wlonmin=T(NaN), wlonmax=T(NaN)),
  ((range=0:0, min=T(NaN), max=T(NaN)),), true, "","")


"""
# struct CloudMetadata

Currently only place holder for remarks available during test phase.
"""
struct CloudMetadata{T}
  ID::String
  date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,T}}
  flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}}
  useLON::Bool
  file::String

  """ Unmodified constructor for `CloudMetadata` """
  function CloudMetadata(
    ID::String,
    date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,T}},
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}},
    useLON::Bool,
    file::String
  ) where T<:AbstractFloat
    # T = typeof(area.latmin)
    new{T}(ID, date, area, flex, useLON, file)
  end #constructor 1 CloudMetadata

  """
  Modified constructor for CloudMetadata with some automated construction of fields
  and variable checks.
  """
  function CloudMetadata(
    ID::Union{Int,AbstractString},
    data::DataFrame,
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}},
    useLON::Bool,
    file::AbstractString
  ) where T<:AbstractFloat
    # T = promote_type(eltype(data.lat), eltype(data.lon))
    elonmax = isempty(data.lon[data.lon.≥0]) ? T(NaN) : maximum(data.lon[data.lon.≥0])
    elonmin = isempty(data.lon[data.lon.≥0]) ? T(NaN) : minimum(data.lon[data.lon.≥0])
    wlonmax = isempty(data.lon[data.lon.<0]) ? T(NaN) : maximum(data.lon[data.lon.<0])
    wlonmin = isempty(data.lon[data.lon.<0]) ? T(NaN) : minimum(data.lon[data.lon.<0])
    area = (latmin=minimum(data.lat), latmax=maximum(data.lat),
      elonmin=elonmin, elonmax=elonmax, wlonmin=wlonmin, wlonmax=wlonmax)
    new{T}(ID, (start=data.time[1], stop=data.time[end]), area, flex, useLON, file)
  end #constructor 2 CloudMetadata
end #struct CloudMetadata


"""
# struct SatMetadata

Immutable struct to hold metadata for `SatData` with fields

- `files::Dict{Int,String}`
- `type::Symbol`
- `date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}`
- `created::Union{DateTime,ZonedDateTime}`
- `loadtime::Dates.CompoundPeriod`
- `remarks`

## files
Dictionary with indices of `fileindex` column in `SatData.data` pointing to the
full file names.

## type
Symbol indicating, whether profile or layer data is stored.

## date
`NamedTuple` with fields `start` and `stop` for start and end time of the monitored
satellite period.

## created
time of creation of database

## loadtime
time it took to read data files and load it to the struct

##remarks
any additional data or comments that can be attached to the database


# Instantiation

`SatMetadata` is constructed automatically, when `SatData` is instatiated using
a modified constructor and `files`, `date`, `loadtime`, and `remarks`.

    function SatMetadata(
      files::Vector{String},
      date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
      loadtime::Dates.CompoundPeriod=Dates.canonicalize(Dates.CompoundPeriod());
      remarks=nothing
    ) -> struct SattMetadata
"""
struct SatMetadata
  files::Dict{Int,String}
  type::Symbol
  date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  created::Union{DateTime,ZonedDateTime}
  loadtime::Dates.CompoundPeriod
  remarks

  function SatMetadata(
    files::Dict{Int,String},
    type::Symbol,
    date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    created::Union{DateTime,ZonedDateTime},
    loadtime::Dates.CompoundPeriod,
    remarks=nothing
  )
    new(files, type, date, created, loadtime, remarks)
  end #constructor 1 SatMetadata

  function SatMetadata(
    files::Vector{String},
    date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    loadtime::Dates.CompoundPeriod=Dates.canonicalize(Dates.CompoundPeriod());
    remarks=nothing
  )
    # Find type of satellite data based on first 50 files (~2 days)
    type = occursin("CLay", files[1]) ≥
      count(occursin.("CPro", files[1:min(length(files), 50)])) ? :CLay : :CPro
    # Create a new instance of SatMetadata
    new(Dict(enumerate(files)), type, date, tz.now(tz.localzone()), loadtime, remarks)
  end #constructor 2 SatMetadata
end #struct SatMetadata


"""
# struct DBMetadata

Immutable struct with additional information of databases:

- `altmin::Real`: Minimum altitude threshold for which flight data is considered
- `date`: NamedTuple with entries `start`/`stop` giving the time range of the database
- `created`: time of creation of database
- `loadtime`: time it took to read data files and load it to the struct
- `remarks`: any additional data or comments that can be attached to the database
"""
struct DBMetadata{T} <: PrimarySet{T}
  altmin::T
  date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  created::Union{DateTime,ZonedDateTime}
  loadtime::Dates.CompoundPeriod
  remarks
end #struct DBMetadata


"""
# struct XMetadata

Immutable struct with additional information of intersection data:

- `maxtimediff`: maximum time difference in minutes allowed between
  satellite overpass and aircraft passing at intersection
- `stepwidth`: stepwidth (in meters, partially internally converted to degrees at equator)
  used to interpolate track data
- `Xradius`: radius in meters around an intersection in which further intersections
  will be removed as duplicates due to the interpolation algorithm
- `lidarrange::NamedTuple{(:top,:bottom),Tuple{Real,Real}}`: user defined level thresholds
  for top/bottom heights for which CALIPSO data is considered
- `lidarprofile::NamedTuple`: CALIPSO lidar height levels used in the current dataset
- `sattype::Symbol`: cloud layer or profile data types used for the satellite date
- `satdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}`: date range of
  satellite data
- `altmin::Real`: minimum threshold for which flight data is considered
- `flightdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}`: date range
  of flight data
- `created`: time of creation of database
- `loadtime`: time it took to find intersections and load it to the struct
- `remarks`: any additional data or comments that can be attached to the database

XMetadata can be instantiated using a `Tuple` or `NamedTuple` for the `lidarrange`
"""
struct XMetadata
  maxtimediff::Int
  stepwidth::Real
  Xradius::Real
  expdist::Real
  lidarrange::NamedTuple{(:top,:bottom),Tuple{Real,Real}}
  lidarprofile::NamedTuple
  sattype::Symbol
  satdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  altmin::Real
  flightdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  created::Union{DateTime,ZonedDateTime}
  loadtime::Dates.CompoundPeriod
  remarks

  function XMetadata(
    maxtimediff::Int,
    stepwidth::Real,
    Xradius::Real,
    expdist::Real,
    lidarrange::NamedTuple{(:top,:bottom),Tuple{Real,Real}},
    lidarprofile::NamedTuple,
    sattype::Symbol,
    satdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    altmin::Real,
    flightdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    created::Union{DateTime,ZonedDateTime},
    loadtime::Dates.CompoundPeriod,
    remarks
  )
    new(maxtimediff, stepwidth, Xradius, expdist, lidarrange, lidarprofile,
      sattype, satdates, altmin, flightdates, created, loadtime, remarks)
  end #constructor 1 XMetaData

  function XMetadata(
    maxtimediff::Int,
    stepwidth::Real,
    Xradius::Real,
    expdist::Real,
    lidarrange::Tuple{Real,Real},
    lidarprofile::NamedTuple,
    sattype::Symbol,
    satdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    altmin::Real,
    flightdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    created::Union{DateTime,ZonedDateTime},
    loadtime::Dates.CompoundPeriod,
    remarks=nothing
  )
    new(maxtimediff, stepwidth, Xradius, expdist, (top=lidarrange[1], bottom=lidarrange[2]),
      lidarprofile, sattype, satdates, altmin, flightdates, created, loadtime, remarks)
  end #constructor 2 XMetaData
end #struct XMetaData


## Define structs for satellite data output

"""
# struct CLay

CALIOP cloud layer `data` stored in a `DataFrame` with columns:
- `time::Vector{DateTime}` (time index)
- `lat::Vector{AbstractFloat}` (latitude position of current time index)
- `lon::Vector{AbstractFloat}` (longitude position of current time index)
- `layer::Vector{NamedTuple{(:top,:base),Tuple{Vector{<:AbstractFloat},Vector{<:AbstractFloat}}}}`
  (layer top/base heights in meters)
- `feature::Vector{Vector{Symbol}}` (symbol features of a layer)
- `OD::Vector{<:Vector{<:AbstractFloat}}` (layer optical depth)
- `IWP::Vector{<:Vector{<:Union{Missing,<:AbstractFloat}}}` (layer ice water path)
- `Ttop::Vector{<:Vector{<:AbstractFloat}}` (layer top temperature)
- `Htropo::Vector{<:AbstractFloat}` (tropopause height at current time index)
- `night::BitVector` (flag for nights (`true`))
- `averaging::Vector{<:Vector{Int}}` (horizontal averaging in km)

# Instantiation

    CLay(ms::mat.MSession, files::Vector{String},
      lidarrange::Tuple{Real,Real}=(15_000,-Inf), altmin::Real=5000, Float::DataType=Float32)

Construct `CLay` from a list of file names (including directories) and a running
MATLAB session `ms` and save data, if layers are within the bounds
of `lidarrange` and above flight `altmin` threshold. By default, all values are
saved as `Float32`, but can be set to any other precision by the `Float` kwarg.

Or construct `CLay` by directly handing over the `DataFrame` where the names, order,
and types of each columns are checked and attempted to correct:

    CLay(data::DataFrame) -> struct CLay
"""
struct CLay
  data::DataFrame

  """ Unmodified constructor for `CLay` """
  function CLay(data::DataFrame)
    standardnames = ["time", "lat", "lon",
      "layer", "feature", "OD", "IWP", "Ttop", "Htropo", "night", "averaging"]
      standardtypes = [Vector{DateTime}, Vector{<:AbstractFloat}, Vector{<:AbstractFloat},
      Vector{NamedTuple{(:top,:base),Tuple{T,T}}} where T<:Vector{<:AbstractFloat},
      Vector{Vector{Symbol}}, Vector{<:Vector{<:AbstractFloat}},
      Vector{<:Vector{<:Union{Missing,<:AbstractFloat}}}, Vector{<:Vector{<:AbstractFloat}},
      Vector{<:AbstractFloat}, BitVector, Vector{<:Vector{Int}}]
    bounds = (:lat => (-90,90), :lon => (-180,180))
    checkcols!(data, standardnames, standardtypes, bounds, "CLay")
    new(data)
  end #constructor 1 CLay

  """
  Modified constructor of `CLay` reading data from hdf `files` using MATLAB session `ms`
  in the `lidarrange` (top to bottom), if data is above `altmin`.
  """
  function CLay(ms::mat.MSession, files::Vector{String},
    lidarrange::Tuple{Real,Real}=(15_000,-Inf), altmin::Real=5000, Float::DataType=Float32)
    # Return default empty struct if files are empty
    isempty(files) && return CLay(Float)
    # Initialise arrays
    # essential data
    utc = Vector{Vector{DateTime}}(undef, length(files))
    lat = Vector{Vector{Float}}(undef, length(files))
    lon = Vector{Vector{Float}}(undef, length(files))
    # non-essential data
    layers = Vector{Vector{NamedTuple{(:top,:base),Tuple{Vector{Float},Vector{Float}}}}}(undef, length(files))
    feature = Vector{Vector{Vector{Symbol}}}(undef,length(files))
    OD = Vector{Vector{Vector{Float}}}(undef,length(files))
    IWP = Vector{Vector{Vector{Union{Missing,Float}}}}(undef,length(files))
    Ttop = Vector{Vector{Vector{Float}}}(undef,length(files))
    Htropo = Vector{Vector{Float}}(undef, length(files))
    night = Vector{BitVector}(undef, length(files))
    averaging = Vector{Vector{Vector{Int}}}(undef,length(files))
    # Convert mininmum flight altitude to meters
    altmin = ft2km(Float(altmin))
    # Loop over files
    for (i, file) in enumerate(files)
      ## Retrieve cloud layer data; assumes faulty files are filtered by SatData
      # Extract time
      mat.put_variable(ms, :file, file)
      mat.eval_string(ms, "clear t\ntry\nt = hdfread(file, 'Profile_UTC_Time');\nend")
      utc[i] = convertUTC.(mat.jarray(mat.get_mvariable(ms, :t))[:,2])
      # Extract lat/lon
      mat.eval_string(ms, "clear longitude\ntry\nlongitude = hdfread(file, 'Longitude');\nend")
      lon[i] = mat.jarray(mat.get_mvariable(ms, :longitude))[:,2]
      mat.eval_string(ms, "clear latitude\ntry\nlatitude = hdfread(file, 'Latitude');\nend")
      lat[i] = mat.jarray(mat.get_mvariable(ms, :latitude))[:,2]
      # Save time converted to UTC and lat/lon
      # utc[i], lon[i], lat[i] = convertUTC.(t), longitude, latitude

      ## Extract layer top/base, layer features and optical depth from hdf files
      mat.eval_string(ms, "clear basealt\ntry\nbasealt = hdfread(file, 'Layer_Base_Altitude');\nend")
      mat.eval_string(ms, "clear topalt\ntry\ntopalt = hdfread(file, 'Layer_Top_Altitude');\nend")
      basealt = mat.jarray(mat.get_mvariable(ms, :basealt))
      topalt = mat.jarray(mat.get_mvariable(ms, :topalt))
      mat.eval_string(ms, "clear FCF\ntry\nFCF = hdfread(file, 'Feature_Classification_Flags');\nend")
      FCF = mat.jarray(mat.get_mvariable(ms, :FCF))
      mat.eval_string(ms, "clear FOD\ntry\nFOD = hdfread(file, 'Feature_Optical_Depth_532');\nend")
      FOD = mat.jarray(mat.get_mvariable(ms, :FOD))
      mat.eval_string(ms, "clear IWPath\ntry\nIWPath = hdfread(file, 'Ice_Water_Path');\nend")
      IWPath = mat.jarray(mat.get_mvariable(ms, :IWPath))
      mat.eval_string(ms, "clear LTT\ntry\nLTT = hdfread(file, 'Layer_Top_Temperature');\nend")
      LTT = mat.jarray(mat.get_mvariable(ms, :LTT))
      mat.eval_string(ms, "clear Htropo\ntry\nHtropo = hdfread(file, 'Tropopause_Height');\nend")
      Htropo[i] = 1000vec(mat.jarray(mat.get_mvariable(ms, :Htropo)))
      mat.eval_string(ms, "clear daynight\ntry\ndaynight = hdfread(file, 'Day_Night_Flag');\nend")
      night[i] = Bool.(vec(mat.jarray(mat.get_mvariable(ms, :daynight))))
      mat.eval_string(ms, "clear average\ntry\naverage = hdfread(file, 'Horizontal_Averaging');\nend")
      horav = mat.jarray(mat.get_mvariable(ms, :average))
      # Loop over data and convert to TrackMatcher format
      layer = Vector{NamedTuple{(:top,:base),Tuple{Vector{Float},Vector{Float}}}}(undef,length(utc[i]))
      feat = Vector{Vector{Symbol}}(undef,length(utc[i]))
      optdepth = Vector{Vector{Float}}(undef,length(utc[i]))
      icewater = Vector{Vector{Union{Missing,Float}}}(undef,length(utc[i]))
      toptemp = Vector{Vector{Float}}(undef,length(utc[i]))
      average = Vector{Vector{Int}}(undef,length(utc[i]))
      for n = 1:length(utc[i])
        l = findall((basealt[n,:] .> 0) .& (topalt[n,:] .> 0) .& (basealt[n,:] .< lidarrange[1]) .&
          (topalt[n,:] .> lidarrange[2]) .& (topalt[n,:] .> altmin))
        layer[n], feat[n], optdepth[n], toptemp[n], icewater[n], average[n] = if isempty(l)
          (top = Float[], base = Float[]), Symbol[],
          Float[], Float[], Float[], Int[]
        else
          l = findall((basealt[n,:] .> 0) .& (topalt[n,:] .> 0) .& (basealt[n,:] .< lidarrange[1]) .&
            (topalt[n,:] .> lidarrange[2]))
          (top = [1000topalt[n, m] for m in l] , base = [1000basealt[n, m] for m in l]),
          [feature_classification(classification(FCF[n,m])...) for m in l],
          [FOD[n,m] for m in l],
          [LTT[n,m] for m in l],
          [IWPath[n,m] == -9999 ? missing : IWPath[n,m] for m in l],
          [1000horav[n,m] for m in l]
        end
      end # loop over time steps in current file
      layers[i], feature[i], OD[i], IWP[i], Ttop[i], averaging[i] =
        layer, feat, optdepth, icewater, toptemp, average
    end #loop over files

    # Construct and standardise data
    data = DataFrame(time=[utc...;], lat=[lat...;], lon=[lon...;],
      layer=[layers...;], feature=[feature...;], OD=[OD...;],
      IWP=[IWP...;], Ttop=[Ttop...;], Htropo = [Htropo...;],
      night = [night...;], averaging = [averaging...;])
    # Save time, lat/lon arrays in CLay struct
    new(data)
  end #constructor 2 CLay
end #struct CLay


""" External constructor for emtpy CLay struct """
CLay(Float::DataType=Float32) = CLay(DataFrame(time = DateTime[], lat = Float[], lon = Float[],
  layer = NamedTuple{(:top,:base),Tuple{Vector{Float},Vector{Float}}}[],
  feature = Vector{Symbol}[], OD = Vector{Float}[], IWP = Vector{Float}[],
  Ttop = Vector{Float}[], Htropo = Float[], night = BitVector(),
  averaging = Vector{Int}[]))


"""
# struct CPro

CALIOP cloud profile `data` stored in a `DataFrame` with columns:
- `time::Vector{DateTime}` (current time index)
- `lat::Vector{AbstractFloat}` (latitude coordinate for current time index)
- `lon::Vector{AbstractFloat}` (lonitude coordinate for current time index)
- `feature::Vector{<:Vector{<:Union{Missing,Symbol}}}`
  (feature symbols for every height level at current time index)
- `EC532::Vector{<:Vector{<:Union{Missing,AbstractFloat}}}`
  (extinction coefficient at 532nm at every height level in current time index)

# Instantiation

    CPro(ms::mat.MSession, files::Vector{String}, sattime::Vector{DateTime}, lidarprofile::NamedTuple,
      Float::DataType=Float32)
      -> struct CPro

Construct `CPro` from a list of file names (including directories) and a running
MATLAB session `ms`. CPro data is only stored in the vicinity of intersections for
the designated `sattime`. Column data is stored height-resolved as defined by the
`lidarprofile`. By default, all values are stored as `Float32`, but can be set to
any other precision by the `Float` kwarg.

Or construct `CPro` by directly handing over the `DataFrame` where the names, order,
and types of each columns are checked and attempted to correct:

    CPro(data::DataFrame) -> struct CPro
"""
struct CPro
  data::DataFrame

  """ unmodified constructor """
  function CPro(data::DataFrame)
    standardnames = ["time", "lat", "lon", "feature", "EC532", "Htropo", "temp",
      "pressure", "rH", "IWC", "deltap", "CADscore", "night"]
    standardtypes = [Vector{DateTime}, Vector{<:AbstractFloat}, Vector{<:AbstractFloat},
      Vector{<:Vector{<:Union{Missing,Symbol}}}, Vector{<:Vector{<:Union{Missing,<:AbstractFloat}}},
      Vector{<:AbstractFloat}, Vector{<:Vector{<:Union{Missing,<:AbstractFloat}}},
      Vector{<:Vector{<:Union{Missing,<:AbstractFloat}}}, Vector{<:Vector{<:Union{Missing,<:AbstractFloat}}},
      Vector{<:Vector{<:Union{Missing,<:AbstractFloat}}}, Vector{<:Vector{<:Union{Missing,<:AbstractFloat}}},
      Vector{<:Vector{<:Union{Missing,Int8}}}, BitVector]
    bounds = (:lat => (-90,90), :lon => (-180,180), :Htropo => (4000,22_000),
      :temp => (-120,60), :pressure => (1,1086), :rH => (0,1.5), :IWC => (0,0.54),
      :deltap => (0,1), :CADscore => (-101,106))
    checkcols!(data, standardnames, standardtypes, bounds, "CPro")
    new(data)
  end #constructor 1 CPro

  """
  Modified constructor of `CPro` reading data from hdf `files` for all given `sattime` indices
  using MATLAB session `ms` and `lidarprofile` data, if data is above `altmin`.
  """
  function CPro(ms::mat.MSession, files::Vector{String}, sattime::Vector{DateTime},
    lidarprofile::NamedTuple, Float::DataType=Float32)
    # Return default empty struct if files are empty
    isempty(files) && return CPro(Float)
    # Initialise arrays
    # essential data
    utc = Vector{Vector{DateTime}}(undef, length(files))
    lat = Vector{Vector{Float}}(undef, length(files))
    lon = Vector{Vector{Float}}(undef, length(files))
    fcf = Vector{Vector{Vector{<:Union{Missing,UInt16}}}}(undef, length(files))
    # non-essential data
    ec532 = Vector{Vector{Vector{<:Union{Missing,Float}}}}(undef, length(files))
    Htropo = Vector{Vector{Float}}(undef, length(files))
    temp = Vector{Vector{Vector{<:Union{Missing,Float}}}}(undef, length(files))
    pres = Vector{Vector{Vector{<:Union{Missing,Float}}}}(undef, length(files))
    rH = Vector{Vector{Vector{<:Union{Missing,Float}}}}(undef, length(files))
    iwc = Vector{Vector{Vector{<:Union{Missing,Float}}}}(undef, length(files))
    deltap = Vector{Vector{Vector{<:Union{Missing,Float}}}}(undef, length(files))
    cad = Vector{Vector{Vector{<:Union{Missing,Int8}}}}(undef, length(files))
    night = Vector{BitVector}(undef, length(files))
    # Loop over files with cloud profile data
    for (i, file) in enumerate(files)
      ## Retrieve cloud profile data; assumes faulty files are filtered by SatData
      # Extract time
      mat.put_variable(ms, :file, file)
      mat.eval_string(ms, "clear t\ntry\nt = hdfread(file, 'Profile_UTC_Time');\nend")
      utc[i] = convertUTC.(mat.jarray(mat.get_mvariable(ms, :t))[:,2])
      # Extract lat/lon
      mat.eval_string(ms, "clear longitude\ntry\nlongitude = hdfread(file, 'Longitude');\nend")
      lon[i] = mat.jarray(mat.get_mvariable(ms, :longitude))[:,2]
      mat.eval_string(ms, "clear latitude\ntry\nlatitude = hdfread(file, 'Latitude');\nend")
      lat[i] = mat.jarray(mat.get_mvariable(ms, :latitude))[:,2]
      fcf[i] = get_lidarcolumn(UInt16, ms, "Atmospheric_Volume_Description", lidarprofile,
        coarse=false)
      # Extract non-essential data
      ec532[i] = get_lidarcolumn(Float, ms, "Extinction_Coefficient_532", lidarprofile,
        missingvalues = -9999)
      mat.eval_string(ms, "clear Htropo\ntry\nHtropo = hdfread(file, 'Tropopause_Height');\nend")
      Htropo[i] = 1000vec(mat.jarray(mat.get_mvariable(ms, :Htropo)))
      temp[i] = get_lidarcolumn(Float, ms, "Temperature", lidarprofile, missingvalues = -9999)
      pres[i] = get_lidarcolumn(Float, ms, "Pressure", lidarprofile, missingvalues = -9999)
      rH[i] = get_lidarcolumn(Float, ms, "Relative_Humidity", lidarprofile, missingvalues = -9999)
      iwc[i] = get_lidarcolumn(Float, ms, "Ice_Water_Content_Profile", lidarprofile,
        missingvalues = -9999)
      deltap[i] = get_lidarcolumn(Float, ms, "Particulate_Depalarization_Ratio_Profile_532",
        lidarprofile, missingvalues = -9999)
      cad[i] = get_lidarcolumn(Int8, ms, "CAD_Score", lidarprofile, coarse=false,
        missingvalues = -127)
      mat.eval_string(ms, "clear daynight\ntry\ndaynight = hdfread(file, 'Day_Night_Flag');\nend")
      night[i] = Bool.(vec(mat.jarray(mat.get_mvariable(ms, :daynight))))
    end #loop over files

    # Rearrange time vector and get time range
    utc = [utc...;]
    idx = [findfirst(utc .== t) for t in sattime]
    # Rearrange FCF vector and convert to symbols
    fcf = [fcf...;]
    avd =  Vector{Vector{Union{Missing,Symbol}}}(undef, length(fcf))
    for i = 1:length(fcf)
      vect = Vector{Union{Missing,Symbol}}(undef, length(fcf[i]))
      for j = 1:length(fcf[i])
        vect[j] = ismissing(fcf[i][j]) ? missing :
          feature_classification(classification(fcf[i][j])...)
      end
      avd[i] = vect
    end
    # Construct and standardise data
    data = DataFrame(time=utc[idx], lat=[lat...;][idx], lon=[lon...;][idx],
      feature=avd[idx], EC532=[ec532...;][idx], Htropo = [Htropo...;][idx],
      temp=[temp...;][idx], pressure = [pres...;][idx], rH = [rH...;][idx],
      IWC = [iwc...;][idx], deltap = [deltap...;][idx],
      CADscore = [cad...;][idx], night = [night...;][idx])
    # Save time, lat/lon arrays, and feature classification flags (FCF) in CPro struct
    new(data)
  end #constructor 2 CPro
end #struct CPro


""" External constructor for emtpy CPro struct """
CPro(Float::DataType=Float32) = CPro(DataFrame(time = DateTime[], lat = Float[], lon = Float[],
  feature = Vector{Symbol}[], EC532 = Vector{Float}[], Htropo = Float[], temp = Vector{Float}[],
  pressure = Vector{Float}[], rH = Vector{Float}[], IWC = Vector{Float}[],
  deltap = Vector{Float}[], CADscore = Vector{Int8}[], night = BitVector()))
