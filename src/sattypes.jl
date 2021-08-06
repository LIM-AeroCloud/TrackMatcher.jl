## Define Metadata structs

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

`SatMetadata` is constructed automatically, when `SatData` is instantiated using
a modified constructor and `files`, `date`, `loadtime`, and `remarks`.

    function SatMetadata(
      files::Vector{String},
      date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
      loadtime::Dates.CompoundPeriod=Dates.canonicalize(Dates.CompoundPeriod());
      remarks=nothing
    ) -> struct SattMetadata

Or hand over the individual fields to the unmodified constructor:

    function SatMetadata{T}(
      files::Dict{Int,String},
      type::Symbol,
      date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
      created::Union{DateTime,ZonedDateTime},
      loadtime::Dates.CompoundPeriod,
      remarks=nothing
    ) where T
"""
struct SatMetadata{T} <: SatTrack{T}
  files::Dict{Int,String}
  type::Symbol
  date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  created::Union{DateTime,ZonedDateTime}
  loadtime::Dates.CompoundPeriod
  remarks

  function SatMetadata{T}(
    files::Dict{Int,String},
    type::Symbol,
    date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    created::Union{DateTime,ZonedDateTime},
    loadtime::Dates.CompoundPeriod,
    remarks=nothing
  ) where T
    new{T}(files, type, date, created, loadtime, remarks)
  end #constructor 1 SatMetadata

  function SatMetadata{T}(
    files::Vector{String},
    date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    loadtime::Dates.CompoundPeriod=Dates.canonicalize(Dates.CompoundPeriod());
    remarks=nothing
  ) where T
    # Find type of satellite data based on first 50 files (~2 days)
    type = occursin("CLay", files[1]) ≥
      count(occursin.("CPro", files[1:min(length(files), 50)])) ? :CLay : :CPro
    # Create a new instance of SatMetadata
    new(Dict(enumerate(files)), type, date, tz.now(tz.localzone()), loadtime, remarks)
  end #constructor 2 SatMetadata
end #struct SatMetadata

"""
    SatMetadata{T}() where T

External constructor for empty `SatMetadata` with floating point precision `T`.
"""
SatMetadata{T}() where T = SatMetadata{T}(
  Dict{Int,String}(), :undef, (start=Dates.now(), stop=Dates.now()),
  Dates.now(), Dates.CompoundPeriod(), nothing
)

"""
    SatMetadata(args...)

Default SatMetadata constructor for single floating point precision.
"""
SatMetadata(args...) = SatMetadata{Float32}(args...)

"""
    SatMetadata{T}(meta::SatMetadata) where T

External SatMetadata constructor for floating point conversions.
"""
SatMetadata{T}(meta::SatMetadata) where T = SatMetadata{T}(meta.files, meta.type,
  meta.date, meta.created, meta.loadtime, meta.remarks)


## Define structs related to track data

"""
# struct SatData

Satellite data with fields
- `data::DataFrame`
- `metadata::SatMetadata`

The `DataFrame` of `data` has columns for the satellite time and position and indices
for data files with additional information:
`
- `time::Vector{DateTime}`                        (time stamp of measurement)
- `lat::Vector{<:Union{Missing,T}}`   (latitude)
- `lon::Vector{<:Union{Missing,T}}`   (longitude)
- `fileindex::Vector{Int}`                        (index for filenames in `SatMetadata`)


# Instantiation

    function SatData{T}(
      folders::String...;
      type::Symbol=:undef,
      remarks=nothing
    ) where T -> struct SatData

Construct `SatData{T}` from any number of absolute or relative folder paths given as string.
SatData searches for hdf files in all folders recursively and determines the data type
(`CLay` or `CPro`) from the file names of the first 50 files unless the type is specified
with a Symbol `:CLay` or `:CPro` by the keyword argument `type`. Only one type of satellite
data can be stored in `SatData`. If `{T}` is omitted, it is set to `Float32`.

Alternatively, use handover all fields to the unmodified constructor, where basic
data validity checks will be performed:

    SatData{T}(data::DataFrame, metadata::SatMetadata) where T
"""
struct SatData{T} <: SatTrack{T}
  data::DataFrame
  metadata::SatMetadata{T}

  """ Unmodified constructor for `SatData` with basic checks for correct `data`"""
  function SatData{T}(data::DataFrame, metadata::SatMetadata) where T
    # Ensure floats of correct precision
    convertFloats!(data, T)
    # Check for correct column names and data types
    standardnames = ["time", "lat", "lon", "fileindex"]
    standardtypes = [Vector{DateTime}, Vector{T}, Vector{T}, Vector{Int}]
    bounds = (:time => (DateTime(2006), Dates.now()), :lat => (-90,90),
      :lon => (-180,180), :fileindex => (1, Inf))
    checkcols!(data, standardnames, standardtypes, bounds, "CLay")
    new{T}(data, metadata)
  end #constructor 1 SatData

  """
  Modified constructor creating the database from mat files in the given `folders`
  or any subfolder using the floating point precision given by `Float`. The sat data
  `type` is determined from the first 50 files in the database unless directly
  specified `type`. Any `remarks` can be added to the metadata.
  """
  function SatData{T}(
    folders::String...;
    type::Symbol=:undef,
    savedir::Union{String,Bool}="abs",
    remarks=nothing
  ) where T
    tstart = Dates.now()
    # Scan folders for HDF4 files
    files = String[]
    for folder in folders
      try findfiles!(files, folder, ".hdf")
      catch
        @warn "read error; data skipped" folder
      end
    end
    files = convertdir.(files, savedir)
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
    isempty(satfiles) && return SatData{T}()
    # Start MATLAB session
    ms = mat.MSession()
    # Initialise arrays
    utc = Vector{Vector{DateTime}}(undef, length(satfiles))
    lat = Vector{Vector{T}}(undef, length(satfiles))
    lon = Vector{Vector{T}}(undef, length(satfiles))
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
        DateTime[], T[], T[], Int[]
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
    SatData{T}(satdata, SatMetadata{T}(satfiles, (start=tmin, stop=tmax), loadtime, remarks=remarks))
  end #constructor 2 SatData
end #struct SatData


"""
    SatData{T}() where T

External constructor for empty `SatData` with floating point precision `T`.
"""
SatData{T}() where T = SatData{T}(
  DataFrame(time = DateTime[], lat = T[], lon=T[], fileindex=Int[]),
  SatMetadata{T}()
)

"""
    SatData(args...; kwargs...)

Default `Float32` constructor for `SatData`.
"""
SatData(args...; kwargs...) = SatData{Float32}(args...; kwargs...)

"""
    SatData{T}(sat::SatData) where T

External constructor for floating point precision conversions.
"""
SatData{T}(sat::SatData) where T = SatData{T}(
  DataFrame(
    time = sat.data.time,
    lat = T.(sat.data.lat),
    lon = T.(sat.data.lon),
    fileindex = sat.data.fileindex
  ), SatMetadata{T}(sat.metadata)
)

"""
    SatTrack{T}(args...; kwargs...) where T

Alias constructor for `SatData`.
"""
SatTrack{T}(args...; kwargs...) where T = SatData{T}(args...; kwargs...)

"""
    SatTrack(args...; kwargs...)

Alias constructor for default `Float32` `SatData`.
"""
SatTrack(args...; kwargs...) = SatData{Float32}(args...; kwargs...)


## ObservationSet

"""
# struct CLay

CALIOP cloud layer `data` stored in a `DataFrame` with columns:
- `time::Vector{DateTime}` (time index)
- `lat::Vector{AbstractFloat}` (latitude position of current time index)
- `lon::Vector{AbstractFloat}` (longitude position of current time index)
- `layer::Vector{NamedTuple{(:top,:base),Tuple{Vector{<:AbstractFloat},Vector{<:AbstractFloat}}}}`
  (layer top/base heights in meters)
- `atmos_state::Vector{Vector{Symbol}}` (symbols describing the atmospheric conditions at the intersection)
- `OD::Vector{<:Vector{<:AbstractFloat}}` (layer optical depth)
- `IWP::Vector{<:Vector{<:Union{Missing,<:AbstractFloat}}}` (layer ice water path)
- `Ttop::Vector{<:Vector{<:AbstractFloat}}` (layer top temperature)
- `Htropo::Vector{<:AbstractFloat}` (tropopause height at current time index)
- `night::BitVector` (flag for nights (`true`))
- `averaging::Vector{<:Vector{Int}}` (horizontal averaging in km)

# Instantiation

    function CLay{T}(
      ms::mat.MSession,
      files::Vector{String},
      timespan::NamedTuple{(:min,:max), Tuple{DateTime,DateTime}},
      lidarrange::Tuple{Real,Real}=(15_000,-Inf),
      altmin::Real=5000
    ) where T

Construct `CLay` from a list of file names (including directories) and a running
MATLAB session `ms` and save data, if layers are within the bounds
of `lidarrange` and above flight `altmin` threshold and time is within `timesapn`.
If `T<:AbstractFloat` is not set, `Float32` will be used as default precision.

Or construct `CLay` by directly handing over the `DataFrame` where the names, order,
and types of each columns are checked and attempted to correct:

    CLay{T}(data::DataFrame) where T -> struct CLay
"""
struct CLay{T} <: ObservationSet{T}
  data::DataFrame

  """ Unmodified constructor for `CLay` """
  function CLay{T}(data::DataFrame) where T
    # Ensure floats of correct precision
    convertFloats!(data, T)
    # Column checks and warnings
    standardnames = ["time", "lat", "lon", "layer_top", "layer_base", "atmos_state",
      "OD", "IWP", "Ttop", "Htropo", "night", "averaging"]
    standardtypes = [Vector{DateTime}, Vector{T}, Vector{T},
      Vector{<:Vector{<:T}}, Vector{<:Vector{<:T}},
      Vector{<:Vector{Symbol}}, Vector{<:Vector{<:T}},
      Vector{<:Vector{<:Union{Missing,<:T}}}, Vector{<:Vector{<:T}},
      Vector{<:T}, BitVector, Vector{<:Vector{<:Int}}]
    bounds = (:lat => (-90,90), :lon => (-180,180))
    checkcols!(data, standardnames, standardtypes, bounds, "CLay")
    new{T}(data)
  end #constructor 1 CLay

  """
  Modified constructor of `CLay` reading data from hdf `files` using MATLAB session `ms`
  in the `lidarrange` (top to bottom), if data is above `altmin`.
  """
  function CLay{T}(
    ms::mat.MSession,
    files::Vector{String},
    timespan::NamedTuple{(:min,:max), Tuple{DateTime,DateTime}},
    lidarrange::Tuple{Real,Real}=(15_000,-Inf),
    altmin::Real=5000,
    saveobs::Union{String,Bool}="abs"
  ) where T
    # Return default empty struct if files are empty
    (isempty(files) || saveobs === false || isempty(saveobs)) && return CLay{T}()
    # Initialise arrays
    # essential data
    utc = Vector{Vector{DateTime}}(undef, length(files))
    lat = Vector{Vector{T}}(undef, length(files))
    lon = Vector{Vector{T}}(undef, length(files))
    # non-essential data
    LayTop = Vector{Vector{Vector{T}}}(undef, length(files))
    LayBase = Vector{Vector{Vector{T}}}(undef, length(files))
    Atmosph = Vector{Vector{Vector{Symbol}}}(undef,length(files))
    OD = Vector{Vector{Vector{T}}}(undef,length(files))
    IWP = Vector{Vector{Vector{Union{Missing,T}}}}(undef,length(files))
    Ttop = Vector{Vector{Vector{T}}}(undef,length(files))
    Htropo = Vector{Vector{T}}(undef, length(files))
    night = Vector{BitVector}(undef, length(files))
    averaging = Vector{Vector{Vector{Int}}}(undef,length(files))
    # Loop over files
    for (i, file) in enumerate(files)
      ## Retrieve cloud layer data; assumes faulty files are filtered by SatData
      # Extract time
      mat.put_variable(ms, :file, file)
      mat.eval_string(ms, "clear t\ntry\nt = hdfread(file, 'Profile_UTC_Time');\nend")
      utc[i] = convertUTC.(mat.jarray(mat.get_mvariable(ms, :t))[:,2])
      timeindex = findall(timespan.min .≤ utc[i] .≤ timespan.max)
      utc[i] = utc[i][timeindex]
      # Extract lat/lon
      mat.eval_string(ms, "clear longitude\ntry\nlongitude = hdfread(file, 'Longitude');\nend")
      lon[i] = mat.jarray(mat.get_mvariable(ms, :longitude))[:,2][timeindex]
      mat.eval_string(ms, "clear latitude\ntry\nlatitude = hdfread(file, 'Latitude');\nend")
      lat[i] = mat.jarray(mat.get_mvariable(ms, :latitude))[:,2][timeindex]
      # Save time converted to UTC and lat/lon
      # utc[i], lon[i], lat[i] = convertUTC.(t), longitude, latitude

      ## Extract layer top/base, layer features and optical depth from hdf files
      mat.eval_string(ms, "clear basealt\ntry\nbasealt = hdfread(file, 'Layer_Base_Altitude');\nend")
      Lbase = 1000mat.jarray(mat.get_mvariable(ms, :basealt))[timeindex,:]
      mat.eval_string(ms, "clear topalt\ntry\ntopalt = hdfread(file, 'Layer_Top_Altitude');\nend")
      Ltop = 1000mat.jarray(mat.get_mvariable(ms, :topalt))[timeindex,:]
      mat.eval_string(ms, "clear FCF\ntry\nFCF = hdfread(file, 'Feature_Classification_Flags');\nend")
      FCF = mat.jarray(mat.get_mvariable(ms, :FCF))[timeindex,:]
      mat.eval_string(ms, "clear FOD\ntry\nFOD = hdfread(file, 'Feature_Optical_Depth_532');\nend")
      FOD = mat.jarray(mat.get_mvariable(ms, :FOD))[timeindex,:]
      mat.eval_string(ms, "clear IWPath\ntry\nIWPath = hdfread(file, 'Ice_Water_Path');\nend")
      IWPath = mat.jarray(mat.get_mvariable(ms, :IWPath))[timeindex,:]
      mat.eval_string(ms, "clear LTT\ntry\nLTT = hdfread(file, 'Layer_Top_Temperature');\nend")
      LTT = mat.jarray(mat.get_mvariable(ms, :LTT))[timeindex,:]
      mat.eval_string(ms, "clear Htropo\ntry\nHtropo = hdfread(file, 'Tropopause_Height');\nend")
      Htropo[i] = 1000vec(mat.jarray(mat.get_mvariable(ms, :Htropo)))[timeindex]
      mat.eval_string(ms, "clear daynight\ntry\ndaynight = hdfread(file, 'Day_Night_Flag');\nend")
      night[i] = Bool.(vec(mat.jarray(mat.get_mvariable(ms, :daynight))))[timeindex]
      mat.eval_string(ms, "clear average\ntry\naverage = hdfread(file, 'Horizontal_Averaging');\nend")
      horav = mat.jarray(mat.get_mvariable(ms, :average))[timeindex,:]
      # Loop over data and convert to TrackMatcher format
      Lt = Vector{Vector{T}}(undef,length(utc[i]))
      Lb = Vector{Vector{T}}(undef,length(utc[i]))
      atm = Vector{Vector{Symbol}}(undef,length(utc[i]))
      optdepth = Vector{Vector{T}}(undef,length(utc[i]))
      icewater = Vector{Vector{Union{Missing,T}}}(undef,length(utc[i]))
      toptemp = Vector{Vector{T}}(undef,length(utc[i]))
      average = Vector{Vector{Int}}(undef,length(utc[i]))
      for n = 1:length(utc[i])
        l = findall((Lbase[n,:] .> 0) .& (Ltop[n,:] .> 0) .& (Lbase[n,:] .< lidarrange[1]) .&
          (Ltop[n,:] .> lidarrange[2]) .& (Ltop[n,:] .> altmin))
        Lt[n], Lb[n], atm[n], optdepth[n], toptemp[n], icewater[n], average[n] =
          if isempty(l)
            T[], T[], Symbol[], T[], T[], T[], Int[]
          else
            l = findall((Lbase[n,:] .> 0) .& (Ltop[n,:] .> 0) .& (Lbase[n,:] .< lidarrange[1]) .&
              (Ltop[n,:] .> lidarrange[2]))
            [Ltop[n, m] for m in l] , [Lbase[n, m] for m in l],
            [feature_classification(classification(FCF[n,m])...) for m in l],
            [FOD[n,m] for m in l],
            [LTT[n,m] for m in l],
            [IWPath[n,m] == -9999 ? missing : IWPath[n,m] for m in l],
            [horav[n,m] for m in l]
          end
      end # loop over time steps in current file
      LayTop[i], LayBase[i], Atmosph[i], OD[i], IWP[i], Ttop[i], averaging[i] =
        Lt, Lb, atm, optdepth, icewater, toptemp, average
    end #loop over files

    # Construct and standardise data
    data = DataFrame(time=[utc...;], lat=[lat...;], lon=[lon...;],
      layer_top=[LayTop...;], layer_base=[LayBase...;], atmos_state=[Atmosph...;],
      OD=[OD...;], IWP=[IWP...;], Ttop=[Ttop...;], Htropo = [Htropo...;],
      night = [night...;], averaging = [averaging...;])
    # Save time, lat/lon arrays in CLay struct
    new{T}(data)
  end #constructor 2 CLay
end #struct CLay


"""
    CLay{T}() where T

External constructor for empty `CLay` struct.
"""
function CLay{T}() where T
  data = DataFrame(time = DateTime[], lat = T[], lon = T[],
  layer_top = Vector{T}[], layer_base = Vector{T}[],
  atmos_state = Vector{Symbol}[], OD = Vector{T}[], IWP = Vector{T}[],
  Ttop = Vector{T}[], Htropo = T[], night = BitVector(), averaging = Vector{Int}[])
  CLay{T}(data)
end

"""
    CLay{T}(clay::CLay) where T

External `CLay` constructor for conversion of floating point precision.
"""
function CLay{T}(clay::CLay) where T
  convertFloats!(clay.data, T)
  CLay{T}(clay.data)
end

""" Default CLay constructor for Float32 """
CLay(args...; kwargs...) = CLay{Float32}(args...; kwargs...)


"""
# struct CPro

CALIOP cloud profile `data` stored in a `DataFrame` with columns:
- `time::Vector{DateTime}` (current time index)
- `lat::Vector{AbstractFloat}` (latitude coordinate for current time index)
- `lon::Vector{AbstractFloat}` (lonitude coordinate for current time index)
- `atmos_state::Vector{<:Vector{<:Union{Missing,Symbol}}}`
  (symbols describing the atmospheric conditions for every height level at current time index)
- `EC532::Vector{<:Vector{<:Union{Missing,AbstractFloat}}}`
  (extinction coefficient at 532nm at every height level in current time index)

# Instantiation

    function CPro{T}(
      ms::mat.MSession,
      files::Vector{String},
      timespan::NamedTuple{(:min,:max), Tuple{DateTime,DateTime}},
      lidarprofile::NamedTuple
    ) where T -> struct CPro

Construct `CPro` from a list of file names (including directories) and a running
MATLAB session `ms`. CPro data is only stored in the vicinity of intersections for
the designated `timespan`. Column data is stored height-resolved as defined by the
`lidarprofile`. If `T<:AbstractFloat` is not set, `Float32` will be used as
default precision.

Or construct `CPro` by directly handing over the `DataFrame` where the names, order,
and types of each columns are checked and attempted to correct:

    CPro{T}(data::DataFrame) where T -> struct CPro
"""
struct CPro{T} <: ObservationSet{T}
  data::DataFrame

  """ unmodified constructor """
  function CPro{T}(data::DataFrame) where T
    # Ensure floats of correct precision
    convertFloats!(data, T)
    # Column checks and warnings
    standardnames = ["time", "lat", "lon", "atmos_state", "EC532", "Htropo", "temp",
      "pressure", "rH", "IWC", "deltap", "CADscore", "night"]
    standardtypes = [Vector{DateTime}, Vector{<:T}, Vector{<:T},
      Vector{<:Vector{<:Union{Missing,Symbol}}}, Vector{<:Vector{<:Union{Missing,<:T}}},
      Vector{<:T}, Vector{<:Vector{<:Union{Missing,<:T}}},
      Vector{<:Vector{<:Union{Missing,<:T}}}, Vector{<:Vector{<:Union{Missing,<:T}}},
      Vector{<:Vector{<:Union{Missing,<:T}}}, Vector{<:Vector{<:Union{Missing,<:T}}},
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
  function CPro{T}(
    ms::mat.MSession,
    files::Vector{String},
    timespan::NamedTuple{(:min,:max), Tuple{DateTime,DateTime}},
    lidarprofile::NamedTuple,
    saveobs::Union{String,Bool}="abs"
  ) where T
    # Return default empty struct if files are empty
    (isempty(files) || saveobs === false || isempty(saveobs)) && return CPro{T}()
    # Initialise arrays
    # essential data
    utc = Vector{Vector{DateTime}}(undef, length(files))
    lat = Vector{Vector{T}}(undef, length(files))
    lon = Vector{Vector{T}}(undef, length(files))
    fcf = Vector{Vector{Vector{<:Union{Missing,UInt16}}}}(undef, length(files))
    # non-essential data
    ec532 = Vector{Vector{Vector{<:Union{Missing,T}}}}(undef, length(files))
    Htropo = Vector{Vector{T}}(undef, length(files))
    temp = Vector{Vector{Vector{<:Union{Missing,T}}}}(undef, length(files))
    pres = Vector{Vector{Vector{<:Union{Missing,T}}}}(undef, length(files))
    rH = Vector{Vector{Vector{<:Union{Missing,T}}}}(undef, length(files))
    iwc = Vector{Vector{Vector{<:Union{Missing,T}}}}(undef, length(files))
    deltap = Vector{Vector{Vector{<:Union{Missing,T}}}}(undef, length(files))
    cad = Vector{Vector{Vector{<:Union{Missing,Int8}}}}(undef, length(files))
    night = Vector{BitVector}(undef, length(files))
    # Loop over files with cloud profile data
    for (i, file) in enumerate(files)
      ## Retrieve cloud profile data; assumes faulty files are filtered by SatData
      # Extract time
      mat.put_variable(ms, :file, file)
      mat.eval_string(ms, "clear t\ntry\nt = hdfread(file, 'Profile_UTC_Time');\nend")
      utc[i] = convertUTC.(mat.jarray(mat.get_mvariable(ms, :t))[:,2])
      timeindex = findall(timespan.min .≤ utc[i] .≤ timespan.max)
      utc[i] = utc[i][timeindex]
      # Extract lat/lon
      mat.eval_string(ms, "clear longitude\ntry\nlongitude = hdfread(file, 'Longitude');\nend")
      lon[i] = mat.jarray(mat.get_mvariable(ms, :longitude))[:,2][timeindex]
      mat.eval_string(ms, "clear latitude\ntry\nlatitude = hdfread(file, 'Latitude');\nend")
      lat[i] = mat.jarray(mat.get_mvariable(ms, :latitude))[:,2][timeindex]
      fcf[i] = get_lidarcolumn(UInt16, ms, "Atmospheric_Volume_Description", lidarprofile,
        coarse=false)[timeindex]
      # Extract non-essential data
      ec532[i] = get_lidarcolumn(T, ms, "Extinction_Coefficient_532", lidarprofile,
        missingvalues = -9999)[timeindex]
      mat.eval_string(ms, "clear Htropo\ntry\nHtropo = hdfread(file, 'Tropopause_Height');\nend")
      Htropo[i] = 1000vec(mat.jarray(mat.get_mvariable(ms, :Htropo)))[timeindex]
      temp[i] = get_lidarcolumn(T, ms, "Temperature", lidarprofile, missingvalues = -9999)[timeindex]
      pres[i] = get_lidarcolumn(T, ms, "Pressure", lidarprofile, missingvalues = -9999)[timeindex]
      rH[i] = get_lidarcolumn(T, ms, "Relative_Humidity", lidarprofile, missingvalues = -9999)[timeindex]
      iwc[i] = get_lidarcolumn(T, ms, "Ice_Water_Content_Profile", lidarprofile,
        missingvalues = -9999)[timeindex]
      deltap[i] = get_lidarcolumn(T, ms, "Particulate_Depalarization_Ratio_Profile_532",
        lidarprofile, missingvalues = -9999)[timeindex]
      cad[i] = get_lidarcolumn(Int8, ms, "CAD_Score", lidarprofile, coarse=false,
        missingvalues = -127)[timeindex]
      mat.eval_string(ms, "clear daynight\ntry\ndaynight = hdfread(file, 'Day_Night_Flag');\nend")
      night[i] = Bool.(vec(mat.jarray(mat.get_mvariable(ms, :daynight))))[timeindex]
    end #loop over files

    # Rearrange time vector and get time range
    utc = [utc...;]
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
    data = DataFrame(time=utc, lat=[lat...;], lon=[lon...;],
      atmos_state=avd, EC532=[ec532...;], Htropo = [Htropo...;],
      temp=[temp...;], pressure = [pres...;], rH = [rH...;],
      IWC = [iwc...;], deltap = [deltap...;],
      CADscore = [cad...;], night = [night...;])
    # Save time, lat/lon arrays, and feature classification flags (FCF) in CPro struct
    new{T}(data)
  end #constructor 2 CPro
end #struct CPro


"""
    CPro{T}() where T

External constructor for empty `CPro` struct.
"""
CPro{T}() where T = CPro{T}(DataFrame(time = DateTime[], lat = T[], lon = T[],
  atmos_state = Vector{Symbol}[], EC532 = Vector{T}[], Htropo = T[], temp = Vector{T}[],
  pressure = Vector{T}[], rH = Vector{T}[], IWC = Vector{T}[],
  deltap = Vector{T}[], CADscore = Vector{Int8}[], night = BitVector()))

"""
    CPro{T}(cpro::CPro) where T

External `CPro` constructor for conversion of floating point precision.
"""
function CPro{T}(cpro::CPro) where T
  convertFloats!(cpro.data, T)
  CPro{T}(cpro.data)
end

"""
    CPro(args...; kwargs...)

Default CPro constructor for Float32.
"""
CPro(args...; kwargs...) = CPro{Float32}(args...; kwargs...)
