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

## granules
DataFrame with `filename`s, `start` and `stop` time.

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
struct SecondaryMetadata{T} <: SecondarySet{T}
  granules::DataFrame
  type::Symbol
  date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  created::Union{DateTime,ZonedDateTime}
  loadtime::Dates.CompoundPeriod
  remarks
end #struct SecondaryMetadata

"""
    SatMetadata{T}() where T

External constructor for empty `SatMetadata` with floating point precision `T`.
"""
SecondaryMetadata{T}() where T = SecondaryMetadata{T}(
  DataFrame(file = String[], tstart = DateTime[], tstop = DateTime[],
      latmin = T[], latmax = T[], elonmin = T[], elonmax = T[], wlonmin = T[],
      wlonmax = T[]), :undef, (start=Dates.now(), stop=Dates.now()),
  Dates.now(), Dates.CompoundPeriod(), nothing
)

"""
    SatMetadata(args...)

Default SatMetadata constructor for single floating point precision.
"""
SecondaryMetadata(args...) = SecondaryMetadata{Float32}(args...)

"""
    SatMetadata{T}(meta::SatMetadata) where T

External SatMetadata constructor for floating point conversions.
"""
SecondaryMetadata{T}(meta::SecondaryMetadata) where T = SecondaryMetadata{T}(
  DataFrame(file = meta.granules.file,
    tstart = meta.granules.tstart, tstop = meta.granules.tstop,
    elonmin = T.(meta.granules.elonmin), elonmax = T.(meta.granules.elonmax),
    wlonmin = T.(meta.granules.wlonmin), wlonmax = T.(meta.granules.wlonmax)),
  meta.type, meta.date, meta.created, meta.loadtime, meta.remarks
)


## Define structs related to track data

struct SatData{T} <: SatTrack{T}
  data::DataFrame

  """ Unmodified constructor for `SatData` with basic checks for correct `data`"""
  function SatData{T}(data::DataFrame) where T
    # Ensure floats of correct precision
    convertFloats!(data, T)
    # Check for correct column names and data types
    standardnames = ["time", "lat", "lon"]
    standardtypes = [Vector{DateTime}, Vector{T}, Vector{T}]
    bounds = (:time => (DateTime(2006), Dates.now()), :lat => (-90,90),
      :lon => (-180,180))
    checkcols!(data, standardnames, standardtypes, bounds, "CLay")
    new{T}(data)
  end #constructor 1 SatData

  """
  Modified constructor creating the database from mat files in the given `folders`
  or any subfolder using the floating point precision given by `Float`. The sat data
  `type` is determined from the first 50 files in the database unless directly
  specified `type`. Any `remarks` can be added to the metadata.
  """
  function SatData{T}(ms::mat.MSession, file::String) where T
    # Send file name to MATLAB
    mat.put_variable(ms, :file, file)
    # Read hdf file with MATLAB
    mat.eval_string(ms,
    """
    clear t longitude latitude
    try
      t = hdfread(file, 'Profile_UTC_Time');
      longitude = hdfread(file, 'Longitude');
      latitude = hdfread(file, 'Latitude');
    end
    data = struct('time', t(:, 2), 'lat', latitude(:, 2), 'lon', longitude(:, 2));
    """)
    # Retrieve data from MATLAB
    data = mat.jdict(mat.get_mvariable(ms, :data))
    # Convert time to UTC DateTime
    utc = convertUTC.(data["time"])
    # instantiate new struct
    new{T}(DataFrame(time = utc, lat = data["lat"], lon = data["lon"]))
  end #constructor 2 SatData
end #struct SatData


"""
    SatData{T}() where T

External constructor for empty `SatData` with floating point precision `T`.
"""
SatData{T}() where T = SatData{T}(
  DataFrame(time = DateTime[], lat = T[], lon=T[])
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
    lon = T.(sat.data.lon)
  )
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


"""
# struct SatSet

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

    function SaSet{T}(
      folders::String...;
      type::Symbol=:undef,
      remarks=nothing
    ) where T -> struct SatData

Construct `SatSet{T}` from any number of absolute or relative folder paths given as string.
SatData searches for hdf files in all folders recursively and determines the data type
(`CLay` or `CPro`) from the file names of the first 50 files unless the type is specified
with a Symbol `:CLay` or `:CPro` by the keyword argument `type`. Only one type of satellite
data can be stored in `SatData`. If `{T}` is omitted, it is set to `Float32`.

Alternatively, use handover all fields to the unmodified constructor, where basic
data validity checks will be performed:

    SatData{T}(data::DataFrame, metadata::SatMetadata) where T
"""
struct SatSet{T} <: SecondarySet{T}
  granules::Vector{SatData{T}}
  metadata::SecondaryMetadata{T}

  """ Unmodified constructor for `SatData` with basic checks for correct `data`"""
  function SatSet{T}(granules::Vector{SatData}, metadata::SecondaryMetadata{T}) where T
    # Ensure floats of correct precision
    convertFloats!.(granules, T)
    # Check for correct column names and data types
    standardnames = ["time", "lat", "lon"]
    new{T}(granules, metadata)
  end #constructor 1 SatData

  """
  Modified constructor creating the database from mat files in the given `folders`
  or any subfolder using the floating point precision given by `Float`. The sat data
  `type` is determined from the first 50 files in the database unless directly
  specified `type`. Any `remarks` can be added to the metadata.
  """
  function SatSet{T}(
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
    # Apply format for saved file paths
    files = convertdir.(files, savedir)
    # If type of satellite data is not defined, find it based on first 50 file names (~2 days)
    sattype = type == :CLay || type == :CPro ? string(type) :
      count(occursin.("CLay", files[1:min(length(files), 50)])) ≥
      count(occursin.("CPro", files[1:min(length(files), 50)])) ? "CLay" : "CPro"
    type == :undef && (type = Symbol(type, sattype))
    # Select files of type based on file name
    satfiles = files[occursin.(sattype, basename.(files))]
    wrongtype = length(files) - length(satfiles)
    wrongtype > 0 &&
      @warn "$wrongtype files with wrong satellite type detected; data skipped"
    # Create empty struct, if no data files were found
    isempty(satfiles) && return SatSet{T}()
    # Start MATLAB session
    ms = mat.MSession()
    # Initialise array of granules
    granules = SatData[]
    metadata = DataFrame(file = String[], tstart = DateTime[], tstop = DateTime[],
      latmin = T[], latmax = T[], elonmin = T[], elonmax = T[], wlonmin = T[],
      wlonmax = T[])
    # Loop over files
    prog = pm.Progress(length(satfiles), "load sat data...")
    for file in satfiles
      # Find files with cloud layer data
      try
        granule = SatData(ms, file)
        elonmin, elonmax = lonextrema(granule.data.lon, ≥)
        wlonmin, wlonmax = lonextrema(granule.data.lon, <)
        push!(metadata, (file=file,
          tstart = granule.data.time[1], tstop = granule.data.time[end],
          latmin = minimum(granule.data.lat), latmax = maximum(granule.data.lat),
          elonmin, elonmax, wlonmin, wlonmax)
        )
        push!(granules, granule)
      catch
        # Skip data on failure and warn
        @warn "read error in CALIPSO granule; data skipped"  granule = splitext(basename(file))[1]
      end
      # Monitor progress for progress bar
      pm.next!(prog, showvalues = [(:date,Dates.Date(splitdir(dirname(file))[2], "y_m_d"))])
    end #loop over files
    pm.finish!(prog)

    # Close MATLAB session
    mat.close(ms)
    # Save computing times.
    tend = Dates.now()
    tc = tz.ZonedDateTime(tend, tz.localzone())
    loadtime = Dates.canonicalize(Dates.CompoundPeriod(tend - tstart))
    tmin = minimum(granule.data.time[1] for granule in granules)
    tmax = maximum(granule.data.time[end] for granule in granules)

    # Instantiate new struct
    @info string("SatSet loaded in ",
      "$(join(loadtime.periods[1:min(2,length(loadtime.periods))], ", ")) to",
      "\n▪ granules ($(sum(size(granule.data, 1) for granule in granules)) data rows)\n  – time\n  – lat\n  – lon",
      "\n▪ metadata")
    new{T}(granules, SecondaryMetadata{T}(metadata, type, (start=tmin, stop=tmax),
      tc, loadtime, remarks))
  end #constructor 2 SatData
end #struct SatSet


"""
    SatSet{T}() where T

External constructor for empty `SatSet` with floating point precision `T`.
"""
SatSet{T}() where T = SatSet{T}(
  SatData{T}[],
  SecondaryMetadata{T}()
)

"""
    SatSet(args...; kwargs...)

Default `Float32` constructor for `SatSet`.
"""
SatSet(args...; kwargs...) = SatSet{Float32}(args...; kwargs...)

"""
    SatSet{T}(sat::SatSet) where T

External constructor for floating point precision conversions.
"""
SatSet{T}(sat::SatSet) where T = SatSet{T}(
  SatData{T}.(sat), SecondaryMetadata{T}(sat.metadata)
)


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
    timeindex::Vector{UnitRange{Int}},
    lidarrange::Tuple{Real,Real}=(15_000,-Inf),
    altmin::Real=5000,
    saveobs::Union{String,Bool}="abs"
  ) where T
    # Return default empty struct if files are empty
    (isempty(files) || saveobs === false || isempty(saveobs)) && CLay()
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
    averaging = Vector{Vector{Int}}(undef,length(files))
    # Loop over files
    for (i, file) in enumerate(files)
      # Send file name to MATLAB
      mat.put_variable(ms, :file, file)
      # Read hdf file with MATLAB
      mat.eval_string(ms,
        """
        clear t longitude latitude basealt topalt ...
          FCF FOD IWPath LTT Htropo daynight average
        try
          t = hdfread(file, 'Profile_UTC_Time');
          t = t(:, 2);
          longitude = hdfread(file, 'Longitude');
          longitude = longitude(:, 2);
          latitude = hdfread(file, 'Latitude');
          latitude = latitude(:, 2);
          basealt = hdfread(file, 'Layer_Base_Altitude');
          topalt = hdfread(file, 'Layer_Top_Altitude');
          FCF = hdfread(file, 'Feature_Classification_Flags');
          FOD = hdfread(file, 'Feature_Optical_Depth_532');
          IWPath = hdfread(file, 'Ice_Water_Path');
          LTT = hdfread(file, 'Layer_Top_Temperature');
          Htropo = hdfread(file, 'Tropopause_Height');
          daynight = hdfread(file, 'Day_Night_Flag');
          average = hdfread(file, 'Horizontal_Averaging');
        end

        data = struct('time', t, 'lat', latitude, 'lon', longitude, ...
          'LB', basealt, 'LT', topalt, 'FCF', FCF, 'FOD', FOD, 'IWP', IWPath, ...
          'LTT', LTT, 'Htropo', Htropo, 'daynight', daynight, 'average', average);
        """
      )

      # Get data from MATLAB
      data = mat.jdict(mat.get_mvariable(ms, :data))
      # Convert time to UTC DateTime
      utc[i] = convertUTC.(data["time"])
      utc[i] = utc[i][timeindex[i]]
      # Extract lat/lon
      lat[i], lon[i] = data["lat"][timeindex[i]], data["lon"][timeindex[i]]

      ## Extract non-essential information from satellite files
      Ltop, Lbase = 1000data["LT"][timeindex[i],:], 1000data["LB"][timeindex[i],:]
      FCF = data["FCF"][timeindex[i],:]
      FOD = data["FOD"][timeindex[i],:]
      IWPath = data["IWP"][timeindex[i],:]
      LTT = data["LTT"][timeindex[i],:]
      Htropo[i] = 1000data["Htropo"][timeindex[i]]
      night[i] = Bool.(data["daynight"][timeindex[i]])
      averaging[i] = data["average"][timeindex[i]]

      # Loop over data and convert to TrackMatcher format
      Lt = Vector{Vector{T}}(undef,length(utc[i]))
      Lb = Vector{Vector{T}}(undef,length(utc[i]))
      atm = Vector{Vector{Symbol}}(undef,length(utc[i]))
      optdepth = Vector{Vector{T}}(undef,length(utc[i]))
      icewater = Vector{Vector{Union{Missing,T}}}(undef,length(utc[i]))
      toptemp = Vector{Vector{T}}(undef,length(utc[i]))
      for n = 1:length(utc[i])
        l = findall((Lbase[n,:] .> 0) .& (Ltop[n,:] .> 0) .& (Lbase[n,:] .< lidarrange[1]) .&
          (Ltop[n,:] .> lidarrange[2]) .& (Ltop[n,:] .> altmin))
        Lt[n], Lb[n], atm[n], optdepth[n], toptemp[n], icewater[n] =
          if isempty(l)
            T[], T[], Symbol[], T[], T[], T[]
          else
            l = findall((Lbase[n,:] .> 0) .& (Ltop[n,:] .> 0) .& (Lbase[n,:] .< lidarrange[1]) .&
              (Ltop[n,:] .> lidarrange[2]))
            [Ltop[n, m] for m in l] , [Lbase[n, m] for m in l],
            [feature_classification(classification(FCF[n,m])...) for m in l],
            [FOD[n,m] for m in l],
            [LTT[n,m] for m in l],
            [IWPath[n,m] == -9999 ? missing : IWPath[n,m] for m in l]
          end
      end # loop over time steps in current file
      LayTop[i], LayBase[i], Atmosph[i], OD[i], IWP[i], Ttop[i] =
        Lt, Lb, atm, optdepth, icewater, toptemp
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
    new{T}(data)
  end #constructor 1 CPro

  """
  Modified constructor of `CPro` reading data from hdf `files` for all given `sattime` indices
  using MATLAB session `ms` and `lidarprofile` data, if data is above `altmin`.
  """
  function CPro{T}(
    ms::mat.MSession,
    files::Vector{String},
    timeindex::Vector{UnitRange{Int}},
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
      mat.put_variable(ms, :file, file)
      mat.eval_string(ms,
        """
        clear t longitude latitude
        try
          t = hdfread(file, 'Profile_UTC_Time');
          t = t(:, 2);
          longitude = hdfread(file, 'Longitude');
          longitude = longitude(:, 2);
          latitude = hdfread(file, 'Latitude');
          latitude = latitude(:, 2);
          FCF = hdfread(file, 'Atmospheric_Volume_Description');
          EC532 = hdfread(file, 'Extinction_Coefficient_532');
          Htropo = hdfread(file, 'Tropopause_Height');
          temp = hdfread(file, 'Temperature');
          pres = hdfread(file, 'Pressure');
          rH = hdfread(file, 'Relative_Humidity');
          IWC = hdfread(file, 'Ice_Water_Content_Profile');
          deltap = hdfread(file, 'Particulate_Depolarization_Ratio_Profile_532');
          CAD = hdfread(file, 'CAD_Score');
          daynight = hdfread(file, 'Day_Night_Flag');
        end

        data = struct('time', t, 'lat', latitude, 'lon', longitude, ...
          'FCF', FCF, 'EC532', EC532, 'Htropo', Htropo, 'temp', temp, 'pres', pres, ...
          'rH', rH, 'IWC', IWC, 'deltap', deltap, 'CAD', CAD, 'daynight', daynight);
        """
      )
      # Get data from MATLAB
      data = mat.jdict(mat.get_mvariable(ms, :data))
      # Convert time to UTC DateTime
      utc[i] = convertUTC.(data["time"])
      utc[i] = utc[i][timeindex[i]]
      # Extract essential data
      lat[i], lon[i] = data["lat"][timeindex[i]], data["lon"][timeindex[i]]
      fcf[i] = get_lidarcolumn(UInt16, data["FCF"], lidarprofile, coarse=false)[timeindex[i]]
      # Extract non-essential data
      ec532[i] = get_lidarcolumn(T, data["EC532"], lidarprofile, missingvalues = -9999)[timeindex[i]]
      Htropo[i] = 1000data["Htropo"][timeindex[i]]
      temp[i] = get_lidarcolumn(T, data["temp"], lidarprofile, missingvalues = -9999)[timeindex[i]]
      pres[i] = get_lidarcolumn(T, data["pres"], lidarprofile, missingvalues = -9999)[timeindex[i]]
      rH[i] = get_lidarcolumn(T, data["rH"], lidarprofile, missingvalues = -9999)[timeindex[i]]
      iwc[i] = get_lidarcolumn(T, data["IWC"], lidarprofile, missingvalues = -9999)[timeindex[i]]
      deltap[i] = get_lidarcolumn(T, data["deltap"], lidarprofile, missingvalues = -9999)[timeindex[i]]
      cad[i] = get_lidarcolumn(Int8, data["CAD"], lidarprofile, coarse=false,
        missingvalues = -127)[timeindex[i]]
      night[i] = Bool.(data["daynight"])[timeindex[i]]
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
    # Instantiate new struct
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
