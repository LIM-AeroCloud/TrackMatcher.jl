## Types related to Intersection

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
struct XMetadata{T} <: Intersection{T}
  maxtimediff::Int
  stepwidth::T
  Xradius::T
  expdist::T
  lidarrange::NamedTuple{(:top,:bottom),Tuple{T,T}}
  lidarprofile::NamedTuple
  sattype::Symbol
  satdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  altmin::T
  flightdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  created::Union{DateTime,ZonedDateTime}
  loadtime::Dates.CompoundPeriod
  remarks

  function XMetadata{T}(
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
    remarks
  ) where T
    new{T}(maxtimediff, stepwidth, Xradius, expdist,
      (top=T(lidarrange[1]), bottom=T(lidarrange[2])), lidarprofile,
      sattype, satdates, altmin, flightdates, created, loadtime, remarks)
  end #constructor 1 XMetaData

  function XMetadata{T}(
    maxtimediff::Int,
    stepwidth::Real,
    Xradius::Real,
    expdist::Real,
    lidarrange::NamedTuple{(:top,:bottom),Tuple{T,T}},
    lidarprofile::NamedTuple,
    sattype::Symbol,
    satdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    altmin::Real,
    flightdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
    created::Union{DateTime,ZonedDateTime},
    loadtime::Dates.CompoundPeriod,
    remarks=nothing
  ) where T
    new{T}(maxtimediff, stepwidth, Xradius, expdist, lidarrange, lidarprofile,
      sattype, satdates, altmin, flightdates, created, loadtime, remarks)
  end #constructor 2 XMetaData
end #struct XMetaData

""" Default XMetadata constructor for single floating point precision """
XMetadata(args...) = XMetadata{Float32}(args...)

""" External XMetadata constructor for floating point conversions """
XMetadata{T}(meta::XMetadata) where T = XMetadata{T}(
  meta.maxtimediff,
  T(meta.stepwidth),
  T(meta.Xradius),
  T(meta.expdist),
  (top = T(meta.lidarrange.top), bottom = T(meta.lidarrange.bottom)),
  (coarse = T.(meta.lidarprofile.coarse), fine = T.(meta.lidarprofile.fine),
    ibottom = meta.lidarprofile.ibottom, itop = meta.lidarprofile.itop,
    i30 = meta.lidarprofile.i30),
  meta.sattype,
  meta.satdates,
  T(meta.altmin),
  meta.flightdates,
  meta.created,
  meta.loadtime,
  meta.remarks
)


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
      flights::FlightSet,
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

Construct `Intersection` from the preloaded `FlightSet` and `SatData` with the option
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
struct XData{T} <: Intersection{T}
  data::DataFrame
  tracked::DataFrame
  accuracy::DataFrame
  metadata::XMetadata


  """ Unmodified constructor for `Intersection` """
  function XData{T}(
    data::DataFrame,
    tracked::DataFrame,
    accuracy::DataFrame,
    metadata::XMetadata{T}
  ) where T
    # Ensure floats of correct precision
    convertFloats!(data, T)
    convertFloats!(accuracy, T)
    tracked.flight = FlightData{T}.(tracked.flight)
    tracked.CPro = CPro{T}.(tracked.CPro)
    tracked.CLay = CLay{T}.(tracked.CLay)
    # Check data
    standardnames = ["id", "lat", "lon", "alt", "tdiff", "tflight", "tsat", "atmos_state"]
    standardtypes = [Vector{String}, Vector{<:T}, Vector{<:T}, Vector{<:T},
      Vector{Dates.CompoundPeriod}, Vector{DateTime}, Vector{DateTime},
      Vector{<:Union{Missing,Symbol}}]
    bounds = (:lat => (-90,90), :lon => (-180,180), :alt => (0, Inf))
    checkcols!(data, standardnames, standardtypes, bounds, "Intersection.data")
    # Check tracked (measured data)
    standardnames = ["id", "flight", "CPro", "CLay"]
    standardtypes = [Vector{String}, Vector{FlightData{T}}, Vector{CPro{T}}, Vector{CLay{T}}]
    bounds = ()
    checkcols!(tracked, standardnames, standardtypes, bounds, "Intersection.tracked",
      essentialcols = [1])
    # Check accuracy
    standardnames = ["id", "intersection", "flightcoord", "satcoord", "flighttime", "sattime"]
    standardtypes = [Vector{String}, Vector{<:T}, Vector{<:T}, Vector{<:T},
      Vector{Dates.CompoundPeriod}, Vector{Dates.CompoundPeriod}]
    bounds = ()
    checkcols!(accuracy, standardnames, standardtypes, bounds, "Intersection.accuracy",
      essentialcols = [1])
    new{T}(data, tracked, accuracy, metadata)
  end #constructor 1 XData


  """ Modified constructor with some automated calculations of the flight intersection data. """
  function XData{T}(
    tracks::PrimarySet,
    sat::SatData,
    savesecondsattype::Bool=false;
    maxtimediff::Int=30,
    primspan::Int=0,
    secspan::Int=15,
    lidarrange::Tuple{Real,Real}=(15_000,-Inf),
    stepwidth::Real=0.01,
    Xradius::Real=20_000,
    expdist::Real=Inf,
    remarks=nothing
  ) where T
    # Initialise DataFrames with Intersection data and monitor start time
    tstart = Dates.now()
    Xdata = DataFrame(id=String[], lat=T[], lon=T[], alt=T[],
      tdiff=Dates.CompoundPeriod[], tflight = DateTime[],
      tsat = DateTime[], atmos_state = Union{Missing,Symbol}[])
    tracked = DataFrame(id=String[], flight=FlightData{T}[], CPro=CPro{T}[], CLay=CLay{T}[])
    accuracy = DataFrame(id=String[], intersection=T[], flightcoord=T[],
      satcoord=T[], flighttime=Dates.CompoundPeriod[], sattime=Dates.CompoundPeriod[])
    # Combine all flight datasets and find intersections
    trackdata = tracks isa FlightSet ?
      [[getfield(tracks, f) for f in fieldnames(FlightSet)[1:end-1]]...;] : tracks.tracks
    # Get lidar altitude levels
    lidarprofile = get_lidarheights(lidarrange, T)
    # New MATLAB session
    ms = mat.MSession()
    # Loop over data from different datasets and interpolate track data and time, throw error on failure
    prog = pm.Progress(length(trackdata), "find intersections...")
    for (i, track) in enumerate(trackdata)
      # Get dataset source and ID
      dataset = track isa FlightTrack ? trackdata[i].metadata.source : "CloudTrack"
      ID = track isa FlightTrack ? trackdata[i].metadata.dbID : trackdata[i].metadata.ID
      try
        # Find sat tracks in the vicinity of flight tracks, where intersections are possible
        overlap = findoverlap(track, sat, maxtimediff, ID)
        if isempty(overlap)
          pm.next!(prog, showvalues = [(:hits, length(Xdata.id)),
            (:featured, length(Xdata.id[.!ismissing.(Xdata.atmos_state) .&
            (Xdata.atmos_state .≠ :no_signal) .& (Xdata.atmos_state .≠ :clear)]))])
          continue
        end
        # Interpolate trajectories with PCHIP method
        primtracks = interpolate_trackdata(track)
        sectracks = interpolate_satdata(sat, overlap, trackdata[i].metadata.useLON)
        # Calculate intersections and store data and metadata in DataFrames
        currdata, currtrack, curraccuracy = find_intersections(ms, track,
          primtracks, tracks.metadata.altmin, sat, sectracks, dataset, ID, maxtimediff,
          stepwidth, Xradius, lidarprofile, lidarrange, primspan, secspan,
          expdist, savesecondsattype, T)
        append!(Xdata, currdata); append!(tracked, currtrack)
        append!(accuracy, curraccuracy)
      catch err
        @debug begin
          @show ID
          rethrow(err)
        end
        # Issue warning on failure of interpolating track or time data
        @warn("Track data and/or time could not be interpolated. Data ignored.",
          dataset, ID)
      end
      # Monitor progress for progress bar
      pm.next!(prog, showvalues = [(:hits, length(Xdata.id)),
        (:featured, length(Xdata.id[.!ismissing.(Xdata.atmos_state) .&
        (Xdata.atmos_state .≠ :no_signal) .& (Xdata.atmos_state .≠ :clear)]))])
    end #loop over flights
    pm.finish!(prog)
    # Close MATLAB session after looping over all data
    mat.close(ms)
    # Calculate load time
    tend = Dates.now()
    tc = tz.ZonedDateTime(tend, tz.localzone())
    loadtime = Dates.canonicalize(Dates.CompoundPeriod(tend - tstart))
    # Return Intersections after completion
    @info string("Intersection data ($(length(Xdata[!,1])) matches) loaded in ",
      "$(join(loadtime.periods[1:min(2,length(loadtime.periods))], ", ")) to",
      "\n▪ data\n▪ tracked\n▪ accuracy\n▪ metadata")
    new{T}(Xdata, tracked, accuracy, XMetadata{T}(maxtimediff, stepwidth, Xradius,
      expdist, lidarrange, lidarprofile, sat.metadata.type, sat.metadata.date,
      tracks.metadata.altmin, tracks.metadata.date, tc, loadtime, remarks))
  end #constructor 2 XData
end #struct XData


""" External constructor for conversion of floating point precision """
XData{T}(X::XData) where T = XData{T}(X.data, X.tracked, X.accuracy, XMetadata{T}(X.metadata))

""" Default FlightData constructor for Float32 """
function XData(
  tracks::PrimarySet{T1},
  sat::SatData{T2},
  savesecondsattype::Bool=false;
  kwargs...
) where {T1, T2}
  T = promote_type(T1, T2)
  tracks isa PrimarySet{T} || (tracks = PrimarySet{T}(tracks))
  sat isa SatData{T} || (sat = SatData{T}(sat))
  XData{T}(tracks, sat, savesecondsattype; kwargs...)
end
""" Default FlightTrack constructor for Float32 FlightData """
Intersection(args...; kwargs...) = XData(args...; kwargs...)

""" FlightTrack constructor for FlightData """
Intersection{T}(args...; kwargs...) where T = XData{T}(args...; kwargs...)



## ObservationSet

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
    altmin::Real=5000
  ) where T
    # Return default empty struct if files are empty
    isempty(files) && return CLay{T}()
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


""" External constructor for emtpy CLay struct """
function CLay{T}() where T
  data = DataFrame(time = DateTime[], lat = T[], lon = T[],
  layer_top = Vector{T}[], layer_base = Vector{T}[],
  atmos_state = Vector{Symbol}[], OD = Vector{T}[], IWP = Vector{T}[],
  Ttop = Vector{T}[], Htropo = T[], night = BitVector(), averaging = Vector{Int}[])
  CLay{T}(data)
end

""" External CLay constructor for conversion of floating point precision """
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
    lidarprofile::NamedTuple
  ) where T
    # Return default empty struct if files are empty
    isempty(files) && return CPro{T}()
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


""" External constructor for emtpy CPro struct """
CPro{T}() where T = CPro{T}(DataFrame(time = DateTime[], lat = T[], lon = T[],
  atmos_state = Vector{Symbol}[], EC532 = Vector{T}[], Htropo = T[], temp = Vector{T}[],
  pressure = Vector{T}[], rH = Vector{T}[], IWC = Vector{T}[],
  deltap = Vector{T}[], CADscore = Vector{Int8}[], night = BitVector()))

""" External CPro constructor for conversion of floating point precision """
function CPro{T}(cpro::CPro) where T
  convertFloats!(cpro.data, T)
  CPro{T}(cpro.data)
end

""" Default CPro constructor for Float32 """
CPro(args...; kwargs...) = CPro{Float32}(args...; kwargs...)


struct Data{T} <: DataSet{T}
  flight::Union{Nothing,FlightSet{T}}
  cloud::Union{Nothing,CloudSet{T}}
  sat::Union{Nothing,SatData{T}}
  intersection::Union{Nothing,NamedTuple{(:flight,:cloud), Tuple{XData{T},XData{T}}}}

  function Data{T}(
    flight::Union{Nothing,FlightSet{T}},
    cloud::Union{Nothing,CloudSet{T}},
    sat::Union{Nothing,SatData{T}},
    intersection::Union{Nothing,NamedTuple{(:flight,:cloud), Tuple{XData{T},XData{T}}}}
  ) where T
    new{T}(flight, cloud, sat, intersection)
  end

  function Data{T}(
    folders::Vector{<:Pair{String,<:Any}},
    savesecondsattype::Bool=false;
    sattype::Symbol=:undef,
    altmin::Real=5000,
    odelim::Union{Nothing,Char,String}=nothing,
    maxtimediff::Int=30,
    primspan::Int=0,
    secspan::Int=15,
    lidarrange::Tuple{Real,Real}=(15_000,-Inf),
    stepwidth::Real=0.01,
    Xradius::Real=20_000,
    expdist::Real=Inf,
    remarks::Vector{<:Pair{String,<:Any}}=Pair{String,Any}[]
  ) where T
    # Process function arguments that need to be distributed to several structs
    folders = init_dict(folders, String[])
    remarks = init_dict(remarks, nothing)

    # Load data
    flights = FlightSet{T}(;
      inventory = folders["inventory"],
      archive = folders["archive"],
      onlineData = folders["onlineData"],
      altmin, odelim, remarks=remarks["flights"]
    )
    @debug trim_vec!.([flights.inventory, flights.archive, flights.onlineData], 300)
    clouds = CloudSet{T}(folders["cloudtracks"]...; remarks = remarks["clouds"])
    sat = SatTrack{T}(
      folders["sat"]...;
      type = sattype,
      remarks = remarks["sat"]
    )

    # Calculate Intersections
    intersections = (
    flight=Intersection{T}(
      flights, sat, savesecondsattype;
      maxtimediff, primspan, secspan, lidarrange,
      stepwidth, Xradius, expdist,
      remarks = remarks["Xflight"]
    ),
    cloud = Intersection{T}(
      clouds, sat, savesecondsattype;
      maxtimediff, primspan, secspan, lidarrange,
      stepwidth, Xradius, expdist,
      remarks = remarks["Xflight"]
    ))
    new{T}(flights, clouds, sat, intersections)
  end #constructor 2 for Data
end


Data(args...; kwargs...) = Data{Float32}(args...; kwargs...)

Data{T}(data::Data) where T = Data{T}(
  FlightSet{T}(data.flight),
  CloudSet{T}(data.cloud),
  SatData{T}(data.sat),
  (flight = Intersection{T}(data.intersection.flight),
    cloud = Intersection{T}(data.intersection.cloud))
)

DataSet{T}(args...; kwargs...) where T = Data{T}(args...; kwargs...)

DataSet(args...; kwargs...) = DataSet{Float32}(args...; kwargs...)
