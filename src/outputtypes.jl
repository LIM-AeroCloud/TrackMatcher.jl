## Types related to Intersection

"""
# struct XMetadata

Immutable struct with additional information of intersection data:

- `maxtimediff::Int`: maximum time difference in minutes allowed between
  satellite overpass and aircraft passing at intersection
- `stepwidth::T`: stepwidth (in meters, partially internally converted to degrees at equator)
  used to interpolate track data
- `Xradius::T`: radius in meters around an intersection in which further intersections
  will be removed as duplicates due to the interpolation algorithm
- `expdist::T`: threshold for allowed minimum distance of a calculated intersection
  to the nearest measured track point
- `lidarrange::NamedTuple{(:top,:bottom),Tuple{Real,Real}}`: user defined level thresholds
  for top/bottom heights for which CALIPSO data is considered
- `lidarprofile::NamedTuple`: CALIPSO lidar height levels used in the current dataset
- `sattype::Symbol`: cloud layer or profile data types used for the satellite date
- `satdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}`: date range of
  satellite data
- `altmin::T`: minimum threshold for which flight data is considered
- `flightdates::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}`: date range
  of flight data
- `created::Union{DateTime,ZonedDateTime}`: time of creation of database
- `loadtime::Dates.CompoundPeriod`: time it took to find intersections and load it to the struct
- `remarks::Any`: any additional data or comments that can be attached to the database

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

"""
    XMetadata(args...)

Default `XMetadata` constructor for single floating point precision.
"""
XMetadata(args...) = XMetadata{Float32}(args...)

"""
    XMetadata{T}(meta::XMetadata) where T

External `XMetadata` constructor for floating point conversions.
"""
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
# struct XData

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
- `atmos_state::Vector{<:Union{Missing,Symbol}}`: atmospheric conditions at intersection


## tracked

Original track data in the vicinity of the intersection.

`DataFrame` columns are:
- `id: Vector{String}`: unique intersection identifier
- `flight::Vector{FlightData}`: `FlightData` in the vicinity of the intersection
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
    ) where T -> struct Intersection

Construct `XData` from the preloaded `PrimarySet` and `SatData` with the option
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
- `remarks=nothing`: any data or remarks attached to the metadata

If `T<:AbstractFloat` is omitted, the default `Float32` precision is used.

Or construct `Intersection` by directly handing over the different `DataFrame`s where the names, order,
and types of each columns are checked and attempted to correct, together with the metadata:

    function XData{T}(
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
    savedir::Union{String,Bool}="abs",
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
          expdist, savedir, savesecondsattype, T)
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


"""
    XData{T}(X::XData) where T

External constructor for conversion of floating point precision.
"""
XData{T}(X::XData) where T = XData{T}(X.data, X.tracked, X.accuracy, XMetadata{T}(X.metadata))

"""
    function XData(
      tracks::PrimarySet{T1},
      sat::SatData{T2},
      savesecondsattype::Bool=false;
      kwargs...
    ) where {T1, T2}

Default `XData` constructor where `T1` and `T2` floating point precisions of the flight
and sat data are promoted to the higher precision.
"""
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

"""
    Intersection(args...; kwargs...)

Alias constructor for default `XData` construction with promoted types of flight
and sat data.
"""
Intersection(args...; kwargs...) = XData(args...; kwargs...)

"""
    Intersection{T}(args...; kwargs...) where T

Alias constructor for `XData{T}`.
"""
Intersection{T}(args...; kwargs...) where T = XData{T}(args...; kwargs...)


## Overall combined data for one-step data loading and intersection calculation

"""
# struct MeasuredData

Store all relevant primary and secondary track data depending on the primary source to fields:

- `flight::Union{Nothing,FlightSet{T}}`
- `cloud::Union{Nothing,CloudSet{T}}`
- `sat::Union{Nothing,SatData{T}}`

# Instantiate

Construct `Data` from the individual fields or use a modified constructor to load
all necessary data from the file names given in a vector of pairs with the following
`String` keywords for the different databases:
- `"inventory"`: VOLPE AEDT database
- `"archive"`: FlightAware commercial data
- `"onlineData"`: FlightAware web content
- `"cloudtracks"`: cloud track data
- `"sat"`: CALIPSO satellite track data

Furthermore, the keyword arguments for `FlightSet`, `CloudSet`, and `SatData`are passed on.

- `sattype::Symbol=:unde`,
- `altmin::Real=5000`
- `odelim::Union{Nothing,Char,String}=nothing`
- `remarks::Vector{<:Pair{String,<:Any}}=Pair{String,Any}[]`

For remarks, a vector of pairs is used again to differentiate
between the different datasets:

- `"flights"`: flight data (all sources within `FlightSet`)
- `"clouds"`: cloud track data
- `"sat"`: CALIPSO satellite data
- `"Xflight"`: calculated intersections using flight data as primary source
- `"Xcloud"`: calculated intersections using cloud data as primary source
"""
struct MeasuredData{T} <: MeasuredSet{T}
  flight::Union{Nothing,FlightSet{T}}
  cloud::Union{Nothing,CloudSet{T}}
  sat::Union{Nothing,SatData{T}}

  """ unmodified constructor for `MeasuredData` """
  function MeasuredData{T}(
    flight::Union{Nothing,FlightSet{T}},
    cloud::Union{Nothing,CloudSet{T}},
    sat::Union{Nothing,SatData{T}},
  ) where T
    new{T}(flight, cloud, sat)
  end # unmodified constructor for MeasuredData

  """ modified constructor for `MeasuredData` """
  function MeasuredData{T}(
    folders::Vector{<:Pair{String,<:Any}};
    sattype::Symbol=:undef,
    altmin::Real=5000,
    odelim::Union{Nothing,Char,String}=nothing,
    savedir::Union{String,Bool}="abs",
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
      altmin, odelim, savedir, remarks=remarks["flights"]
    )
    @debug trim_vec!.([flights.inventory, flights.archive, flights.onlineData], 300)
    clouds = CloudSet{T}(folders["cloudtracks"]...; savedir, remarks = remarks["clouds"])
    sat = SatTrack{T}(
      folders["sat"]...;
      type = sattype,
      savedir,
      remarks = remarks["sat"]
    )

    # Instantiate
    new{T}(flights, clouds, sat)
  end # modified constructor for MeasuredData
end #struct MeasuredData

"""
    MeasuredData(args...; kwargs...)

Default `MeasuredData` constructor for single floating point precision.
"""
MeasuredData(args...; kwargs...) = MeasuredData{Float32}(args...; kwargs...)

"""
    MeasuredData{T}(data::MeasuredData) where T

Constructor for floating point conversions.
"""
MeasuredData{T}(data::MeasuredData) where T = MeasuredData{T}(
  FlightSet{T}(data.flight),
  CloudSet{T}(data.cloud),
  SatData{T}(data.sat)
)

"""
    MeasuredSet{T}(args...; kwargs...) where T

Alias constructor for `MeasuredData` of type `T`.
"""
MeasuredSet{T}(args...; kwargs...) where T = MeasuredData{T}(args...; kwargs...)

"""
    MeasuredSet(args...; kwargs...)

Alias default constructor for `MeasuredData` of type `Float32`.
"""
MeasuredSet(args...; kwargs...) = MeasuredSet{Float32}(args...; kwargs...)


"""
# struct Data

Store all relevant primary and secondary track data and calculated intersections
depending on the primary source to fields:

- `trackdata::MeasuredData{T}`
- `intersection::NamedTuple{(:flight,:cloud), Tuple{XData{T},XData{T}}}`

# Instantiate

Construct `Data` from the individual fields or use a modified constructor to load
all necessary data from the file names given in a vector of pairs with the following
`String` keywords for the different databases:
- `"inventory"`: VOLPE AEDT database
- `"archive"`: FlightAware commercial data
- `"onlineData"`: FlightAware web content
- `"cloudtracks"`: cloud track data
- `"sat"`: CALIPSO satellite track data

Furthermore, the keyword arguments for `FlightSet`, `CloudSet`, `SatData`, and
`XData` are passed on.

- `savesecondsattype::Bool=false`
- `sattype::Symbol=:unde`,
- `altmin::Real=5000`
- `odelim::Union{Nothing,Char,String}=nothing`
- `maxtimediff::Int=30`
- `primspan::Int=0`
- `secspan::Int=15`
- `lidarrange::Tuple{Real,Real}=(15_000,-Inf)`
- `stepwidth::Real=0.01`
- `Xradius::Real=20_000`
- `expdist::Real=Inf`
- `remarks::Vector{<:Pair{String,<:Any}}=Pair{String,Any}[]`

For remarks, a vector of pairs is used again to differentiate
between the different datasets:

- `"flight"`: flight data (all sources within `FlightSet`)
- `"cloud"`: cloud track data
- `"sat"`: CALIPSO satellite data
- `"Xflight"`: calculated intersections using flight data as primary source
- `"Xcloud"`: calculated intersections using cloud data as primary source
"""
struct Data{T} <: DataSet{T}
  trackdata::MeasuredData{T}
  intersection::NamedTuple{(:flight,:cloud), Tuple{XData{T},XData{T}}}

  """ unmodified constructor for `Data` """
  function Data{T}(
    trackdata::MeasuredData{T},
    intersection::NamedTuple{(:flight,:cloud), Tuple{XData{T},XData{T}}}
  ) where T
    new{T}(trackdata, intersection)
  end # unmodified constructor for Data

  """ modified constructor for `Data` """
  function Data{T}(
    folders::Vector{<:Pair{String,<:Any}},
    savesecondsattype::Bool=false;
    sattype::Symbol=:undef,
    altmin::Real=5000,
    odelim::Union{Nothing,Char,String}=nothing,
    savedir::Union{String,Bool}="abs",
    maxtimediff::Int=30,
    primspan::Int=0,
    secspan::Int=15,
    lidarrange::Tuple{Real,Real}=(15_000,-Inf),
    stepwidth::Real=0.01,
    Xradius::Real=20_000,
    expdist::Real=Inf,
    remarks::Vector{<:Pair{String,<:Any}}=Pair{String,Any}[]
  ) where T

    # Load data
    tracks = MeasuredSet{T}(folders; sattype, altmin, odelim, savedir, remarks)

    # Process function arguments that need to be distributed to several structs
    folders = init_dict(folders, String[])
    remarks = init_dict(remarks, nothing)

    # Calculate Intersections
    intersections = (
    flight=Intersection{T}(
      tracks.flight, tracks.sat, savesecondsattype;
      maxtimediff, primspan, secspan, lidarrange,
      stepwidth, Xradius, expdist, savedir,
      remarks = remarks["Xflight"]
    ),
    cloud = Intersection{T}(
      tracks.cloud, tracks.sat, savesecondsattype;
      maxtimediff, primspan, secspan, lidarrange,
      stepwidth, Xradius, expdist, savedir,
      remarks = remarks["Xcloud"]
    ))

    # Instantiate
    new{T}(tracks, intersections)
  end # modified constructor for Data
end #struct Data

"""
    Data(args...; kwargs...)

Default `Data` constructor for single floating point precision.
"""
Data(args...; kwargs...) = Data{Float32}(args...; kwargs...)

"""
    Data{T}(data::Data) where T

Constructor for floating point conversions.
"""
Data{T}(data::Data) where T = Data{T}(
  FlightSet{T}(data.flight),
  CloudSet{T}(data.cloud),
  SatData{T}(data.sat),
  (flight = Intersection{T}(data.intersection.flight),
    cloud = Intersection{T}(data.intersection.cloud))
)

"""
    DataSet{T}(args...; kwargs...) where T

Alias constructor for `Data` of type `T`.
"""
DataSet{T}(args...; kwargs...) where T = Data{T}(args...; kwargs...)

"""
    DataSet(args...; kwargs...)

Alias default constructor for `Data` of type `Float32`.
"""
DataSet(args...; kwargs...) = DataSet{Float32}(args...; kwargs...)
