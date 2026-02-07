## structs for satellite data

"""
# struct SecondaryMetadata{T<:AbstractFloat} <: SecondarySet{T}

Immutable struct to hold metadata for `SatSet{T}` with fields

- `files::Dict{Int,String}`
- `type::Symbol`
- `date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}`
- `created::Union{DateTime,ZonedDateTime}`
- `loadtime::Dates.CompoundPeriod`
- `remarks`

## fields

### granules

DataFrame with complimentary information on each data file, such as file `name` and `root`,
`tstart` and `tstop` time, and lat/lon extrema (overall and per hemisphere).

### roots

OrderedDict connecting root paths with indices for more efficient storage in `granules`.

### type

Symbol indicating, whether profile or layer data is stored.

### date

`NamedTuple` with fields `start` and `stop` for start and end time of the monitored
satellite period.

### created

time of creation of database

### loadtime

time it took to read data files and load it to the struct

### attachments

any additional data or comments that can be attached to the database


## Instantiation

`SecondaryMetadata` is constructed automatically, when `SatSet` is instantiated using
a modified constructor. Any `attachments` can be passed to the metadata with the keyword
argument `attachments` when constructing `SatSet`.

Additional constructors help to instantiate empty `SecondaryMetadata` and to promote the
floating point type of an existing `SecondaryMetadata`. The floating point type can be omitted,
in which case `Float32` will be used by default.
"""
struct SecondaryMetadata{T<:AbstractFloat} <: SecondarySet{T}
    granules::DataFrame
    roots::ds.OrderedDict{Int,String}
    type::Symbol
    date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
    created::AbstractDateTime
    loadtime::Dates.CompoundPeriod
    attachments

    #* Internal constructor with data checks
    function SecondaryMetadata{T}(
        granules::DataFrame,
        roots::ds.OrderedDict{Int,String},
        type::Symbol,
        date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
        created::AbstractDateTime,
        loadtime::Dates.CompoundPeriod,
        attachments
    ) where T<:AbstractFloat
    # Check that required columns are present in granules DataFrame
    colnames = ["file", "root", "tstart", "tstop", "latmin", "latmax", "elonmin", "elonmax", "wlonmin", "wlonmax"]
    coltypes = [Vector{String}, Vector{Int}, Vector{DateTime}, Vector{DateTime},
        Vector{T}, Vector{T}, Vector{T}, Vector{T}, Vector{T}, Vector{T}]
    essentialcols = collect(1:length(colnames))
    bounds = (:tstart => (Date(2000), Dates.now()), :tstop => (Date(2000), Dates.now()),
        :latmin => (-90, 90), :latmax => (-90, 90),
        :elonmin => (0, 180), :elonmax => (0, 180), :wlonmin => (-180, 0), :wlonmax => (-180, 0))
    checkcols!(granules, colnames, coltypes, bounds, "SecondaryMetadata"; essentialcols)
    granule_roots = (sort∘unique)(granules.root)
    root_keys = roots.keys
    granule_roots == root_keys ||
        throw(ArgumentError("values in root column of granules must match keys of the roots OrderedDict\n"*
            "got: $granule_roots\nvs.: $root_keys"))
    all(granules.tstart .≤ granules.tstop) ||
        throw(ArgumentError("tstart must be less than or equal to tstop for all granules"))
    all(granules.latmin .≤ granules.latmax) ||
        throw(ArgumentError("latmin must be less than or equal to latmax for all granules"))
    all(granules.elonmin .≤ granules.elonmax) ||
        throw(ArgumentError("elonmin must be less than or equal to elonmax for all granules"))
    all(granules.wlonmin .≤ granules.wlonmax) ||
        throw(ArgumentError("wlonmin must be less than or equal to wlonmax for all granules"))
    new{T}(granules, roots, type, date, created, loadtime, attachments)
  end
end #struct SecondaryMetadata

#* Constructor for empty SecondaryMetadata
SecondaryMetadata{T}() where T = SecondaryMetadata{T}(
    DataFrame(file = String[], root = Int[], tstart = DateTime[], tstop = DateTime[],
        latmin = T[], latmax = T[], elonmin = T[], elonmax = T[], wlonmin = T[], wlonmax = T[]
    ),
    ds.OrderedDict{Int,String}(),
    :undef, (start=Dates.now(), stop=Dates.now()), Dates.now(), Dates.CompoundPeriod(), nothing
)

#* Constructor for default Float32 SecondaryMetadata
SecondaryMetadata(args...) = SecondaryMetadata{Float32}(args...)

#* Constructor for floating point type promotion
SecondaryMetadata{T}(meta::SecondaryMetadata) where T = SecondaryMetadata{T}(
    DataFrame(file = meta.granules.file, root = meta.granules.root,
        tstart = meta.granules.tstart, tstop = meta.granules.tstop,
        latmin = T.(meta.granules.latmin), latmax = T.(meta.granules.latmax),
        elonmin = T.(meta.granules.elonmin), elonmax = T.(meta.granules.elonmax),
        wlonmin = T.(meta.granules.wlonmin), wlonmax = T.(meta.granules.wlonmax)
    ),
    meta.roots, meta.type, meta.date, meta.created, meta.loadtime, meta.attachments
)


"""
# struct SatData{T<:AbstractFloat} <: SatTrack{T}

Satellite track data storing a single granule from a individual satellite data file with fields

- `time::Vector{DateTime}`: UTC time of satellite measurements
- `lat::Vector{T}`: latitudes of satellite measurements
- `lon::Vector{T}`: longitudes  of satellite measurements

see also [`SatTrack`](@ref) and [`SatSet`](@ref)
"""
struct SatData{T<:AbstractFloat} <: SatTrack{T}
    time::Vector{DateTime}
    lat::Vector{T}
    lon::Vector{T}

    #* Internal constructor with data checks
    function SatData{T}(time::Vector{DateTime}, lat::Vector{<:AbstractFloat}, lon::Vector{<:AbstractFloat}) where T
        length(time) == length(lat) == length(lon) ||
            throw(DimensionMismatch("SatData vectors must have equal length"))
        checklimits(time, DateTime(2000), Dates.now())
        checklimits(lat, T(-90), T(90))
        checklimits(lon, T(-180), T(180))
        new{T}(time, lat, lon)
    end
end


#* Main constructor for SatData from file
function SatData{T}(file::AbstractString) where T<:AbstractFloat
    # Open HDF5 file and read relevant datasets
    t, lat, lon = h5.h5open(file, "r") do data
        read(data, "Profile_UTC_Time")[2,:],
        read(data, "Latitude")[2,:],
        read(data, "Longitude")[2,:]
    end
    # Convert time to DateTime UTC
    time = convert_utc.(t)
    # CReturn SatData
    return SatData{T}(time, lat, lon)
end #constructor 2 SatData

#* Constructor for empty SatData
SatData{T}() where T<:AbstractFloat = SatData{T}(DateTime[], T[], T[])

#* Constructor for default Float32 SatData
SatData(args...; kwargs...) = SatData{Float32}(args...; kwargs...)

#* Constructor for floating point type promotion
SatData{T}(sat::SatData) where T<:AbstractFloat = SatData{T}(
  sat.time, T.(sat.lat), T.(sat.lon)
)

"""
    SatTrack{T}(args...; kwargs...) where T<:AbstractFloat

Alias constructor for [`SatData`](@ref). If `T` is omitted, `Float32` will be used by default.
"""
function SatTrack end

SatTrack{T}(args...; kwargs...) where T<:AbstractFloat = SatData{T}(args...; kwargs...)
SatTrack(args...; kwargs...) = SatData{Float32}(args...; kwargs...)


"""
# struct SatSet{T<:AbstractFloat} <: SecondarySet{T}

Satellite data with fields
- `granules::StructArray{SatData{T}}`
- `metadata::SecondaryMetadata{T}`

The `granules` contains data from individual satellite data files as `StructArray{SatData{T}}`
with each `SatData{T}` storing continuous satellite track data from one file.

see also [`SatTrack`](@ref), [`SatData`](@ref), and [`SecondaryMetadata`](@ref)

## Instantiation

    function SaSet{T}(
      folders::String...;
      type::Symbol=:undef,
      attachments = nothing
    ) where T<:AbstractFloat -> struct SatData

Construct `SatSet{T}` from any number of absolute or relative folder paths given as string.
SatData searches for `.h5` files in all folders recursively and determines the data type
(`CLay` or `CPro`) from the majority of found files in the first valid folder automatically
unless the type is specified with a Symbol `:CLay` or `:CPro` by the keyword argument `type`.
Only one type of satellite data can be stored in `SatData`. All floating point data will be
stored with the precision of `T`. If omitted, `Float32` will be used by default.

Alternatively, handover all fields to the unmodified constructor, where basic data validity
checks will be performed:

    SatData{T}(granules::StructArray{SatData{T}}, metadata::SecondaryMetadata{T}) where T<:AbstractFloat

Floating point type promotion is possible with:

    SatSet{T}(sat::SatSet) where T<:AbstractFloat
"""
struct SatSet{T<:AbstractFloat} <: SecondarySet{T}
    granules::StructArray{SatData{T}}
    metadata::SecondaryMetadata{T}

    #* Internal constructor with data checks
    function SatSet{T}(
        granules::StructArray{SatData{T}},
        metadata::SecondaryMetadata{T}
    ) where T<:AbstractFloat
        # Check that granules and metadata lengths match
        length(granules) == df.nrow(metadata.granules) ||
            throw(DimensionMismatch("SatSet granules and metadata length mismatch"))
        new{T}(granules, metadata)
    end
end


#* Main constructor for SatSet
function SatSet{T}(
    folders::String...;
    type::Symbol = :undef,
    attachments = nothing
)::SatSet{T} where T<:AbstractFloat
    # Scan folders for satellite data files
    t0 = Dates.now()
    paths = []
    for folder in folders
        files = scandir(folder, ".h5")
        type = sat_datafiles!(files, type)
        push!(paths, (; root = folder, files))
    end
    # Loop over found files and load data
    granules = StructArray{SatData{T}}(undef, 0)
    metadata = DataFrame(file = String[], root = Int[], tstart = DateTime[], tstop = DateTime[],
        latmin = T[], latmax = T[], elonmin = T[], elonmax = T[], wlonmin = T[], wlonmax = T[]
    )
    roots = ds.OrderedDict{String,Int}()
    pm.@showprogress dt=1 desc="load sat data..." for data in paths
        root = realpath(data.root)
        haskey(roots, root) || (roots[root] = length(roots) + 1)
        for file in data.files
            datafile = joinpath(root, file)
            try
                granule = SatData{T}(datafile)
                elonmin, elonmax = lonextrema(granule.lon, ≥)
                wlonmin, wlonmax = lonextrema(granule.lon, <)
                push!(metadata, (;file, root = roots[root],
                    tstart = granule.time[1], tstop = granule.time[end],
                    latmin = minimum(granule.lat), latmax = maximum(granule.lat),
                    elonmin, elonmax, wlonmin, wlonmax
                ))
                push!(granules, granule)
            catch e
                @warn "read error; data skipped" file exception=(e, catch_backtrace())
            end
        end
    end
    # Return SatSet constructor
    tend = Dates.now()
    if isempty(metadata)
        @warn "No satellite data files successfully loaded"
        date = (start = Date(9999), stop = Date(9999))
    else
        date = (start = minimum(metadata.tstart), stop = maximum(metadata.tstop))
    end
    return SatSet{T}(granules, SecondaryMetadata{T}(
        metadata,
        ds.OrderedDict(v => k for (k,v) in roots),
        type,
        date,
        tz.ZonedDateTime(tend, tz.localzone()),
        tend - t0,
        attachments
    ))
end

#* Constructor for emtpy SatSet
SatSet{T}() where T<:AbstractFloat = SatSet{T}(
  StructArray{SatData{T}}(undef, 0),
  SecondaryMetadata{T}()
)

#* Constructor for default Float32 SatSet
SatSet(args...; kwargs...) = SatSet{Float32}(args...; kwargs...)

#* Constructor for floating point type promotion
SatSet{T}(sat::SatSet) where T<:AbstractFloat = SatSet{T}(
  SatData{T}.(sat.granules), SecondaryMetadata{T}(sat.metadata)
)


## stucts for observations
