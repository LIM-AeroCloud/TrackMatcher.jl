## Define Metadata structs

"""
# struct CloudMetadata{T<:AbstractFloat} <: CloudTrack{T}

Metadata struct for cloud track data with fields `id`, `date` (date range `start` to `stop`),
`area`, `flex`, `use_lon`, `root` and `file`.

See also: [`FlightMetadata`](@ref), [`CloudData`](@ref), [`CloudTrack`](@ref), and [`CloudSet`](@ref)

## Instantiation

Use modified constructor to calculate `area` and `date` range from `time` and
`lat`/`lon` vectors in `data`.

    CloudMetadata{T}(
        data::DataFrame,
        id::Int32,
        flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}},
        use_lon::Bool,
        root::UInt16,
        file::UInt16
    ) where T<:AbstractFloat

Or hand over individual fields to the unmodified constructor.
"""
struct CloudMetadata{T<:AbstractFloat} <: CloudTrack{T}
    id::Int32
    date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
    area::NamedTuple{(:latmin,:latmax,:elonmin,:elonmax,:wlonmin,:wlonmax),NTuple{6,T}}
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}}
    use_lon::Bool
    root::UInt16
    file::UInt16
end

#* Main constructor for CloudMetadata with auto
function CloudMetadata{T}(
    data::DataFrame,
    id::Int32,
    flex::Tuple{Vararg{NamedTuple{(:range, :min, :max),Tuple{UnitRange{Int},T,T}}}},
    use_lon::Bool,
    root::UInt16,
    file::UInt16
) where T<:AbstractFloat
    # T = promote_type(eltype(data.lat), eltype(data.lon))
    elonmax = isempty(data.lon[data.lon.≥0]) ? T(NaN) : T(maximum(data.lon[data.lon.≥0]))
    elonmin = isempty(data.lon[data.lon.≥0]) ? T(NaN) : T(minimum(data.lon[data.lon.≥0]))
    wlonmax = isempty(data.lon[data.lon.<0]) ? T(NaN) : T(maximum(data.lon[data.lon.<0]))
    wlonmin = isempty(data.lon[data.lon.<0]) ? T(NaN) : T(minimum(data.lon[data.lon.<0]))
    area = (latmin=T(minimum(data.lat)), latmax=T(maximum(data.lat)),
    elonmin=elonmin, elonmax=elonmax, wlonmin=wlonmin, wlonmax=wlonmax)
    CloudMetadata{T}(id, (start=data.time[1], stop=data.time[end]), area, flex, use_lon, root, file)
end #constructor 2 CloudMetadata

#* Default constructor for type promotion to Float32
CloudMetadata(args...) = CloudMetadata{Float32}(args...)

#* Constructor for empty CloudMetadata
CloudMetadata{T}() where T<:AbstractFloat = CloudMetadata{T}(
    Int32(0), (start=Dates.now(), stop=Dates.now()), (latmin=T(NaN), latmax=T(NaN),
    elonmin=T(NaN), elonmax=T(NaN), wlonmin=T(NaN), wlonmax=T(NaN)),
    ((range=0:0, min=T(NaN), max=T(NaN)),), false, 0x0000, 0x0000
)

#* Constructor for type promotion from CloudMetadata with different float precision
CloudMetadata{T}(meta::CloudMetadata) where T<:AbstractFloat = CloudMetadata{T}(
    meta.id, meta.date, (latmin = T(meta.area.latmin), latmax = T(meta.area.latmax),
    elonmin = T(meta.area.elonmin), elonmax = T(meta.area.elonmax),
    wlonmin = T(meta.area.wlonmin), wlonmax = T(meta.area.wlonmax)),
    Tuple([(range = m.range, min = T.(m.min), max = T.(m.max)) for m in meta.flex]),
    meta.use_lon, meta.root, meta.file
)


## Define structs related to cloud data

"""
# struct CloudData{T<:AbstractFloat} <: PrimaryTrack{T}

Store Data related to a single cloud track.

See also: [`CloudMetadata`](@ref), [`CloudTrack`](@ref), [`CloudSet`](@ref), and [`FlightData`](@ref)

## Fields

- `time::Vector{DateTime}`
- `lat::Vector{T}`
- `lon::Vector{T}`
- `metadata::CloudMetadata{T}`

Default floating point precision is `Float32`. There are basic validity checks
during instantiation.

## Instantiation

The main constructor constructs `CloudData` from a `DataFrame` with columns `time`, `lat`, and
`lon`:

    CloudData{T}(data::DataFrame, metadata::CloudMetadata) where T
"""
struct CloudData{T<:AbstractFloat} <: PrimaryTrack{T}
    time::Vector{DateTime}
    lat::Vector{T}
    lon::Vector{T}
    metadata::CloudMetadata{T}

    #* Unmodified constructor for `CloudData` with basic checks for correct `data``
    function CloudData{T}(
        time::Vector{DateTime}, lat::Vector{T}, lon::Vector{T}, metadata::CloudMetadata{T}
        ) where T<:AbstractFloat
        # Column checks and warnings
        length(time) == length(lat) == length(lon) ||
            throw(DimensionMismatch("all data vectors must have the same length in CloudData"))
        checklimits(time, DateTime(2000), Dates.now(), "time")
        checklimits(lat, -90, 90, "latitude")
        checklimits(lon, -180, 180, "longitude")
        new{T}(time, lat, lon, metadata)
    end #constructor 1 CloudTrack
end #struct CloudTrack

#* Main constructor for `CloudData`
function CloudData{T}(
    data::DataFrame,
    metadata::CloudMetadata
)::CloudData{T} where T<:AbstractFloat
    CloudData{T}(data.time, data.lat, data.lon, CloudMetadata{T}(metadata))
end #constructor 2 CloudData

#* Default constructor for type promotion to Float32
CloudData(args...) = CloudData{Float32}(args...)

#* Constructor for empty CloudData
CloudData{T}() where T = CloudData{T}(DateTime[], T[], T[], CloudMetadata{T}())

#* Constructor for type promotion from CloudData with different float precision
CloudData{T}(cloud::CloudData) where T<:AbstractFloat = CloudData{T}(
    cloud.time, T.(cloud.lat), T.(cloud.lon), CloudMetadata{T}(cloud.metadata)
)

"""
    CloudTrack{T}(args...) where T
    CloudTrack(args...)

Alias constructor for `CloudData{T}`.
"""
function CloudTrack end

CloudTrack{T}(args...) where T<:AbstractFloat = CloudData{T}(args...)
CloudTrack(args...) = CloudData{Float32}(args...)


"""
# struct CloudSet{T<:AbstractFloat} <: PrimarySet{T}

Database for cloud track data with fields:

- `tracks::StructArray{CloudData{T}}`
- `metadata::PrimaryMetadata{T}`

# Instantiation

Instantiate by giving the folder paths to all `mat` files with cloud track data
and optional attachments.

    CloudSet{T}(folders::String...; attachments=nothing) where T

Or use the unmodified constructor.

    CloudSet{T}(tracks::StructArray{CloudData{T}}, metadata::PrimaryMetadata{T}) where T
"""
struct CloudSet{T<:AbstractFloat} <: PrimarySet{T}
    tracks::StructArray{CloudData{T}}
    metadata::PrimaryMetadata{T}
end

#* Main constructor for `CloudSet` with folder paths to mat files
function CloudSet{T}(
    folders::AbstractString...;
    structname::AbstractString="filtered_trajectories",
    attachments=nothing
) where T<:AbstractFloat
    # Return empty CloudSet, if no folders are passed to constructor
    isempty(folders) && return CloudSet{T}()
    # Track computing time
    tstart = Dates.now()
    # Scan folders for mat files
    pathdict = ds.OrderedDict{String,ds.OrderedDict}(
        "roots" => ds.OrderedDict{String,UInt16}(),
        "files" => ds.OrderedDict{String,UInt16}()
    )
    paths = findfiles!(pathdict, collect(folders), ".mat")

    #* Load cloud tracks from mat files into TrackMatcher in Julia format
    tracks = load_cloudtracks(paths, pathdict, structname, T)
    # Calculate load time
    tend = Dates.now()
    tc = tz.ZonedDateTime(tend, tz.localzone())
    loadtime = Dates.canonicalize(Dates.CompoundPeriod(tend - tstart))

    # For now find min/max times in CloudTracks
    tmin = minimum(t.time[1] for t in tracks)
    tmax = maximum(t.time[end] for t in tracks)
    # Define lookup dictionary for metadata
    roots = ds.OrderedDict{UInt16,String}(v => k for (k, v) in pathdict["roots"])
    files = ds.OrderedDict{UInt16,String}(v => k for (k, v) in pathdict["files"])
    pathdict = ds.OrderedDict{String,ds.OrderedDict}("roots" => roots, "files" => files)

    @info string("CloudSet loaded in ",
        "$(join(loadtime.periods[1:min(2,length(loadtime.periods))], ", ")) to",
        "\n▪ tracks ($(length(tracks)) entries)\n▪ metadata")

    # Instantiate CloudSet
    CloudSet{T}(tracks, PrimaryMetadata{T}(0, (start=tmin, stop=tmax), pathdict, tc, loadtime, attachments))
end #modified constructor 2

#* Constructor for empty CloudSet
CloudSet{T}() where T<:AbstractFloat = CloudSet{T}(CloudData{T}[], PrimaryMetadata{T}())

#* Default constructor for type promotion to Float32
CloudSet(args...; kwargs...) = CloudSet{Float32}(args...; kwargs...)

#* Constructor for type promotion from CloudSet with different float precision
CloudSet{T}(cloud::CloudSet) where T<:AbstractFloat = CloudSet{T}(
    CloudData{T}.(cloud.tracks),
    PrimaryMetadata{T}(cloud.metadata)
)

#* Alias constructor for `CloudSet` (documented in flighttypes.jl)
PrimarySet{T}(tracks::StructArray{CloudData{T}}, metadata::PrimaryMetadata{T}) where T<:AbstractFloat =
    CloudSet{T}(tracks, metadata)

PrimarySet{T}(folders::String...; kwargs...) where T<:AbstractFloat =
    CloudSet{T}(folders...; kwargs...)
