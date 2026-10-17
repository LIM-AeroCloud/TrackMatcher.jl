## structs for satellite data

struct SatMetadata{T<:AbstractFloat} <: SatTrack{T}
    date::NamedTuple{(:start, :stop),Tuple{DateTime,DateTime}}
    area::NamedTuple{(:latmin, :latmax, :elonmin, :elonmax, :wlonmin, :wlonmax),NTuple{6,T}}
    # REFACTOR source::SatProduct > move to secondary metadata
    file::String
    root::UInt16

    #* Internal constructor with data checks
    function SatMetadata{T}(
        date::NamedTuple{(:start, :stop),Tuple{DateTime,DateTime}},
        area::NamedTuple{(:latmin, :latmax, :elonmin, :elonmax, :wlonmin, :wlonmax),NTuple{6,T}},
        file::AbstractString,
        root::UInt16
    ) where T<:AbstractFloat
        date.start ≤ date.stop ||
            throw(ArgumentError("Unordered SatMetadata in $file: $date"))
        # NaN is valid for hemispheres without data, but only for both limits at once
        inrange(a, b, lo, hi) = isnan(a) ? isnan(b) : lo ≤ a ≤ b ≤ hi
        inrange(area.latmin, area.latmax, T(-90), T(90)) &&
            inrange(area.elonmin, area.elonmax, T(0), T(180)) &&
            inrange(area.wlonmin, area.wlonmax, T(-180), T(0)) ||
            throw(ArgumentError("SatMetadata area is out of range or inconsistent in $file: $area"))
        new{T}(date, area, file, root)
    end
end

#* Main constructor with limited tests and limited automated data construction
function SatMetadata{T}(
    date::Vector{DateTime},
    lat::Vector{T},
    lon::Vector{T},
    file::AbstractString,
    root::UInt16
) where T<:AbstractFloat
    elonmin, elonmax, wlonmin, wlonmax = lonextrema(lon)
    area = (latmin=minimum(lat), latmax=maximum(lat), elonmin, elonmax, wlonmin, wlonmax)
    SatMetadata{T}((start=date[begin], stop=date[end]), area, file, root)
end

#* Constructor for empty SatMetadata
SatMetadata{T}() where T<:AbstractFloat = SatMetadata{T}(
    (start=Dates.now(), stop=Dates.now()),
    (latmin=T(NaN), latmax=T(NaN), elonmin=T(NaN), elonmax=T(NaN), wlonmin=T(NaN), wlonmax=T(NaN)),
    "", 0x0000
)

#* Constructor for default Float32 SatMetadata
SatMetadata(args...) = SatMetadata{Float32}(args...)

#* Constructor for floating point type promotion
SatMetadata{T}(meta::SatMetadata) where T<:AbstractFloat =
    SatMetadata{T}(meta.date, map(T, meta.area), meta.file, meta.root)

"""
# struct SecondaryMetadata{T<:AbstractFloat} <: SecondarySet{T}

Immutable struct to hold metadata for `SatSet{T}` with fields
- `granules::DataFrame`
- `roots::OrderedDict{UInt16,String}`
- `type::Symbol`
- `date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}`
- `created::AbstractDateTime`
- `loadtime::Dates.CompoundPeriod`
- `attachments`


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
    roots::OrderedDict{UInt16,String}
    # TODO mission::Enum{Caliop,EarthCare}
    product::SatProduct
    version::Union{String,SortedSet{String}}
    date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
    created::AbstractDateTime
    loadtime::Dates.CompoundPeriod
    attachments

    #* Internal constructor with data checks
    function SecondaryMetadata{T}(
        roots::OrderedDict{UInt16,String},
        type::Symbol,
        date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}},
        created::AbstractDateTime,
        loadtime::Dates.CompoundPeriod,
        attachments
    ) where T<:AbstractFloat
    # Check that required columns are present in granules DataFrame
    colnames = ["file", "root", "tstart", "tstop", "latmin", "latmax", "elonmin", "elonmax", "wlonmin", "wlonmax"]
    coltypes = [Vector{String}, Vector{UInt16}, Vector{DateTime}, Vector{DateTime},
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
        throw(ArgumentError("tstart must be less than or equal to tstop for all granules\n"*
            "got: $(granules.tstart)\nvs.: $(granules.tstop)"))
    all(zip(granules.latmin, granules.latmax)) do (a, b)
        isnan(a) || isnan(b) || a ≤ b
    end || throw(ArgumentError("latmin must be less than or equal to latmax for all granules\n"*
        "got: $(granules.latmin)\nvs.: $(granules.latmax)"))
    all(zip(granules.elonmin, granules.elonmax)) do (a, b)
        isnan(a) || isnan(b) || a ≤ b
    end || throw(ArgumentError("elonmin must be less than or equal to elonmax for all granules\n"*
        "got: $(granules.elonmin)\nvs.: $(granules.elonmax)"))
    all(zip(granules.wlonmin, granules.wlonmax)) do (a, b)
        isnan(a) || isnan(b) || a ≤ b
    end || throw(ArgumentError("wlonmin must be less than or equal to wlonmax for all granules\n"*
        "got: $(granules.wlonmin)\nvs.: $(granules.wlonmax)"))
    new{T}(granules, roots, type, date, created, loadtime, attachments)
  end
end #struct SecondaryMetadata

#* Constructor for empty SecondaryMetadata
SecondaryMetadata{T}() where T = SecondaryMetadata{T}(
    DataFrame(file = String[], root = UInt16[], tstart = DateTime[], tstop = DateTime[],
        latmin = T[], latmax = T[], elonmin = T[], elonmax = T[], wlonmin = T[], wlonmax = T[]
    ),
    OrderedDict{UInt16,String}(),
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
    copy_for_promotion(meta.roots), meta.type, meta.date, meta.created, meta.loadtime, copy_for_promotion(meta.attachments)
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
    metadata::SatMetadata{T}

    #* Internal constructor with data checks
    function SatData{T}(
        time::Vector{DateTime},
        lat::Vector{T},
        lon::Vector{T},
        metadata::SatMetadata{T}
    ) where T<:AbstractFloat
        length(time) == length(lat) == length(lon) ||
            throw(DimensionMismatch("SatData vectors must have equal length"))
        if !isempty(time)
            issorted(time) || throw(ArgumentError("SatData time must be sorted"))
            checklimits(time, DateTime(2000), Dates.now(), "time")
            checklimits(lat, T(-90), T(90), "latitude")
            checklimits(lon, T(-180), T(180), "longitude")
            metadata.date == (start=time[begin], stop=time[end]) ||
                throw(ArgumentError("sat metadata does not seem to match sat data for granule $(metadata.file)"))
        end
        new{T}(time, lat, lon, metadata)
    end
end


#* Main constructor for SatData from file
function SatData{T}(
    file::AbstractString,
    root::AbstractString,
    roots::OrderedDict{String,UInt16}=OrderedDict{String,UInt16}()
) where T<:AbstractFloat
    # Open HDF5 file and read relevant datasets
    t, lat, lon = h5.h5open(joinpath(root, file), "r") do data
        read(data, "Profile_UTC_Time")[2,:],
        read(data, "Latitude")[2,:],
        read(data, "Longitude")[2,:]
    end
    # Convert data to TrackMatcher format
    time = convert_utc.(t)
    lat, lon = convert(Vector{T}, lat), convert(Vector{T}, lon)
    # Construct metadata
    rootid = get(roots, root, UInt16(length(roots) + 1))
    metadata = SatMetadata{T}(time, lat, lon, file, rootid)
    # Return SatData
    return SatData{T}(time, lat, lon, metadata)
end #constructor 2 SatData

#* Constructor for empty SatData
SatData{T}() where T<:AbstractFloat = SatData{T}(DateTime[], T[], T[], SatMetadata{T}())

#* Constructor for default Float32 SatData
SatData(args...) = SatData{Float32}(args...)

#* Constructor for floating point type promotion
SatData{T}(sat::SatData) where T<:AbstractFloat = SatData{T}(
  copy_for_promotion(sat.time), T.(sat.lat), T.(sat.lon), SatMetadata{T}(sat.metadata)
)

"""
    SatTrack{T}(args...) where T<:AbstractFloat

Alias constructor for [`SatData`](@ref). If `T` is omitted, `Float32` will be used by default.

See also: [`SatData`](@ref), [`SatSet`](@ref), and [`SecondaryMetadata`](@ref)
"""
function SatTrack end

SatTrack{T}(args...) where T<:AbstractFloat = SatData{T}(args...)
SatTrack(args...) = SatData{Float32}(args...)


"""
# struct SatSet{T<:AbstractFloat} <: SecondarySet{T}

Satellite data with fields

- `granules::StructArray{SatData{T}}`
- `metadata::SecondaryMetadata{T}`

The `granules` contains data from individual satellite data files as `StructArray{SatData{T}}`
with each `SatData{T}` storing continuous satellite track data from one file.

See also: [`SatTrack`](@ref), [`SatData`](@ref), and [`SecondaryMetadata`](@ref)

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
    # REFACTOR
    function SatSet{T}(
        granules::StructArray{SatData{T}},
        metadata::SecondaryMetadata{T}
    ) where T<:AbstractFloat
        # Granules must be ordered in time without duplicates (required for granule selection)
        tstart = [m.date.start for m in granules.metadata]
        tstop = [m.date.stop for m in granules.metadata]
        issorted(tstart; lt = ≤) && issorted(tstop) ||
            throw(ArgumentError("SatSet granules must be sorted by time without duplicates"))
        new{T}(granules, metadata)
    end
end


#* Main constructor for SatSet
function SatSet{T}(
    folders::String...;
    type::SatProduct = undef,
    attachments = nothing
)::SatSet{T} where T<:AbstractFloat
    # Scan folders for satellite data files
    t0 = Dates.now()
    paths = Vector{@NamedTuple{root::String,files::Vector{String}}}()
    for root in unique(realpath.(folders)) # ℹ skip duplicate folders
        files = scandir(root, ".h5")
        type = sat_datafiles!(files, type)
        push!(paths, (;root, files))
    end
    # Loop over found files and load data
    granules = StructArray{SatData{T}}(undef, 0)
    #// metadata = DataFrame(file = String[], root = UInt16[], tstart = DateTime[], tstop = DateTime[],
    #//     latmin = T[], latmax = T[], elonmin = T[], elonmax = T[], wlonmin = T[], wlonmax = T[])
    roots = OrderedDict{String,UInt16}()
    total_files = sum(length(data.files) for data in paths)
    progress = pm.Progress(total_files; desc="Loading satellite data... ", enabled=progress_enabled())
    for (root, files) in paths
        for file in files
            #// datafile = joinpath(root, file)
            try
                granule = SatData{T}(file, root, roots)
                #// elonmin, elonmax = lonextrema(granule.lon, ≥)
                #// wlonmin, wlonmax = lonextrema(granule.lon, <)
                #// push!(metadata, (;file, root = roots[root],
                #//     tstart = granule.time[1], tstop = granule.time[end],
                #//     latmin = minimum(granule.lat), latmax = maximum(granule.lat),
                #//     elonmin, elonmax, wlonmin, wlonmax
                #// ))
                push!(granules, granule)
                haskey(roots, root) || (roots[root] = granule.metadata.root)
            catch e
                @warn "read error; data skipped" file exception=(e, catch_backtrace())
            end
            pm.next!(progress)
        end
    end
    pm.finish!(progress)
    # Sort granules by time and drop duplicates (same start and stop time), e.g., from copies
    # of the same data in different folders
    by(m) = (m.date.start, m.date.stop, basename(m.file), m.root, m.file)
    issorted(granules.metadata; by) || (granules = granules[sortperm(granules.metadata; by)])
    meta = granules.metadata
    keep = (reverse∘unique)(i -> meta[i].date, Iterators.reverse(eachindex(meta))) # ℹ keep last (highest version)
    if length(keep) < length(meta)
        @warn "$(length(meta) - length(keep)) duplicate granules skipped"
        # Same time but different file names point to different product versions or renamed data files
        kept = Dict(meta[i].date => basename(meta[i].file) for i in keep)
        any(i -> basename(meta[i].file) ≠ kept[meta[i].date], setdiff(eachindex(meta), keep)) &&
            @warn "duplicate granules differ in file names (different product versions?); last one kept"
        granules = granules[keep]
    end
    # Return SatSet constructor
    # REFACTOR
    tend = Dates.now()
    date = if isempty(granules)
        @warn "no satellite data files successfully loaded"
        (start = DateTime(9999), stop = DateTime(9999))
    else
        tstart = [granule.metadata.date.start for granule in granules]
        tstop = [granule.metadata.date.stop for granule in granules]
        (start = minimum(tstart), stop = maximum(tstop))
    end
    # TODO: Secondary metadata and instantiation
    #=
    roots::OrderedDict{UInt16,String}
    # TODO mission::Enum{Caliop,EarthCare}
    product::SatProduct
    version::Union{String,SortedSet{String}}
    date::NamedTuple{(:start, :stop),Tuple{DateTime,DateTime}}
    created::AbstractDateTime
    loadtime::Dates.CompoundPeriod
    attachments
    =#
    return SatSet{T}(granules, SecondaryMetadata{T}(
        metadata,
        OrderedDict(v => k for (k,v) in roots),
        type,
        date,
        tz.ZonedDateTime(tend, tz.localzone()),
        Dates.canonicalize(tend - t0),
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


## Alias constructors for SatSet
SecondarySet{T}(args...; kwargs...) where T<:AbstractFloat = SatSet{T}(args...; kwargs...)
SecondarySet(args...; kwargs...) = SatSet{Float32}(args...; kwargs...)
