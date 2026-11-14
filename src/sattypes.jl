## structs for satellite data

"""
# struct SatMetadata{T<:AbstractFloat} <: SatTrack{T}

Metadata for satellite data ([`SatData{T}`](@ref)), including date and area coverage, source file,
and root identifier.

see also [`SatData`](@ref), [`SatTrack`](@ref), [`SatSet`](@ref), and [`SecondaryMetadata`](@ref)

## Fields

- `date::NamedTuple{(:start, :stop),Tuple{DateTime,DateTime}}`: start and stop time of the satellite data coverage
- `area::NamedTuple{(:latmin, :latmax, :elonmin, :elonmax, :wlonmin, :wlonmax),NTuple{6,T}}`:
  latitude and longitude extrema (longitude is devided in eastern and western hemisphere for easier processing)
- `file::String`:` source file of the satellite data
- `root::UInt16`: root identifier for the source file
"""
struct SatMetadata{T<:AbstractFloat} <: SatTrack{T}
    date::NamedTuple{(:start, :stop),Tuple{DateTime,DateTime}}
    area::NamedTuple{(:latmin, :latmax, :elonmin, :elonmax, :wlonmin, :wlonmax),NTuple{6,T}}
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
    (start=DateTime(9999), stop=DateTime(9999)),
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

- `roots::OrderedDict{UInt16,String}`: mapping of root indices to root paths
- `resolution::Union{String,SortedSet{String}}`: the spatial resolution(s) of the satellite data
- `product::SatProduct`: the satellite product type
- `version::Union{String,SortedSet{String}}`: the version(s) of the satellite data
- `date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}`: the start and stop time of the monitored satellite period
- `created::AbstractDateTime`: the time of creation of the [`SatSet`](@ref)
- `loadtime::Dates.CompoundPeriod`: the time it took to read data files and load them into the struct
- `attachments`: any additional data or comments attached to the metadata

## Instantiation

`SecondaryMetadata` is constructed automatically, when `SatSet` is instantiated using
a modified constructor. Any `attachments` can be passed to the metadata with the keyword
argument `attachments` when constructing [`SatSet`](@ref).

Additional constructors help to instantiate empty `SecondaryMetadata` and to promote the
floating point type of an existing `SecondaryMetadata`. The floating point type can be omitted,
in which case `Float32` will be used by default.
"""
struct SecondaryMetadata{T<:AbstractFloat} <: SecondarySet{T}
    roots::OrderedDict{UInt16,String}
    # TODO mission::Enum{Caliop,EarthCare}
    resolution::Union{String,SortedSet{String}}
    product::SatProduct
    version::Union{String,SortedSet{String}}
    date::NamedTuple{(:start,:stop),Tuple{Date,Date}}
    created::AbstractDateTime
    loadtime::Dates.CompoundPeriod
    attachments

    #* Internal constructor with data checks
    function SecondaryMetadata{T}(
        roots::OrderedDict{UInt16,String},
        resolution::Union{String,SortedSet{String}},
        product::SatProduct,
        version::Union{String,SortedSet{String}},
        date::NamedTuple{(:start,:stop),Tuple{Date,Date}},
        created::AbstractDateTime,
        loadtime::Dates.CompoundPeriod,
        attachments
    ) where T<:AbstractFloat
        # Root ids are assigned sequentially from 1 when satellite files are loaded
        collect(keys(roots)) == UInt16.(1:length(roots)) ||
            throw(ArgumentError("SecondaryMetadata root ids must be sequential from 1\n"*
                "got: $(collect(keys(roots)))"))
        # ℹ Empty sets use a placeholder date in year 9999, so only the lower limit is checked
        Date(2000) ≤ date.start ≤ date.stop ||
            throw(ArgumentError("SecondaryMetadata date must be ordered and newer than 2000\n"*
                "got: $date"))
        new{T}(roots, resolution, product, version, date, created, loadtime, attachments)
    end
end #struct SecondaryMetadata

#* Constructor for empty SecondaryMetadata
SecondaryMetadata{T}() where T<:AbstractFloat = SecondaryMetadata{T}(
    OrderedDict{UInt16,String}(), "", unknown, "",
    (start=Dates.now(), stop=Dates.now()), Dates.now(), Dates.CompoundPeriod(), nothing
)

#* Constructor for default Float32 SecondaryMetadata
SecondaryMetadata(args...) = SecondaryMetadata{Float32}(args...)

#* Constructor for floating point type promotion
SecondaryMetadata{T}(meta::SecondaryMetadata) where T<:AbstractFloat = SecondaryMetadata{T}(
    copy_for_promotion(meta.roots), copy_for_promotion(meta.resolution), meta.product,
    copy_for_promotion(meta.version), meta.date, meta.created, meta.loadtime,
    copy_for_promotion(meta.attachments)
)


"""
# struct SatData{T<:AbstractFloat} <: SatTrack{T}

Satellite track data storing a single granule from an individual satellite data file with fields

- `time::Vector{DateTime}`: UTC time of satellite measurements
- `lat::Vector{T}`: latitudes of satellite measurements
- `lon::Vector{T}`: longitudes  of satellite measurements
- `metadata::SatMetadata{T}`: metadata associated with the granule

see also [`SatTrack`](@ref), [`SatMetadata`](@ref), and [`SatSet`](@ref)
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
)::SatData{T} where T<:AbstractFloat
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

See also: [`SatData`](@ref), [`SatSet`](@ref), and [`SatMetadata`](@ref)
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

    function SatSet{T}(
        folders::String...;
        type::SatProduct = unknown,
        attachments = nothing
    ) -> SatSet{T} where T<:AbstractFloat

Construct `SatSet{T}` from any number of absolute or relative folder paths given as string.
SatData searches for `.h5` files in all folders recursively and determines the data type
(`CLay` or `CPro`) from the file name. If multiple files are found, the type will be determined
from the majority of found files unless the `type` is specified directly.
Only one type of satellite data can be stored in `SatSet`. All floating point data will be
stored with the precision of `T`. If omitted, `Float32` will be used by default.
"""
struct SatSet{T<:AbstractFloat} <: SecondarySet{T}
    granules::StructArray{SatData{T}}
    metadata::SecondaryMetadata{T}

    #* Internal constructor with data checks
    function SatSet{T}(
        granules::StructArray{SatData{T}},
        metadata::SecondaryMetadata{T}
    ) where T<:AbstractFloat
        # Granules must be strictly ordered in time (required for granule selection)
        tstart = [m.date.start for m in granules.metadata]
        tstop = [m.date.stop for m in granules.metadata]
        all(start <= stop for (start, stop) in zip(tstart, tstop)) ||
            throw(ArgumentError("SatSet granule stop times must not precede start times"))
        all(>(Dates.Millisecond(0)), diff(tstart)) && all(>(Dates.Millisecond(0)), diff(tstop)) ||
            throw(ArgumentError("SatSet granules must be sorted by time without duplicates"))
        new{T}(granules, metadata)
    end
end


#* Main constructor for SatSet
function SatSet{T}(
    folders::String...;
    type::SatProduct = unknown,
    attachments = nothing
)::SatSet{T} where T<:AbstractFloat
    # Scan folders for satellite data files
    # Note: parse_metadata is only needed with #112, currently the file name can be used
    t0 = Dates.now()
    product, filedict = init_satimport(type, folders...)
    # Init loop over files
    total_files = filedict["filecount"]
    nfiles = 0
    granules = StructArray{SatData{T}}(undef, total_files)
    roots = OrderedDict{String,UInt16}()
    progress = pm.Progress(total_files; desc="Loading satellite data... ", enabled=progress_enabled())
    # Loop over found files and load data
    for root in keys(filedict["sources"]), file in filedict["sources"][root]
        try
            nfiles = load_satdata!(granules, file, root, roots, nfiles)
        catch e
            @warn "read error; data skipped" file exception=(e, catch_backtrace())
            pop!(granules)
        end
        pm.next!(progress)
    end
    pm.finish!(progress)
    # Create and return SatSet from loaded granules
    create_satset(granules, product, filedict, roots, t0, attachments)
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
