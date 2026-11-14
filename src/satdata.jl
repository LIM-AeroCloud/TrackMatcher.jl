## Helper functions to import satellite data

## Sat set parsing

"""
    load_satdata!(
        granules::StructArray{SatData{T<:AbstractFloat}},
        file::AbstractString,
        root::AbstractString,
        roots::OrderedDict{String,UInt16}=OrderedDict{String,UInt16}(),
        nfiles::Int=0
    ) -> Int

Load satellite data from an HDF5 `file` in a specified `root` directory and append it to the
`granules` array. Save new `roots` to the `roots` dictionary and increment and return the
file count (`nfiles`) accordingly.
"""
function load_satdata!(
    granules::StructArray{SatData{T}},
    file::AbstractString,
    root::AbstractString,
    roots::OrderedDict{String,UInt16}=OrderedDict{String,UInt16}(),
    nfiles::Int=0
)::Int where T<:AbstractFloat
    # Open HDF5 file and read relevant datasets
    t, lat, lon = h5.h5open(joinpath(root, file), "r") do data
        read(data, "Profile_UTC_Time")[2, :],
        read(data, "Latitude")[2, :],
        read(data, "Longitude")[2, :]
    end
    # Convert data to TrackMatcher format
    time = convert_utc.(t)
    lat, lon = convert(Vector{T}, lat), convert(Vector{T}, lon)
    # Construct metadata
    rootid = get(roots, root, UInt16(length(roots) + 1))
    metadata = SatMetadata{T}(time, lat, lon, file, rootid)
    haskey(roots, root) || (roots[root] = rootid)
    nfiles += 1
    try
        granules[nfiles] = SatData{T}(time, lat, lon, metadata)
    catch e
        @warn "read error; data skipped" file exception=(e, catch_backtrace())
        nfiles -= 1
        pop!(granules)
    end
    return nfiles
end #constructor 2 SatData


"""
    init_satimport(type::SatProduct, folders::String...) -> SatProduct, Dict

Initialize satellite import by scanning the specified `folders` for HDF5 files and
categorizing them based on the satellite product `type`. Returns a dictionary containing
file information for each product type.
"""
function init_satimport(type::SatProduct, folders::String...)::Tuple{SatProduct,Dict}
    # Init file dictionary for each product type
    valuedict = Dict(
        "version"=>SortedSet{String}(),
        "resolution"=>SortedSet{String}(),
        "sources" => OrderedDict{String, Vector{String}}(),
        "filecount" => 0
    )
    file_dict = Dict(
        unknown => deepcopy(valuedict),
        cpro => deepcopy(valuedict),
        clay => deepcopy(valuedict),
        apro => deepcopy(valuedict),
        alay => deepcopy(valuedict)
    )
    # Scan unique folders with removed duplicates for satellite data files
    for root in unique(realpath.(folders)), file in scandir(root, ".h5")
        add_file!(file_dict, file, root)
    end
    # Find and return product type with most files
    filter_files(file_dict, type)
end


"""
    add_file!(filedict::Dict, file::String, root::String)

Add the given `file` from the `root` directory to the `filedict` under the appropriate
satellite product type, determined from the file name.
"""
function add_file!(filedict::Dict, file::String, root::String)
    regex = r"^(?:\w+_)*?(\d+(?:km|m))?([A-Za-z]+)?(?:-[A-Za-z]+)*(?:[-_]?[Vv])?(\d+(?:[-.]\d+)?)?(?:\..*)?$"
    matches = match(regex, file)
    if isnothing(matches)
        haskey(filedict[unknown]["sources"], root) || (filedict[unknown]["sources"][root] = String[])
        push!(filedict[unknown]["sources"][root], file)
        filedict[unknown]["filecount"] += 1
        return
    end
    type = isnothing(matches.captures[2]) ? nothing : Symbol(lowercase(matches.captures[2]))
    product = if !isnothing(type) && type in Symbol.(instances(TrackMatcher.SatProduct))
        eval(type)
    else
        unknown
    end
    captures = map(x -> isnothing(x) ? nothing : String(x), matches.captures)
    !isnothing(captures[1]) && push!(filedict[product]["resolution"], captures[1])
    !isnothing(captures[3]) && push!(filedict[product]["version"], captures[3])

    haskey(filedict[product]["sources"], root) || (filedict[product]["sources"][root] = String[])
    push!(filedict[product]["sources"][root], file)
    filedict[product]["filecount"] += 1
    return
end


"""
    create_satset(
        granules::StructArray{SatData{T}},
        product::SatProduct,
        filedict::Dict,
        roots::OrderedDict{String,UInt16},
        t0::DateTime,
        attachments
    ) -> SatSet{T} where T<:AbstractFloat

Create a `SatSet` from the provided `granules` of satellite data and the respective sat set
metadata from the satellite `product`, the file information `filedict`, root directory mapping
`roots`, start time `t0`, and any additional `attachments`. Return the final `SatSet` object.
"""
function create_satset(
    granules::StructArray{SatData{T}},
    product::SatProduct,
    filedict::Dict,
    roots::OrderedDict{String,UInt16},
    t0::DateTime,
    attachments
)::SatSet{T} where T<:AbstractFloat
    # Ensure correct order of granules
    granules = sort(granules)
    # Return SatSet constructor
    tend = Dates.now()
    date = if isempty(granules)
        (start=Date(9999), stop=Date(9999))
    else
        tstart = [granule.metadata.date.start for granule in granules]
        tstop = [granule.metadata.date.stop for granule in granules]
        (start=minimum(tstart), stop=maximum(tstop))
    end
    # Return SatSet
    return SatSet{T}(granules, SecondaryMetadata{T}(
        OrderedDict(v => k for (k, v) in roots),
        unbox(filedict["resolution"]),
        product,
        unbox(filedict["version"]),
        date,
        tz.ZonedDateTime(tend, tz.localzone()),
        Dates.canonicalize(tend - t0),
        attachments
    ))
end


## Data processing and management helpers

"""
    Base.sort(granules::StructVector{<:SatData})

Sort the provided `granules` by their start and stop times, removing duplicates (same start and stop time).
Keeps the last occurrence of each unique time interval, which typically corresponds to the highest version.
"""
function Base.sort(granules::StructVector{<:SatData})
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
    return granules
end


"""
    filter_files(filedict::Dict, type::SatProduct) -> SatProduct, Dict

Filter the provided `filedict` for files of the given satellite `type`. If `type` is `unknown`,
the function determines the product type with the most files and returns it along with the
corresponding file information.
"""
function filter_files(filedict::Dict, type::SatProduct)::Tuple{SatProduct, Dict}
    # Filter for files of given product type
    if type !== unknown
        if filedict[type]["filecount"] > 0
            @info "loading $(filedict[type]["filecount"]) $type satellite files"
        else
            @warn "no valid $type satellite files found"
        end
        return type, filedict[type]
    end
    # Determine the product type with the most files if type is unknown
    prod, m = unknown, 0
    for product in instances(SatProduct)
        if product === unknown
            if filedict[unknown]["filecount"] > 0
                @warn "$(filedict[unknown]["filecount"]) files of unknown product type ignored"
            end
            continue
        end
        if filedict[product]["filecount"] > m
            if m > 0
                @warn "$m files for product $prod ignored"
            end
            m = filedict[product]["filecount"]
            prod = product
        elseif filedict[product]["filecount"] > 0
            @warn "$(filedict[product]["filecount"]) files for product $product ignored"
        end
    end
    # Check and fiter for the determined product type
    m > 0 ? @info("loading $m $prod files") : @warn("no valid satellite product files found")
    return prod, filedict[prod]
end


"""
    unbox(v::Union{AbstractVector,AbstractSet})

Return the single element contained in the collection `v`, if it has exactly one element,
otherwise return the collection itself.
"""
unbox(v::Union{AbstractVector,AbstractSet}) = length(v) == 1 ? only(v) : v



## Parsing satellite metadata

"""
    parse_metadata(f::HDF5.File, attr::String) -> OrderedDict

Parse the metadata from the `HDF5.File` object `f` stored in the attribute `attr` and
return it as an `OrderedDict`.
"""
function parse_metadata(f::h5.File, attr::String)::OrderedDict
    # Read metadata from the HDF5 attribute and store in array of Strings without whitespace
    meta = read_attribute(f, attr)
    lines = split(replace(meta, " " => ""), '\n')
    # Init loop over metadata
    i = 1
    metadata, keys = OrderedDict(), String[]
    key, value = "", ""
    while i ≤ length(lines)
        # Get current line, skip empty exit on END
        line = lines[i]
        isempty(line) && (i += 1; continue)
        line == "END" && break
        # Parse key, value pairs
        data = split(line, "=")
        key, value = length(data) == 2 ? data : ("", data)
        if key == "GROUP"
            # Add new section to dict for groups
            push!(keys, value)
            setnested!(metadata, keys, OrderedDict())
        elseif key == "END_GROUP"
            # Close dict section at end of group
            pop!(keys)
        elseif key == "OBJECT"
            # Add new objects as sub-dict entries
            save_object(metadata, value, keys, lines, i)
        end
        i += 1
    end
    # Remove outer dict layer and return dictionary data
    length(metadata) == 1 && (metadata = first(values(metadata)))
    return metadata
end


"""
    save_object(
        metadata::OrderedDict,
        obj_value::String,
        keys::Vector{String},
        lines::Vector{String},
        i::Int
    ) -> Int

Parse and save an object data (`obj_value`) from the metadata `lines` into the nested
`metadata` `OrderedDict`. Index `i` points to the current line in `lines` and `keys` are used
to navigate the nested `metadata` dict.
"""
function save_object(
    metadata::OrderedDict,
    obj_value::String,
    keys::Vector{String},
    lines::Vector{String},
    i::Int
)::Int
    # Init loop over object
    key, value = "", ""
    last_index = findfirst(isequal("END_OBJECT="*obj_value), lines)
    while key ≠ "END_OBJECT" && value ≠ obj_value
        # Set line counter and skip empty lines
        i += 1
        isempty(lines[i]) && continue
        # Process key/value pairs
        key, value = split(lines[i], "=")
        if key == "NUM_VAL"
            # Only parse scalar values
            len = tryparse(Int, value)
            len == 1 || return last_index
        elseif key == "CLASS"
            # Ignore class objects
            return last_index
        elseif key == "VALUE"
            # Parse value and covert to appropriate type
            value = if startswith(value, "\"")
                value = strip(value, '"')
                try
                    DateTime(value[1:(end-1)], Dates.dateformat"y-m-dTH:M:S")
                catch
                    value
                end
            else
                v = tryparse(Int, value)
                isnothing(v) ? parse(Float64, value) : v
            end
            # Save parsed object in metadata dict
            setnested!(metadata, [keys; obj_value], value)
            return last_index
        end
    end
end


"""
    getnested(dict, keys) -> Any

Retrieve the value from the nested `dict` using the hierarchy from the `keys` vector.
"""
getnested(dict, keys) = foldl(getindex, keys; init=dict)


"""
    setnested!(dict, keys, value) -> Nothing

Set the `value` in the nested `dict` using the hierarchy from the `keys` vector.
"""
function setnested!(dict, keys, value)
    d = getnested(dict, keys[1:(end-1)])
    d[keys[end]] = value
end
