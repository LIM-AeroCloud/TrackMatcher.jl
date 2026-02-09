## Helper functions for filesystem scans

"""
    findfiles!(
        rootdict::ds.OrderedDict{String,UInt16},
        roots::Union{AbstractString,Vector{<:AbstractString}},
        ext::Union{AbstractString,Vector{<:AbstractString}}
    ) -> Vector{@NamedTuple{root::String,files::Vector{String}}}

Find all files with the given extension(s) `ext` in the given `roots` folder(s) and return a
vector of named tuples with the root and the found files. Assign a unique ID to each new root
in `rootdict`.
"""
function findfiles!(
    rootdict::ds.OrderedDict{String,UInt16},
    roots::Union{AbstractString,Vector{<:AbstractString}},
    ext::Union{AbstractString,Vector{<:AbstractString}}
)::Vector{@NamedTuple{root::String,files::Vector{String}}}
    paths = @NamedTuple{root::String,files::Vector{String}}[]
    roots isa Vector || (roots = [roots])
    for root in roots
        files = scandir(root, ext)
        haskey(rootdict, root) || (rootdict[root] = length(rootdict) + 1)
        push!(paths, (; root = root, files))
    end
    return paths
end


"""
    scandir(root::AbstractString, ext::Union{AbstractString,Vector{<:AbstractString}}) -> Vector{String}

Scan `root` recursively for any files with the given extension(s) `ext` and return a vector of
strings with file paths relative to `root`.
"""
function scandir(root::AbstractString, ext::Union{AbstractString,Vector{<:AbstractString}})::Vector{String}
    # Setup
    ext isa Vector || (ext = [ext])
    root = realpath(root) # ℹ ensures path exists
    datafiles = String[]
    # Scan directory recursively for files with given extension(s)
    for (path, _, files) in walkdir(root)
        rpath = relpath(path, root)
        rpath == "." && (rpath = "")
        append!(datafiles, filter(e -> splitext(e)[2] in ext, joinpath.(rpath, files)))
    end
    # Return found data files
    return datafiles
end


"""
    sat_datafiles!(files::Vector{String}, type::Symbol) -> Symbol

Determine the satellite data type from `type` or, if type is `:undef`, from the given files
names and clean `files` of other data types.
Returns the detected data type as a symbol (`:CPro` or `:CLay`).
"""
function sat_datafiles!(files::Vector{String}, type::Symbol)::Symbol
    # Determine satellite data type
    if type === :undef
        nclay = count(occursin.("CLay", files))
        ncpro = count(occursin.("CPro", files))
        if nclay == ncpro == 0
            @warn "no satellite data files found, TrackMatcher requires CALIOP data"
            empty!(files)
            return type
        end
        type = ncpro ≥ nclay ? :CPro : :CLay
        @info "satellite data type auto-detected as $type"
    end
    # Clean files from other data types
    i = findall(!contains(string(type)), files)
    if length(i) > 0
        @warn "the given folder contains $(length(i)) non-$type files which will be ignored"
        deleteat!(files, i)
    end
    if length(files) == 0
        @warn "no $type data files found, TrackMatcher requires CALIOP data"
    end
    return type
end
