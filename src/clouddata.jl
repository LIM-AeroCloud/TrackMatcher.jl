### Routines related to loading cloud track data

"""
    load_cloudtracks(files::String...; Float::DataType=Float32)

Load cloud track data from mat `files` and store in Julia format as `CloudTrack` structs.
Set the floating point precision to `Float` (default: `Float32`).
"""
function load_cloudtracks(
    paths::Vector{@NamedTuple{root::String,files::Vector{String}}},
    pathdict::ds.OrderedDict{String,ds.OrderedDict},
    structname::String="filtered_trajectories",
    Float::DataType=Float32
)
    # Initialise
    cloudtracks = StructArray{CloudData{Float}}(undef, 0)
    id = Int32(0)
    # Loop over all mat files
    ntracks = sum(length.(getfield.(paths, :files)))
    prog = pm.Progress(ntracks, desc="load cloud tracks...", enabled=progress_enabled())
    for path in paths, file in path.files
        # Read data from mat files
        tracks = matread(file, path.root, structname, Float)
        # Store data in Julia format as Vector of CloudTrack structs
        matsave!(cloudtracks, tracks, id, pathdict, file, path.root, Float)
        pm.next!(prog)
    end #loop over files
    pm.finish!(prog)

    return cloudtracks
end #function load_cloudtracks


"""
    matread(file::AbstractString, root::AbstractString, structname::String="cloud", Float::DataType=Float32)

Read cloud track data from a mat `file` in `root`. The saved struct has the top level `structname`.
For floats, the precision defined by `Float` is used (default: `Float32`).
"""
function matread(file::AbstractString, root::AbstractString, structname::String="cloud", Float::DataType=Float32)
    # Read data from mat file
    data = MAT.matread(joinpath(root, file))
    # Setup named tuple with time and lat/lon data
    records = length(data[structname]["timestamp"])
    tracks = (time = Vector{Vector{DateTime}}(), lat = Vector{Vector{Float}}(), lon = Vector{Vector{Float}}())
    sizehint!(tracks.time, records)
    sizehint!(tracks.lat, records)
    sizehint!(tracks.lon, records)

    # Loop over records and convert data to Julia format
    for i = 1:records
        t = vec(DateTime.(data[structname]["timestamp"][i], "yyyymmddHHMM"))
        lat, lon = collect.(Float, eachcol(data[structname]["centrLatLon"][i]))
        push!(tracks.time, t)
        push!(tracks.lat, lat)
        push!(tracks.lon, lon)
    end
    # Return data
    return tracks
end


"""
    matsave!(
        tracks::StructArray{<:CloudData},
        trackdata::NamedTuple{(:time,:lat,:lon)},
        id::Int32,
        pathdict::ds.OrderedDict,
        file::AbstractString,
        root::AbstractString,
        Float::DataType=Float32
    )

Append the vector with cloud `tracks` by timestamps `t` and coordinates `lat` and `lon` from
`trackdata` using the floating point precision set by `Float` (default: `Float32`). Pass on
`id`, `file`, and `root` to the metadata with the help of the `pathdict`.
"""
function matsave!(
    tracks::StructArray{<:CloudData},
    trackdata::NamedTuple{(:time,:lat,:lon)},
    id::Int32,
    pathdict::ds.OrderedDict,
    file::AbstractString,
    root::AbstractString,
    Float::DataType=Float32
)::Nothing
    # Transform MATLAB data into Julia Format and store as CloudTrack struct in a vector
    for i = 1:length(trackdata.time)
        # Determine predominant trajectory direction, inflection points, and remove duplicate entries
        data = DataFrame(
            time = trackdata.time[i],
            lat = trackdata.lat[i],
            lon = trackdata.lon[i]
        )
        flex, use_lon = preptrack!(data)
        isempty(flex) && continue
        id += Int32(1)
        push!(tracks, CloudData{Float}(data, CloudMetadata{Float}(
            data, id, flex, use_lon, pathdict["roots"][root], pathdict["files"][file]
        )))
    end
    return
end #function matsave!
