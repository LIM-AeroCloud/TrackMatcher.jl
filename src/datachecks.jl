### Helper functions for system scans, data checks and corrections

## File system scans and data processing

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
    pathdict::ds.OrderedDict{String,ds.OrderedDict},
    roots::Union{AbstractString,Vector{<:AbstractString}},
    ext::Union{AbstractString,Vector{<:AbstractString}}
)::Vector{@NamedTuple{root::String,files::Vector{String}}}
    paths = @NamedTuple{root::String,files::Vector{String}}[]
    roots isa Vector || (roots = [roots])
    for root in roots
        files = scandir(root, ext)
        haskey(pathdict["roots"], root) || (pathdict["roots"][root] = length(pathdict["roots"]) + 1)
        for file in files
            haskey(pathdict["files"], file) || (pathdict["files"][file] = length(pathdict["files"]) + 1)
        end
        push!(paths, (; root, files))
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


"""
    remdup!(data::DataFrame, use_lon::Bool)

Remove entries with duplicate x and y (`lat`/`lon` or `lon`/`lat`) values from
`data` or increase `x` by an infinitessimal number if `x` data is identical, but
`y` data is not.
"""
function remdup!(data::DataFrame, use_lon::Bool)
    # Define x and y data
    x, y = use_lon ? (:lon, :lat) : (:lat, :lon)
    # Initialise
    i = 1
    iEnd = size(data, 1)
    # Loop over entries in vector
    while i < iEnd
        j = i + 1 # index for next consecutive line
        while j ≤ iEnd && isapprox(data[i, x], data[j, x], atol = eps(data[j, x]))
            if isapprox(data[j-1, y], data[j, y], atol = eps(data[j, y]))
                # Delete datarow, if x and y are identical
                data[j, x] = data[j-1, x]
                deleteat!(data, j-1)
                # Decrease the counter for the end of the arrays
                iEnd -= 1
            else
                # If x values are equal but y values differ,
                # add infinitessimal small increment to x value
                data[j, x] = nextfloat(data[j-1, x])
                j += 1 # set to next index
            end
        end
        i += 1 # increase line counter
    end
end #function remdup!


"""
    findflex(x::AbstractArray{<:Union{V,T} where V where T})
      -> Vector{ NamedTuple{(:range, :min, :max), Tuple{UnitRange, T, T}}} where T<:AbstractFloat

Find inflection points in `x` and return a vector of named tuples with index ranges
between inflection points and corresponding extrema.
"""
function findflex(x::AbstractArray{<:Union{V,T} where V where T})
    # Save starting index
    flex = Int[1]
    # Loop over data points
    for i = 2:length(x)-1
        # Find and save inflection points
        if x[i-1] > x[i] < x[i+1] || x[i-1] < x[i] > x[i+1]
            push!(flex, i)
        end
    end
    T = eltype(x)
    # Save end index
    push!(flex, length(x))
    # Initialise a vector for ranges between flex points, and min/max values
    ranges = NamedTuple{(:range, :min, :max), Tuple{UnitRange{Int}, T, T}}[]
    # Loop over saved inflection points
    for i = 2:length(flex)
        # Evaluate ascending/descending order of data and save extrema and range,
        # if it range contains more than one data point
        xmin, xmax = x[flex[i-1]] < x[flex[i]] ? (x[flex[i-1]], x[flex[i]]) :
        (x[flex[i]], x[flex[i-1]])
        flex[i] - flex[i-1] > 1 && push!(ranges, (range=flex[i-1]:flex[i], min=xmin, max=xmax))
    end

    # Return all ranges between flex points together with the corresponding extrema
    return Tuple(ranges)
end #function findflex


## Data checks and corrections

"""
    checklimits(data::Vector, min, max, property::AbstractString) -> Bool

Return `true`, if all values in `data` are between `min` and `max` or throw an ArgumentError.
The `property` name is used for the error message.
"""
function checklimits(data::Vector, min, max, property::AbstractString)::Bool
    # Return true if all values are missing (valid case for optional data)
    all(ismissing, data) && return true
    success = all(min .≤ skipmissing(data) .≤ max)
    success || throw(ArgumentError("$property values are out of range\n"*
        "got: [$(minimum(skipmissing(data))) ... $(maximum(skipmissing(data)))], expected: [$min ... $max]"))
end


"""
    checkcols!(
      data::DataFrame,
      standardnames::Vector{String},
      standardtypes::Vector{<:Type},
      bounds::Tuple{Vararg{Pair{<:Union{Int,Symbol},<:Tuple}}},
      dataset::T where T<:Union{Nothing,<:AbstractString}=nothing,
      id::Union{Nothing,Int,AbstractString}=nothing;
      essentialcols::Vector{Int}=[1,2,3]
    ) -> data::DataFrame

Check the columns of the `data` `DataFrame` for the use of correct `standardnames`,
`standardtypes`, and `bounds`. Throw an error for deviations in `essentialcols`,
otherwise issue a warning stating the `dataset` and `id` and try to correct the `data`.
"""
function checkcols!(
    data::DataFrame,
    standardnames::Vector{String},
    standardtypes::Vector{<:Type},
    bounds::Tuple{Vararg{Pair{<:Union{Int,Symbol},<:Tuple}}},
    dataset::T where T<:Union{Nothing,<:AbstractString}=nothing,
    id::Union{Nothing,Int,AbstractString}=nothing;
    essentialcols::Vector{Int}=[1,2,3]
)
    # Warn of non-standardised data
    check = (names(data) == standardnames && all([typeof(d) <: t for (d, t) in zip(df.eachcol(data), standardtypes)])) ||
        @warn "Non-standard names and/or order used for data columns. Trying to correct..." dataset id
    # Bring column bounds into the right format
    colbounds = definebounds(bounds, standardnames)
    # Setup vector holding checked and correct data column indices
    correctcols = zeros(Int, length(standardnames))
    # Find columns by name
    findbyname!(correctcols, data, standardnames, standardtypes, colbounds, check)
    # Find unchecked columns
    opencols = setdiff(collect(1:min(length(names(data)), length(standardnames))),
        filter(x -> x .≠ 0, correctcols))
    # Find columns by position
    findbyposition!(correctcols, opencols, data, standardtypes, colbounds)
    # Find columns in data not associated with a correct column
    remainingcols = setdiff(collect(1:length(names(data))),
        filter(x -> x .≠ 0, correctcols))
    # Find unchecked columns
    opencols = findall(isequal(0), correctcols)
    # Find columns by type and bounds
    findbytype!(data, correctcols, opencols, remainingcols, standardtypes, colbounds)
    # Correct DataFrame to the right format
    correctDF!(data, correctcols, standardnames, essentialcols)
    return data
end #function checkcols!


"""
    definebounds(
        bounds::Tuple{Vararg{Pair{<:Union{Int,Symbol},<:Tuple}}},
        colnames::Vector{Symbol}
    ) -> bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}}

From a Tuple of Pairs with column names or indices pointing to tuples with min/max
bounds, return a vector `bounds` of tuples with the min/max bounds corresponding
to each of the `colnames`. If `bounds` are not defined, they are set to `-Inf`/`Inf`
in the output vector.
"""
function definebounds(
    bounds::Tuple{Vararg{Pair{<:Union{Int,Symbol},<:Tuple}}},
    colnames::Vector{String}
)
    colbounds = Dict(bounds)
    bounds = Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}}(undef,length(colnames))
    for (i, name) in enumerate(colnames)
        lower, upper = if haskey(colbounds, name)
            colbounds[name]
        elseif haskey(colbounds, i)
            colbounds[i]
        else
            (-Inf, Inf)
        end
        bounds[i] = (lower, upper)
    end

    return bounds
end #function definebounds!


"""
    checklimits!(
        correctcols::Vector{Int},
        bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}},
        data::DataFrame,
        col::T where T <: Union{Int,Symbol},
        pos::Int,
        val::Int
    )
Check the `bounds` of a `col`umn in `data` and add the index `val`, if correct,
to `correctcols`. The index is added in `correctcols` at the `pos`ition, the column
will have in the final corrected DataFrame.
"""
function checklimits!(
    correctcols::Vector{Int},
    bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}},
    data::DataFrame,
    col::T where T <: Union{Int,String},
    pos::Int,
    val::Int
)
    if typeof(bounds[pos]) .== Tuple{DateTime,DateTime} &&
        all(bounds[pos][1] .< skipmissing(data[!,col]) .< bounds[pos][2])
        correctcols[pos] = val
    elseif typeof(bounds[pos]) .== Tuple{DateTime,DateTime}
        return
    elseif all(isinf.(bounds[pos])) ||
        all(bounds[pos][1] .< skipmissing(data[!,col]) .< bounds[pos][2])
        correctcols[pos] = val
    end
end #function checklimits!


"""
    findbyname!(
        correctcols::Vector{Int},
        data::DataFrame,
        standardnames::Vector{String},
        standardtypes::Vector{<:Type},
        bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}},
        check::T where T<:Union{Nothing,Bool}
    )

Find columns by `standardnames` in `data` and add the current position as `Int`
at the correct position in `correctcols`, if the type corresponds to the
`standardtypes` and `bounds` are correct. Issue a warning, if the previous `check`
found a mismatch of the column order, but all `data` columns could be identified by
`standardnames`.
"""
function findbyname!(
    correctcols::Vector{Int},
    data::DataFrame,
    standardnames::Vector{String},
    standardtypes::Vector{<:Type},
    bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}},
    check::T where T<:Union{Nothing,Bool}
)
    for (i, name) in enumerate(standardnames)
        col = findfirst(isequal(name), names(data))
        col ≠ nothing && typeof(data[!,name]) <: standardtypes[i] &&
        checklimits!(correctcols, bounds, data, name, i, col)
    end
    isempty(findall(isequal(0), correctcols)) && isnothing(check) &&
        @warn "all columns corrected based on column names"
end #function findbyname!


"""
    findbyposition!(
        correctcols::Vector{Int},
        opencols::Vector{Int},
        data::DataFrame,
        standardtypes::Vector{<:Type},
        bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}}
    )

Find columns by the `standardtypes` and `bounds` at each column index in `data`
and add the current in `correctcols` for positive matches. Issue a warning,
if all data could be corrected by type (and previously name).
"""
function findbyposition!(
    correctcols::Vector{Int},
    opencols::Vector{Int},
    data::DataFrame,
    standardtypes::Vector{<:Type},
    bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}}
)
    isempty(findall(isequal(0), correctcols)) && return
    for pos in opencols
        try typeof(data[!,pos]) <: standardtypes[pos] &&
            checklimits!(correctcols, bounds, data, pos, pos, pos)
        catch
            continue
        end
    end
    isempty(findall(isequal(0), correctcols)) &&
        @warn "all columns corrected based on\n- column names\n- column positions"
end #function findbyposition!


"""
    findbytype!(
        data::DataFrame,
        correctcols::Vector{Int},
        opencols::Vector{Int},
        remainingcols::Vector{Int},
        standardtypes::Vector{<:Type},
        bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}}
    )

Find the missing columns `opencols` in `data` using the `remainingcols` in `data`
which have not yet been assigned and add the current position as `Int`
at the correct position in `correctcols`, if `standardtypes` and `bounds` are
correct in a column. Issue a warning, if all data could be corrected by type
(and previously name and position).
"""
function findbytype!(
    data::DataFrame,
    correctcols::Vector{Int},
    opencols::Vector{Int},
    remainingcols::Vector{Int},
    standardtypes::Vector{<:Type},
    bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}}
)
    isempty(findall(isequal(0), correctcols)) && return
    for i in opencols, j in remainingcols
        try typeof(data[!,j]) <: standardtypes[i] && correctcols[i] == 0 &&
            checklimits!(correctcols, bounds, data, j, i, j)
        catch
            continue
        end
    end

    isempty(findall(isequal(0), correctcols)) &&
        @warn "all columns corrected based on\n- column names\n- column positions\n- column types and bounds"
end #function findbytype!


"""
    correctDF!(
        data::DataFrame,
        correctcols::Vector{Int},
        standardnames::Vector{String},
        essentialcols::Vector{Int}
    )

Correct the column order and/or names in `data` based on the order of indices in
`correctcols` and the `standardnames`. Throw an error for missing columns in
`essentialcols` otherwise issue warnings for missing columns filled with `missing`
values or additional columns being deleted.
"""
function correctDF!(
    data::DataFrame,
    correctcols::Vector{Int},
    standardnames::Vector{String},
    essentialcols::Vector{Int}
)
    for (i, col) in enumerate(findall(isequal(0), correctcols))
        data[!,Symbol("missing$i")] = [missing for i = 1:length(data[!,1])]
        correctcols[col] = length(names(data))
        col in essentialcols ? throw(ArgumentError("column $(standardnames[col]) not found in data")) :
        @warn("column $(standardnames[col]) not found in data; filled with `missing`")
    end
    additionalcols = setdiff(collect(1:length(names(data))), correctcols)
    !isempty(additionalcols) &&
        @warn "additional column $(join(additionalcols, ", ", " and ")) deleted in data"
    df.select!(data, correctcols)
    df.rename!(data, standardnames)
end #function correctDF!


## Prepare trajectory for interpolation and analysis of intersections

"""
    preptrack(flight::DataFrame) -> flight, flex, use_lon

Use the `flight` data to prepare the track for interpolation. Find the predominant
flight direction and inflection (flex) points in the flight track.

Return a tidied flight data from which duplicate entries are removed together
with the `flex` points in the x data for interpolation and a boolean `use_lon`,
which is true for longitude values used x data in the track interpolation.
"""
function preptrack!(trajectory::DataFrame)
    # calculate area covered by flight
    lp = any(trajectory.lon .≥ 0) ? maximum(filter(l -> l ≥ 0, trajectory.lon)) -
        minimum(filter(l -> l ≥ 0, trajectory.lon)) : 0
    ln = any(trajectory.lon .< 0) ? maximum(filter(l -> l < 0, trajectory.lon)) -
        minimum(filter(l -> l < 0, trajectory.lon)) : 0
    # Determine main direction of trajectory (N<>S, E<>W) and use it as x values
    # for trajectory interpolation (info stored as bool use_lon)
    use_lon = maximum(trajectory.lat) - minimum(trajectory.lat) ≤ (lp + ln) *
        cosd(stats.mean(trajectory.lat)) ? true : false
    # Adjust for duplicate entries and
    remdup!(trajectory, use_lon)
    # find flex points to cut data in segments needed for the interpolation
    flex = use_lon ? findflex(trajectory.lon) : findflex(trajectory.lat)

    return flex, use_lon
end #function preptrack


"""
    closest_points(arr::Vector{T}) where T<:AbstractFloat -> Tuple{Int,Int}

Return the indices of the elements in `arr` with the minimum value and the adjacent
next larger value.

Vector `arr` consists of distances between coinciding coordinate pairs of different tracks.
"""
function closest_points(arr::Vector{T})::Tuple{Int,Int} where T<:Real
    m1 = argmin(abs.(arr))
    m2 = if m1 == 1
        2
    elseif m1 == length(arr)
        m1 - 1
    elseif arr[m1-1] < arr[m1+1]
        m1-1
    else
        m1+1
    end
    return m1, m2
end #function closest_points


"""
    lonextrema(lon::Vector{T}, rel::Function) where T

Find extrema in `lon`gitude values for eastern and western hemisphere by defining
`rel` as `>` or `<` 0 filtering for the respective longitude values.
"""
function lonextrema(lon::Vector{T}, rel::Function) where T
    # Return minimum or NaN for non-existing track data in current hemisphere
    try minimum(filter(rel(0), lon))
    catch
        T(NaN)
    end,
    # Return maximum or NaN for non-existing track data in current hemisphere
    try maximum(filter(rel(0), lon))
    catch
        T(NaN)
    end
end #function lonextrema



"""
    withinlimits(area::NamedTuple, track, atol::Real=0.1) -> Vector{Bool}

Compares the coordinate pairs in columns `lat` and `lon` of the `DataFrame` `track`
to the `area` bounds `latmin`, `latmax`, `elonmin`, `elonmax`, `wlonmin`, and `wlonmax`
of the given `area`. Coordinates are allowed to exceed the `area` by `atol` degrees.
"""
function withinlimits(area::NamedTuple, track, atol::Real=0.1)
    if isnan(area.elonmin) || isnan(area.elonmax)
        (area.latmin - atol .≤ track.lat .≤ area.latmax + atol) .&
        (area.wlonmin - atol .≤ track.lon .≤ area.wlonmax + atol)
    elseif isnan(area.wlonmin) || isnan(area.wlonmax)
        (area.latmin - atol .≤ track.lat .≤ area.latmax + atol) .&
        (area.elonmin - atol .≤ track.lon .≤ area.elonmax + atol)
    else
        (area.latmin - atol .≤ track.lat .≤ area.latmax + atol) .&
        ((area.elonmin - atol .≤ track.lon .≤ area.elonmax + atol) .|
        (area.wlonmin - atol .≤ track.lon .≤ area.wlonmax + atol))
    end
end #function withinlimits


"""
    withinlimits(area::NamedTuple, atol::Real=0.1) -> function (track) -> Vector{Bool}

Generates a function that takes a `DataFrame` `track` with columns `lat` and `lon`
and compares the coordinate pairs in columns `lat` and `lon` to the boundaries
`latmin`, `latmax`, `elonmin`, `elonmax`, `wlonmin`, and `wlonmax` of the given
`area`. Coordinates are allowed to exceed the `area` by `atol` degrees.
"""
withinlimits(area::NamedTuple, atol::Real=0.1) = track -> withinlimits(area, track, atol)
