# Functions related to processing CALIOP data

"""
    get_lidarheights(lidarrange::Tuple{Real,Real}, Float::DataType=Float32) -> NamedTuple

Return a `NamedTuple` with the following entries in the `lidarrange` (top, bottom):
- `coarse`: altitude levels of the CALIOP lidar data as defined in the metadata of the hdf files
- `fine`: altitude levels of CALIOP lidar data with 30m intervals below 8.3km
- `itop`: index in the original data array of the first selected top height
- `ibottom`: index in the original data array of the last selected bottom height
- `i30`: vector index of the first altitude level with 30m intervals

Height levels in fiels `coarse` and `fine` are saved in single precision unless
otherwise specified by `Float`.
"""
function get_lidarheights(lidarrange::Tuple{Real,Real}, Float::DataType=Float32)
    # Read CPro lidar altitude profile
    hfile = normpath(@__DIR__, "..", "data", "CPro_Lidar_Altitudes_m.dat")
    hprofile = CSV.read(hfile, DataFrame, copycols = false, types=Float)
    # Validate lidarprofile indices
    if lidarrange[1] < lidarrange[2]
        throw(ArgumentError("invalid lidar height range\naltitude bounds are inverted"))
    end
    if lidarrange[1] < hprofile.CPro[end] || lidarrange[2] > hprofile.CPro[1]
        throw(ArgumentError("invalid lidar height range\nselected values must overlap with " *
            "$(hprofile.CPro[1]) — $(hprofile.CPro[end])"))
    end
    # Find starting point at the top of the atmosphere
    # ℹ min/max in top/bottom search ensure find results and graceful error handling of swapped lidarrange values
    top = findlast(hprofile.CPro .≥ min(lidarrange[1], hprofile.CPro[1]))
    bottom = findfirst(≤(max(lidarrange[2], hprofile.CPro[end])), hprofile.CPro)
    coarse = hprofile.CPro[top:bottom]
    # bottom = length(coarse)
    top, bottom
    # add extralayers for region with double precision
    h30m = findlast(coarse.≥8300)
    h30m = isnothing(h30m) ? 1 : h30m
    fine = [(coarse[i] + coarse[i+1])/2 for i = h30m:length(coarse)-1]
    length(fine) > 1 && (fine = [fine; coarse[end] + (coarse[end] - coarse[end-1])/2])
    fine = [coarse[1:h30m-1]; vec(permutedims(hcat(coarse[h30m:end], fine)))]

    # Get all necessary indices for the fine-grained array splicing
    ftop = length(fine) ≥ 2 && fine[2] ≥ lidarrange[1] ? 2 : 1
    fbottom = length(fine) ≥ 2 && fine[end-2] ≤ lidarrange[2] ? 2 : 1
    fine = fine[ftop:end-fbottom]

    # Return original and refined altitude profiles and important indices
    f = (; top = ftop, bottom = fbottom, h30m)
    return (; coarse, fine, i = (; top, bottom, f))
end #function get_lidarheights


"""
    function get_lidarcolumn(
        T::DataType,
        var::Array{Tv},
        lidarprofile::NamedTuple;
        coarse::Bool = true,
        missingvalues = missing
    ) where Tv

Append the vector `vec` with CALIPSO lidar data of the `variable` using the MATLAB
session `ms` to read the variable from an hdf file already opened in `ms` outside
of `append_lidardata!`. Information about the lidar heigths used in `vec` are stored
in `lidarprofile`. Further information whether to use `coarse` levels (when set to
`true`, otherwise fine levels are used) is needed as input.
Ensure output is a vector with element type `T`. Missing values marked in the input
data with certain values/key words can be flagged by the `missingvalues` keyword
argument and are replaced by `missing` in `vec`.
"""
function get_lidarcolumn(
    T::DataType,
    var::Array,
    lidarprofile::NamedTuple;
    coarse::Bool = true,
    missingvalues = missing
)
    # Convert input variable to array of union type Missing and T
    var = convert(Array{Union{Missing,T}, ndims(var)}, var)
    # Determine number of profiles (last dimension for HDF5 data)
    nprofiles = size(var)[end]
    # Initialise vector to store all row vectors and loop over profiles
    row = Vector{Vector{Union{Missing,T}}}(undef, nprofiles)
    for i = 1:nprofiles
        row[i] = if ndims(var) == 2 && !coarse
            throw(ArgumentError("invalid input: 2D array provided for fine resolution heights; " *
                "set coarse to true to resolve"))
        elseif coarse
            # Save coarse data and transform missing values
            v = ndims(var) == 2 ? var[lidarprofile.i.top:lidarprofile.i.bottom, i] :
                var[1, lidarprofile.i.top:lidarprofile.i.bottom, i]
            ismissing(missingvalues) ? v : replace!(v, missingvalues => missing)
        else
            # 3D array: Save fine resolution heights with optional missing-value transformation
            v = var[:, lidarprofile.i.top:lidarprofile.i.bottom, i]
            !ismissing(missingvalues) && replace!(v, missingvalues => missing)
            v = [v[1,1:lidarprofile.i.f.h30m-1];
                vec(permutedims(hcat(v[1,lidarprofile.i.f.h30m:end], v[2,lidarprofile.i.f.h30m:end])))]
            v[lidarprofile.i.f.top:end-lidarprofile.i.f.bottom]
        end
    end
    return row
end #function get_lidarcolumn


"""
    classification(FCF::UInt16) -> UInt16, UInt16

From the CALIPSO feature classification flag (`FCF`), return the `type` and `subtype`.


# Extended help
## Detailed CALIPSO description

This function accepts an array a feature classification flag (FCF) and
returns the feature type and feature subtype. Feature and subtype are
defined in Table 44 of the CALIPSO Data Products Catalog (Version 3.2).
Note that the definition of feature subtype depends on the feature type!

INPUTS
 FCF, an array of CALIPSO feature classification flags.
     Type: UInt16
OUTPUTS
 ftype, an array the same size as FCF containing the feature type.
     Type: Int
 subtype, an array the same size as FCF containing the feature subtype.
     Type: Int
Modified from M. Vaughn's "extract5kmAerosolLayerDescriptors.m" code by
J. Tackett September 23, 2010. jason.l.tackett@nasa.gov

## Feature type
These are the values returned and their definitions.
 0 = invalid (bad or missing data)
 1 = "clear air"
 2 = cloud
 3 = aerosol
 4 = stratospheric feature; polar stratospheric cloud (PSC) or stratospheric aerosol
 5 = surface
 6 = subsurface
 7 = no signal (totally attenuated)

NOTE: The definition of feature subtype will depend on the feature type!
See the feature classification flags definition table in the CALIPSO Data
Products Catalog for more details.

## Subtype definitions for AEROSOLS
 0 = not determined
 1 = clean marine
 2 = dust
 3 = polluted continental
 4 = clean continental
 5 = polluted dust
 6 = smoke
 7 = other

## Subtype definitions for CLOUDS
 0 = low overcast, transparent
 1 = low overcast, opaque
 2 = transition stratocumulus
 3 = low, broken cumulus
 4 = altocumulus (transparent)
 5 = altostratus (opaque)
 6 = cirrus (transparent)
 7 = deep convective (opaque)
"""
classification(FCF::UInt16)::Tuple{UInt16,UInt16} = FCF & 0b0111, FCF >> 9 & 0b0111


"""
    feature_classification(ftype::UInt16, fsub::UInt16) -> Enum{UInt16}

From the feature classification `ftype` and `fsub`type, return a descriptive `Enum` with the
human-readable feature classification.
"""
function feature_classification(ftype::UInt16, fsub::UInt16)
    feature = SkyCondition(ftype)
    feature = if feature == aerosol
        AerosolType(fsub)
    elseif feature == cloud
        CloudType(fsub)
    else
        feature
    end
    return feature
end


"""
    atmosphericinfo(
        sat::CPro,
        hlevels::Vector{<:AbstractFloat},
        isat::Int,
        flightalt::Union{Missing,Real},
        flight_num::Union{Int,String}
    ) -> Enum{UInt16}

From the `CPro` cloud profile data at data point `isat` in `Intersection`
(index in the `DataFrame` of the intersection), return a `Enum{UInt16}` with a human-readable
feature classification.

Use the `hlevels` in the lidar column data and the `flightalt`itude to determine
the atmospheric conditions (`feature`) at flight level at the intersection.

`atmosphericinfo` returns `invalid` if `flightalt` is missing, no suitable lidar level
is found near flight altitude, or the feature array cannot be accessed.

`flight_num` is used in warning/error messages to identify the source of retrieval issues.
"""
function atmosphericinfo(
    sat::CPro,
    hlevels::Vector{<:AbstractFloat},
    isat::Int,
    flightalt::Union{Missing,Real},
    flight_num::Union{Int,String}
)::Enum{UInt16}
    ismissing(flightalt) && return invalid
    isempty(hlevels) && return invalid
    i = argmin(abs.(hlevels .- flightalt))
    if abs(hlevels[i] - flightalt) > 60
        println(); @warn string("insufficient altitudes for lidar data saved; ",
        "invalid used for feature in intersections of flight $(flight_num)")
        invalid
    else
        state = get(sat.atmos_state, isat, Enum{UInt16}[])
        feature = get(state, i, invalid)
        if feature == invalid && (isat ∉ eachindex(sat.atmos_state) || i ∉ eachindex(state))
            @error "failed to retrieve atmospheric state for flight $(flight_num) "*
                "at altitude $(hlevels[i]); setting to 'invalid'"
        end
        feature
    end
end #function atmosphericinfo


"""
    atmosphericinfo(
        sat::CLay,
        alt::Union{Missing,Real},
        isat::Int
    ) -> Enum{UInt16}

From the `CLay` cloud layer data at data point `isat` in `Intersection`
(index in the `DataFrame` of the intersection), return a `Enum{UInt16}` with a human-readable
feature classification at flight `alt`itude.

`atmosphericinfo` returns `invalid` if `alt` is missing, no feature was found at flight
level, or layer data could not be accessed.
"""
function atmosphericinfo(
    sat::CLay,
    alt::Union{Missing,Real},
    isat::Int
)::Enum{UInt16}
    ismissing(alt) && return invalid
    top = get(sat.layer_top, isat, nothing)
    base = get(sat.layer_base, isat, nothing)
    state = get(sat.atmos_state, isat, nothing)
    if isnothing(top) || isnothing(base) || isnothing(state)
        @error "failed to retrieve atmospheric state at altitude $alt; setting to 'invalid'"
        return invalid
    end

    for i in eachindex(top)
        if base[i] ≤ alt ≤ top[i]
            return state[i]
        end
    end

    invalid
end #function atmosphericinfo
