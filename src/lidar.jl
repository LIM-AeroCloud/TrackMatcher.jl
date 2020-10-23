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
  hfile = normpath(@__DIR__, "../data/CPro_Lidar_Altitudes_m.dat")
  hprofile = CSV.File(hfile, type=Float) |> df.DataFrame!
  # Consider only levels between max/min given in lidarrange
  itop = findfirst(hprofile.CPro .≤ lidarrange[1])
  ibottom = findlast(hprofile.CPro .≥ lidarrange[2])
  levels = try hprofile.CPro[itop:ibottom]
	catch
		AbstractFloat[]
	end

  # add extralayers for region with double precision
  h30m = findlast(levels.≥8300)
  hfine = try hfine = [(levels[i] + levels[i+1])/2 for i = h30m:length(levels)-1]
 		sort([levels; hfine; levels[end]-30], rev=true)
	catch
		AbstractFloat[]
	end

  # Return original and refined altitude profiles and important indices
  return (coarse = levels, fine = hfine, itop = itop === nothing ? 0 : itop,
		ibottom = ibottom === nothing ? 0 : ibottom, i30 = h30m === nothing ? 0 : h30m)
end #function get_lidarheights


"""
    function get_lidarcolumn(
      vect::Vector{<:Vector{<:Vector{<:Union{Missing,T}}}},
      ms::mat.MSession,
      variable::String,
      lidarprofile::NamedTuple,
      coarse::Bool = true;
      missingvalues = missing
    ) where T

Append the vector `vec` with CALIPSO lidar data of the `variable` using the MATLAB
session `ms` to read the variable from an hdf file already opened in `ms` outside
of `append_lidardata!`. Information about the lidar heigths used in `vec` are stored
in `lidarprofile`. Further information whether to use `coarse` levels (when set to
`true`, otherwise fine levels are used) is needed as input.
`missingvalues` can be set to any value, which will be replaced with `missing`
in `vec`.
"""
function get_lidarcolumn(
  T::DataType,
  ms::mat.MSession,
  variable::String,
  lidarprofile::NamedTuple;
  coarse::Bool = true,
  missingvalues = missing
)
  # Read variable from hdf file with MATLAB
	mat.eval_string(ms, "try\nvar = hdfread(file, '$variable');\nend")
	var = mat.jarray(mat.get_mvariable(ms, :var))
  # Initialise vector to store all row vectors and loop over matrix
  row = Vector{Vector{Union{Missing,T}}}(undef, size(var, 1))
	for i = 1:size(var, 1)
    row[i] = if coarse && ismissing(missingvalues)
      # Save row vector for coarse heights data without transforming missing values
      var[i,lidarprofile.itop:lidarprofile.ibottom]
    elseif coarse
      # Save row vector for coarse heights data after transforming missing values
      v = convert(Vector{Union{Missing,T}}, var[i,lidarprofile.itop:lidarprofile.ibottom])
      v[v.==missingvalues] .= missing
      v
    elseif ismissing(missingvalues)
      # Save row vector for coarse heights data after transforming missing values
  		v = Vector{T}(undef,length(lidarprofile.fine))
  		v[1:lidarprofile.i30-1] = var[i,lidarprofile.itop:lidarprofile.itop+lidarprofile.i30-2,1]
  		v[lidarprofile.i30:2:end] = var[i,lidarprofile.itop+lidarprofile.i30-1:lidarprofile.ibottom,1]
  		v[lidarprofile.i30+1:2:end] = var[i,lidarprofile.itop+lidarprofile.i30-1:lidarprofile.ibottom,2]
			v
    else
      # Save row vector for refined heights data after transforming missing values
      v = convert(Matrix{Union{Missing,T}}, var[i,lidarprofile.itop:lidarprofile.ibottom,:])
      v[v.==missingvalues] .= missing
  		v = Vector{T}(undef,length(lidarprofile.fine))
  		v[1:lidarprofile.i30-1] = var[i,lidarprofile.itop:lidarprofile.itop+lidarprofile.i30-2,1]
  		v[lidarprofile.i30:2:end] = var[i,lidarprofile.itop+lidarprofile.i30-1:lidarprofile.ibottom,1]
  		v[lidarprofile.i30+1:2:end] = var[i,lidarprofile.itop+lidarprofile.i30-1:lidarprofile.ibottom,2]
			v
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
classification(FCF::UInt16)::Tuple{UInt16,UInt16} = FCF & 7, FCF << -9 & 7


"""
    feature_classification(ftype::UInt16, fsub::UInt16) -> Symbol

From the feature classification `ftype` and `fsub`type, return a descriptive `Symbol`
explaining the feature classification type.
"""
function feature_classification(ftype::UInt16, fsub::UInt16)
  if ftype == 0
    missing
  elseif ftype == 1
    :clear
  elseif ftype == 2
    if fsub == 0
      :low_trasparent
    elseif fsub == 1
      :low_opaque
    elseif fsub == 2
      :transition_sc
    elseif fsub == 3
      :cu
    elseif fsub == 4
      :ac
    elseif fsub == 5
      :as
    elseif fsub == 6
      :ci
    elseif fsub == 7
      :cb
    end
  elseif ftype == 3
    if fsub == 0
      :aerosol
    elseif fsub == 1
      :marine
    elseif fsub == 2
      :dust
    elseif fsub == 3
      :polluted
    elseif fsub == 4
      :remote
    elseif fsub == 5
      :polluted_dust
    elseif fsub == 6
      :smoke
    elseif fsub == 7
      :other
    end
  elseif ftype == 4
    :stratospheric
  elseif ftype == 5
    :surface
  elseif ftype == 6
    :subsurface
  elseif ftype == 7
    :no_signal
  end
end


"""
		atmosphericinfo(
			sat::CPro,
			hlevels::Vector{<:AbstractFloat},
			isat::Int,
			flightalt::Real,
			flightid::Union{Int,String}
		)::Union{Missing,Symbol}

From the `CPro` cloud profile data at data point `isat` in `Intersection`
(index in the `DataFrame` of the intersection), return a `Symbol` with a human-readable
feature classification.

Use the `hlevels` in the lidar column data and the `flightalt`itude to determine
the atmospheric conditions (`feature`) at flight level at the intersection.

`atmosphericinfo` returns a `missing` value, if no height level overlap between the
flight altitude and the lidar levels was found or the feature array couldn't be accessed.
On errors, `flightid` will be returned to identify the source of the error.
"""
function atmosphericinfo(
  sat::CPro,
  hlevels::Vector{<:AbstractFloat},
	isat::Int,
  flightalt::Real,
	flightid::Union{Int,String}
)::Union{Missing,Symbol}
  i = argmin(abs.(hlevels .- flightalt))
  if abs(hlevels[i] - flightalt) > 60
    println(); @warn string("insufficient altitudes for lidar data saved; ",
      "missing used for feature in intersections of flight $(flightid)")
    missing
  else
    try sat.data.feature[isat][i]
    catch
      missing
    end
  end
end #function atmosphericinfo


"""
    atmosphericinfo(
      sat::CLay,
      alt::AbstractFloat,
      isat::Int
    )::Union{Missing,Symbol}

From the `CLay` cloud layer data at data point `isat` in `Intersection`
(index in the `DataFrame` of the intersection), return a `Symbol` with a human-readable
feature classification at flight `alt`itude.

`atmosphericinfo` returns a `missing` value, if no feature was found at flight level.
"""
function atmosphericinfo(
  sat::CLay,
  alt::AbstractFloat,
  isat::Int
)::Union{Missing,Symbol}
  top, base = sat.data.layer[isat]
  feature = try
    for i = 1:length(top)
      if base[i] ≤ alt ≤ top[i]
        return sat.data.feature[isat][i]
      end
    end
  catch
    return missing
  end
	isnothing(feature) && return missing
end #function atmosphericinfo
