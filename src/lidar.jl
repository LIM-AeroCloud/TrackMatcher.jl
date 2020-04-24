# Functions related to processing CALIOP data

"""
    get_lidarheights() -> NamedTuple

Return a `NamedTuple` with the following entries:
- `coarse`: altitude levels of the CALIOP lidar data as defined in the metadata of the hdf files
- `fine`: altitude levels of CALIOP lidar data with 30m intervals below 8.3km
- `i30`: vector index of the first altitude level with 30m intervals
"""
function get_lidarheights(lidarrange::Tuple{Real,Real})
  # Read CPro lidar altitude profile
  hfile = normpath(@__DIR__, "../data/CPro_Lidar_Altitudes_km.dat")
  hprofile = CSV.read(hfile)
  # Consider only levels between max/min given in lidarrange
  itop = findfirst(hprofile.CPro .≤ lidarrange[1])
  ibottom = findlast(hprofile.CPro .≥ lidarrange[2])
  levels = hprofile.CPro[itop:ibottom]

  # add extralayers for region with double precision
  h30m = findlast(levels.≥8.3)
  hfine=[(levels[i] + levels[i+1])/2 for i = h30m:length(levels)-1]
  hfine = sort([levels; hfine; levels[end]-0.03], rev=true)

  # Return original and refined altitude profiles as Vector{AbstractFloat}
  return (coarse = levels, fine = hfine, itop=itop, ibottom=ibottom, i30=h30m)
end #function get_lidarheights


"""
    function append_lidardata!(
      vec::Vector{Vector{<:Union{Missing,T}}},
      ms::mat.MSession,
      variable::String,
      lidar::NamedTuple,
      coarse::Bool = false;
      missingvalues = missing
    ) where T


"""
function append_lidardata!(
  vec::Vector{Vector{<:Union{Missing,T}}},
  ms::mat.MSession,
  variable::String,
  lidar::NamedTuple,
  coarse::Bool = false;
  missingvalues = missing
) where T
	mat.eval_string(ms, "try\nvar = hdfread(file, '$variable');\nend")
	var = mat.jarray(mat.get_mvariable(ms, :var))
	for i = 1:size(var, 1)
    if coarse && ismissing(missingvalues)
      push!(vec, var[i,lidar.itop:lidar.ibottom])
    elseif coarse
      v = convert(Vector{Union{Missing,T}}, var[i,lidar.itop:lidar.ibottom])
      v[v.==missingvalues] .= missing
      push!(vec, v)
    else
  		v = Vector{T}(undef,length(lidar.fine))
  		v[1:lidar.i30-1] = var[i,lidar.itop:lidar.itop+lidar.i30-2,1]
  		v[lidar.i30:2:end] = var[i,lidar.itop+lidar.i30-1:lidar.ibottom,1]
  		v[lidar.i30+1:2:end] = var[i,lidar.itop+lidar.i30-1:lidar.ibottom,2]
  		push!(vec, v)
    end
	end

	return vec
end #function append_lidardata!


"""
classification(FCF::UInt16) -> UInt16, UInt16

From the CALIPSO feature classification flag (`FCF`), return the `type` and `subtype`.


# Detailed CALIPSO description

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
function feature_classification(ftype::UInt16, fsub::UInt16)::Symbol
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


function atmosphericinfo(
  sat::SatDB,
  sattype::Symbol,
  flight::FlightData,
  index::Tuple{Int,Int}
)::Union{Missing,Symbol}
  hflight = ft2km(flight.data.alt[index[1]])
  if sattype == :CPro
    hlevels = sat.CPro.metadata.lidarlevels.fine
    ifine = argmin(abs.(hlevels .- hflight))
    if hlevels[ifine] - hflight > 0.03
      println(); @warn string("insufficient altitudes for lidar data saved; ",
        "missing used for feature in intersections of flight $(flight.metadata.dbID)")
      missing
    else
      feature_classification(classification(sat.CPro.data.FCF[index[2]][ifine])...)
    end
  else #atmospheric volume information currently available only for CPro data
    missing
  end
end #function atmosphericinfo
