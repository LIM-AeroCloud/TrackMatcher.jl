## ObservationSet

"""
# struct CLay

CALIOP cloud layer `data` stored in a `DataFrame` with columns:
- `time::Vector{DateTime}` (time index)
- `lat::Vector{AbstractFloat}` (latitude position of current time index)
- `lon::Vector{AbstractFloat}` (longitude position of current time index)
- `layer::Vector{NamedTuple{(:top,:base),Tuple{Vector{<:AbstractFloat},Vector{<:AbstractFloat}}}}`
  (layer top/base heights in meters)
- `atmos_state::Vector{Vector{Symbol}}` (symbols describing the atmospheric conditions at the intersection)
- `OD::Vector{<:Vector{<:AbstractFloat}}` (layer optical depth)
- `IWP::Vector{<:Vector{<:Union{Missing,<:AbstractFloat}}}` (layer ice water path)
- `Ttop::Vector{<:Vector{<:AbstractFloat}}` (layer top temperature)
- `Htropo::Vector{<:AbstractFloat}` (tropopause height at current time index)
- `night::BitVector` (flag for nights (`true`))
- `averaging::Vector{<:Vector{Int}}` (horizontal averaging in km)

# Instantiation

    function CLay{T}(
        files::Vector{String},
        timespan::NamedTuple{(:min,:max), Tuple{DateTime,DateTime}},
        lidarrange::Tuple{Real,Real}=(15_000,-Inf),
        altmin::Real=5000
    ) where T

Construct `CLay` from a list of file names (including directories) and a running
MATLAB session `ms` and save data, if layers are within the bounds
of `lidarrange` and above flight `altmin` threshold and time is within `timesapn`.
If `T<:AbstractFloat` is not set, `Float32` will be used as default precision.

Or construct `CLay` by directly handing over the `DataFrame` where the names, order,
and types of each columns are checked and attempted to correct:

    CLay{T}(data::DataFrame) where T -> struct CLay
"""
struct CLay{T} <: ObservationSet{T}
    time::Vector{DateTime}
    lat::Vector{T}
    lon::Vector{T}
    layer_top::Vector{Vector{T}}
    layer_base::Vector{Vector{T}}
    atmos_state::Vector{Vector{Enum{UInt16}}}
    OD::Vector{Vector{T}}
    IWP::Vector{<:Vector{<:Union{Missing,<:T}}}
    Ttop::Vector{Vector{T}}
    Htropo::Vector{T}
    night::BitVector
    averaging::Vector{Vector{Int}}

    """ Unmodified constructor for `CLay` """
    function CLay{T}(
        time::Vector{DateTime},
        lat::Vector{T},
        lon::Vector{T},
        layer_top::Vector{Vector{T}},
        layer_base::Vector{Vector{T}},
        atmos_state::Vector{Vector{Enum{UInt16}}},
        OD::Vector{Vector{T}},
        IWP::Vector{<:Vector{<:Union{Missing,T}}},
        Ttop::Vector{Vector{T}},
        Htropo::Vector{T},
        night::BitVector,
        averaging::Vector{Vector{Int}}
    ) where T
        all(length(time) == length(prop) for prop in
            (lat, lon, layer_top, layer_base, atmos_state, OD, IWP, Ttop, Htropo, night, averaging)) ||
            throw(DimensionMismatch("all input vectors must have the same length in CLay data"))
        checklimits(time, DateTime(2000), Dates.now(), "time")
        checklimits(lat, T(-90), T(90), "latitude")
        checklimits(lon, T(-180), T(180), "longitude")
        checklimits.(OD, 0, 10, "optical depth")
        checklimits.(IWP, 0, 10000, "ice water path")
        checklimits.(Ttop, -120, 60, "layer top temperature")
        checklimits.(Htropo, 4000, 22_000, "tropopause height")
        new{T}(time, lat, lon, layer_top, layer_base, atmos_state, OD, IWP, Ttop, Htropo, night, averaging)
    end #constructor 1 CLay
end #struct CLay

"""
Modified constructor of `CLay` reading data from hdf `files` using MATLAB session `ms`
in the `lidarrange` (top to bottom), if data is above `altmin`.
"""
function CLay{T}(
    files::Vector{String},
    timeindex::Vector{UnitRange{Int}},
    lidarrange::Tuple{Real,Real}=(15_000,-Inf),
    altmin::Real=5000,
    saveobs::Bool=true
) where T
    # Return default empty struct if files are empty
    (isempty(files) || !saveobs) && return CLay{T}()
    # Initialise arrays
    # essential data
    utc = Vector{Vector{DateTime}}(undef, length(files))
    lat = Vector{Vector{T}}(undef, length(files))
    lon = Vector{Vector{T}}(undef, length(files))
    # non-essential data #TODO catch faulty data and allow missing values
    layer_top = Vector{Vector{Vector{T}}}(undef, length(files))
    layer_base = Vector{Vector{Vector{T}}}(undef, length(files))
    atmos_state = Vector{Vector{Vector{Enum{UInt16}}}}(undef,length(files))
    OD = Vector{Vector{Vector{T}}}(undef,length(files))
    IWP = Vector{Vector{Vector{Union{Missing,T}}}}(undef,length(files))
    Ttop = Vector{Vector{Vector{T}}}(undef,length(files))
    h_tropo = Vector{Vector{T}}(undef, length(files))
    night = Vector{BitVector}(undef, length(files))
    averaging = Vector{Vector{Int}}(undef,length(files))
    # Initialise variables for HDF5 reading to make them known outside try-catch block
    ltop, lbase, FCF, FOD, IWPath, LTT = nothing, nothing, nothing, nothing, nothing, nothing
    # Loop over files
    for (i, file) in enumerate(files)
        # Read hdf file
        h5.h5open(file, "r") do fid
            try
                utc[i] = convert_utc.(read(fid, "Profile_UTC_Time")[2,:][timeindex[i]])
                lat[i] = read(fid, "Latitude")[2,:][timeindex[i]]
                lon[i] = read(fid, "Longitude")[2,:][timeindex[i]]
                ltop = 1000read(fid, "Layer_Top_Altitude")[timeindex[i],:]
                lbase = 1000read(fid, "Layer_Base_Altitude")[timeindex[i],:]
                FCF = read(fid, "Atmospheric_Volume_Description")[timeindex[i],:]
                FOD = read(fid, "Feature_Optical_Depth_532")[timeindex[i],:]
                IWPath = read(fid, "Ice_Water_Path")[timeindex[i],:]
                LTT = read(fid, "Layer_Top_Temperature")[timeindex[i],:]
                h_tropo[i] = 1000read(fid, "Tropopause_Height")[timeindex[i]]
                night[i] = Bool.(read(fid, "Day_Night_Flag")[timeindex[i],:])
                averaging[i] = read(fid, "Horizontal_Averaging")[timeindex[i],:]
            catch
                @error "ReadError: layer observations could not be read from file, skipping" file
                return CLay{T}()
            end
        end

        # Loop over data and convert to TrackMatcher format
        l_top = Vector{Vector{T}}(undef,length(utc[i]))
        l_base = Vector{Vector{T}}(undef,length(utc[i]))
        atm = Vector{Vector{Enum{UInt16}}}(undef,length(utc[i]))
        optdepth = Vector{Vector{T}}(undef,length(utc[i]))
        icewater = Vector{Vector{Union{Missing,T}}}(undef,length(utc[i]))
        toptemp = Vector{Vector{T}}(undef,length(utc[i]))
        for n = 1:length(utc[i])
            l = findall((lbase[n,:] .> 0) .& (ltop[n,:] .> 0) .& (lbase[n,:] .< lidarrange[1]) .&
            (ltop[n,:] .> lidarrange[2]) .& (ltop[n,:] .> altmin))
            l_top[n], l_base[n], atm[n], optdepth[n], toptemp[n], icewater[n] =
            if isempty(l)
                T[], T[], Enum{UInt16}[], T[], T[], T[]
            else
                l = findall((lbase[n,:] .> 0) .& (ltop[n,:] .> 0) .& (lbase[n,:] .< lidarrange[1]) .&
                (ltop[n,:] .> lidarrange[2]))
                [ltop[n, m] for m in l] , [lbase[n, m] for m in l],
                Enum{UInt16}[feature_classification(classification(FCF[n,m])...) for m in l],
                [FOD[n,m] for m in l],
                [LTT[n,m] for m in l],
                [IWPath[n,m] == -9999 ? missing : IWPath[n,m] for m in l]
            end
        end # loop over time steps in current file
        layer_top[i], layer_base[i], atmos_state[i], OD[i], IWP[i], Ttop[i] =
            l_top, l_base, atm, optdepth, icewater, toptemp
    end #loop over files

    # Construct and standardise data
    CLay{T}(vcat(utc...), vcat(lat...), vcat(lon...), vcat(layer_top...), vcat(layer_base...),
        vcat(atmos_state...), vcat(OD...), vcat(IWP...), vcat(Ttop...), vcat(h_tropo...),
        vcat(night...), vcat(averaging...))
end #constructor 2 CLay


"""
    CLay{T}() where T

External constructor for empty `CLay` struct.
"""
function CLay{T}() where T<:AbstractFloat
    CLay{T}(DateTime[], T[], T[], Vector{T}[], Vector{T}[], Vector{Enum{UInt16}}[],
        Vector{T}[], Vector{T}[], Vector{T}[], T[], BitVector(), Vector{Int}[])
end

"""
    CLay{T}(clay::CLay) where T

External `CLay` constructor for conversion of floating point precision.
"""
function CLay{T}(clay::CLay) where T<:AbstractFloat
    CLay{T}(clay.time, T.(clay.lat), T.(clay.lon), [T.(layer) for layer in clay.layer_top],
        [T.(layer) for layer in clay.layer_base], clay.atmos_state,
        [T.(OD) for OD in clay.OD],
        [df.passmissing(T).(iwp) for iwp in clay.IWP],
        [T.(Ttop) for Ttop in clay.Ttop], T.(clay.Htropo), clay.night, clay.averaging)
end

""" Default CLay constructor for Float32 """
CLay(args...; kwargs...) = CLay{Float32}(args...; kwargs...)


"""
# struct CPro

CALIOP cloud profile `data` stored in a `DataFrame` with columns:
- `time::Vector{DateTime}` (current time index)
- `lat::Vector{AbstractFloat}` (latitude coordinate for current time index)
- `lon::Vector{AbstractFloat}` (lonitude coordinate for current time index)
- `atmos_state::Vector{<:Vector{<:Union{Missing,Symbol}}}`
  (symbols describing the atmospheric conditions for every height level at current time index)
- `EC532::Vector{<:Vector{<:Union{Missing,AbstractFloat}}}`
  (extinction coefficient at 532nm at every height level in current time index)

# Instantiation

    function CPro{T}(
      ms::mat.MSession,
      files::Vector{String},
      timespan::NamedTuple{(:min,:max), Tuple{DateTime,DateTime}},
      lidarprofile::NamedTuple
    ) where T -> struct CPro

Construct `CPro` from a list of file names (including directories) and a running
MATLAB session `ms`. CPro data is only stored in the vicinity of intersections for
the designated `timespan`. Column data is stored height-resolved as defined by the
`lidarprofile`. If `T<:AbstractFloat` is not set, `Float32` will be used as
default precision.

Or construct `CPro` by directly handing over the `DataFrame` where the names, order,
and types of each columns are checked and attempted to correct:

    CPro{T}(data::DataFrame) where T -> struct CPro
"""
struct CPro{T} <: ObservationSet{T}
    time::Vector{DateTime}
    lat::Vector{T}
    lon::Vector{T}
    atmos_state::Vector{Vector{Enum{UInt16}}}
    EC532::Vector{<:Vector{<:Union{Missing,T}}}
    h_tropo::Vector{T}
    temp::Vector{<:Vector{<:Union{Missing,T}}}
    pressure::Vector{<:Vector{<:Union{Missing,T}}}
    rH::Vector{<:Vector{<:Union{Missing,T}}}
    IWC::Vector{<:Vector{<:Union{Missing,T}}}
    deltap::Vector{<:Vector{<:Union{Missing,T}}}
    CADscore::Vector{<:Vector{<:Union{Missing,Int8}}}
    night::BitVector

    """ unmodified constructor """
    function CPro{T}(
        time::Vector{DateTime},
        lat::Vector{T},
        lon::Vector{T},
        atmos_state::Vector{Vector{Enum{UInt16}}},
        EC532::Vector{<:Vector{<:Union{Missing,T}}},
        h_tropo::Vector{T},
        temp::Vector{<:Vector{<:Union{Missing,T}}},
        pressure::Vector{<:Vector{<:Union{Missing,T}}},
        rH::Vector{<:Vector{<:Union{Missing,T}}},
        IWC::Vector{<:Vector{<:Union{Missing,T}}},
        deltap::Vector{<:Vector{<:Union{Missing,T}}},
        CADscore::Vector{<:Vector{<:Union{Missing,Int8}}},
        night::BitVector
    ) where T
        all(length(time) == length(prop) for prop in
            (lat, lon, atmos_state, EC532, h_tropo, temp, pressure, rH, IWC, deltap, CADscore, night)) ||
            throw(DimensionMismatch("all input vectors must have the same length in CPro data"))
        checklimits(time, DateTime(2000), Dates.now(), "time")
        checklimits(lat, T(-90), T(90), "latitude")
        checklimits(lon, T(-180), T(180), "longitude")
        # checklimits.(EC532, 0, 500, "extinction coefficient")
        checklimits(h_tropo, 4000, 22_000, "tropopause height")
        checklimits.(temp, -120, 60, "temperature")
        checklimits.(pressure, 1, 1086, "pressure")
        # checklimits.(rH, 0, 1.5, "relative humidity")
        # checklimits.(IWC, 0, 0.54, "ice water content")
        # checklimits.(deltap, 0, 1, "depolarization ratio")
        checklimits.(CADscore, -127, 128, "CAD score")
        new{T}(time, lat, lon, atmos_state, EC532, h_tropo, temp, pressure, rH, IWC, deltap, CADscore, night)
    end #constructor 1 CPro
end #struct CPro

"""
Modified constructor of `CPro` reading data from hdf `files` for all given `sattime` indices
and `lidarprofile` data, if data is above `altmin`.
"""
function CPro{T}(
    files::Vector{String},
    timeindex::Vector{UnitRange{Int}},
    lidarprofile::NamedTuple,
    saveobs::Bool=true
) where T
    # Return default empty struct if files are empty
    (isempty(files) || saveobs === false) && return CPro{T}()
    # Initialise arrays
    # essential data
    utc = Vector{Vector{DateTime}}(undef, length(files))
    lat = Vector{Vector{T}}(undef, length(files))
    lon = Vector{Vector{T}}(undef, length(files))
    fcf = Vector{Vector{Vector{<:Union{Missing,UInt16}}}}(undef, length(files))
    # non-essential data
    ec532 = Vector{Vector{Vector{<:Union{Missing,T}}}}(undef, length(files))
    h_tropo = Vector{Vector{T}}(undef, length(files))
    temp = Vector{Vector{Vector{<:Union{Missing,T}}}}(undef, length(files))
    pres = Vector{Vector{Vector{<:Union{Missing,T}}}}(undef, length(files))
    rH = Vector{Vector{Vector{<:Union{Missing,T}}}}(undef, length(files))
    iwc = Vector{Vector{Vector{<:Union{Missing,T}}}}(undef, length(files))
    deltap = Vector{Vector{Vector{<:Union{Missing,T}}}}(undef, length(files))
    cad = Vector{Vector{Vector{<:Union{Missing,Int8}}}}(undef, length(files))
    night = Vector{BitVector}(undef, length(files))
    # Initialise variables for HDF5 reading to make them known outside try-catch block
    FCF, EC532, IWC, Tcol, pressure, relH, Δp, CAD =
        nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing
    # Loop over files with cloud profile data
    for (i, file) in enumerate(files)
        ## Retrieve cloud profile data; assumes faulty files are filtered by SatData
        h5.h5open(file, "r") do fid
            try
                utc[i] = convert_utc.(read(fid, "Profile_UTC_Time")[2,timeindex[i]])
                lat[i] = read(fid, "Latitude")[2,timeindex[i]]
                lon[i] = read(fid, "Longitude")[2,timeindex[i]]
                fcf[i] = get_lidarcolumn(UInt16, read(fid, "Atmospheric_Volume_Description"),
                    lidarprofile, coarse=false)[timeindex[i]]
                ec532[i] = get_lidarcolumn(T, read(fid, "Extinction_Coefficient_532"),
                lidarprofile, missingvalues = -9999)[timeindex[i]]
                temp[i] = get_lidarcolumn(T, read(fid, "Temperature"), lidarprofile,
                    missingvalues = -9999)[timeindex[i]]
                pres[i] = get_lidarcolumn(T, read(fid, "Pressure"), lidarprofile,
                    missingvalues = -9999)[timeindex[i]]
                rH[i] = get_lidarcolumn(T, read(fid, "Relative_Humidity"), lidarprofile,
                    missingvalues = -9999)[timeindex[i]]
                iwc[i] = get_lidarcolumn(T, read(fid, "Ice_Water_Content_Profile"),
                    lidarprofile, missingvalues = -9999)[timeindex[i]]
                deltap[i] = get_lidarcolumn(T, read(fid, "Particulate_Depolarization_Ratio_Profile_532"),
                    lidarprofile, missingvalues = -9999)[timeindex[i]]
                cad[i] = get_lidarcolumn(Int8, read(fid, "CAD_Score"), lidarprofile,
                    coarse=false, missingvalues = -127)[timeindex[i]]
                h_tropo[i] = 1000read(fid, "Tropopause_Height")[timeindex[i]]
                night[i] = Bool.(read(fid, "Day_Night_Flag")[timeindex[i]])
            catch
                @error "ReadError: profile observations could not be read from file, skipping" file
                return CPro{T}()
            end
        end
    end #loop over files

    # Rearrange vector data
    utc, fcf = vcat(utc...), vcat(fcf...)
    # Convert FCF to feature classification enums
    atmos_state = Vector{Vector{Enum{UInt16}}}(undef, length(fcf))
    for i = 1:length(fcf)
        vect = Vector{Enum{UInt16}}(undef, length(fcf[i]))
        for j = 1:length(fcf[i])
            vect[j] = ismissing(fcf[i][j]) ? invalid :
            feature_classification(classification(fcf[i][j])...)
        end
        atmos_state[i] = vect
    end
    # Instantiate new struct
    CPro{T}(utc, vcat(lat...), vcat(lon...), atmos_state, vcat(ec532...), vcat(h_tropo...),
        vcat(temp...), vcat(pres...), vcat(rH...), vcat(iwc...), vcat(deltap...),
        vcat(cad...), vcat(night...))
end #constructor 2 CPro


"""
    CPro{T}() where T

External constructor for empty `CPro` struct.
"""
CPro{T}() where T = CPro{T}(DateTime[], T[], T[], Vector{Enum{UInt16}}[], Vector{T}[], T[],
    Vector{T}[], Vector{T}[], Vector{T}[], Vector{T}[], Vector{T}[], Vector{Int8}[], BitVector())

"""
    CPro{T}(cpro::CPro) where T

External `CPro` constructor for conversion of floating point precision.
"""
function CPro{T}(cpro::CPro) where T
    CPro{T}(cpro.time, T.(cpro.lat), T.(cpro.lon), cpro.atmos_state,
        [df.passmissing(T).(EC) for EC in cpro.EC532],
        T.(cpro.h_tropo),
        [df.passmissing(T).(temp) for temp in cpro.temp],
        [df.passmissing(T).(pres) for pres in cpro.pressure],
        [df.passmissing(T).(rh) for rh in cpro.rH],
        [df.passmissing(T).(iwc) for iwc in cpro.IWC],
        [df.passmissing(T).(dp) for dp in cpro.deltap],
        [df.passmissing(Int8).(cad) for cad in cpro.CADscore],
        cpro.night)
end

"""
    CPro(args...; kwargs...)

Default CPro constructor for Float32.
"""
CPro(args...; kwargs...) = CPro{Float32}(args...; kwargs...)
