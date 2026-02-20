### Helper functions for type, format and unit conversions

## Type conversions

"""
    convert_floats!(data::DataFrame, T::Type{<:AbstractFloat}=Float32) -> DataFrame

Transform all columns of Type `<:Union{Missing,AbstractFloat}` or `<:Vector{<:Union{Missing,AbstractFloat}}`
to the precision of type `T`. Returns the modified `data` DataFrame for piping.

Ensures output columns have type `Union{Missing, T}` only when missing values are present.
"""
function convert_floats!(data::DataFrame, T::Type{<:AbstractFloat}=Float32)::DataFrame
    for (i, col) in enumerate(eachcol(data))
        if eltype(col) <: Union{Missing, AbstractFloat}
            if any(ismissing, col)
                # Explicitly cast to Union{Missing, T} when missings are present
                data[!, i] = Vector{Union{Missing,T}}(df.passmissing(T).(col))
            else
                data[!, i] = Vector{T}(T.(col))
            end
        elseif eltype(col) <: Vector{<:Union{Missing, AbstractFloat}}
            # Handle nested vectors, using non-missing vectors when possible
            if isempty(col)
                data[!, i] = Vector{Vector{T}}()
            else
                if any(c -> any(ismissing, c), col)
                    data[!, i] = [Vector{Union{Missing, T}}(df.passmissing(T).(c)) for c in col]
                else
                    data[!, i] = [Vector{T}(T.(c)) for c in col]
                end
            end
        end
    end
    return data
end


## Time format conversion of satellite data

"""
    convert_utc(t::AbstractFloat) -> DateTime

Convert the CALIOP Profile UTC time (`t`) to a `DateTime`.
"""
function convert_utc(t::AbstractFloat)::DateTime
    # Extract date from Float before decimal point and convert to Date
    t, d = modf(t)
    d, t = Date(lpad(Int(d), 8, "0"), "yyyymmdd"), t * 86400
    d += Dates.Year(2000)
    # From overall time, calculate hours, minutes, and seconds
    h = floor(Int, t/3600)
    m = floor(Int, t - 3600h)÷60
    s = floor(Int, t - 3600h - 60m)
    ms = round(Int, 1000 * (t - 3600h - 60m - s), digits=3)

    # Return a DateTime from date and time (h/m/s) with timezone UTC
    return DateTime(Dates.yearmonthday(d)..., h, m, s, ms)
end #function convert_utc


## Data/unit conversion

"""
    earthradius(lat::T<:AbstractFloat) -> Float64

Calculate the Earth's radius `R` as function of the current `lat`itude
considering the ellipsoidal shape of the Earth due to the rotational flattening.
"""
function earthradius(lat::T)::Float64 where T<:AbstractFloat
    lat = Float64(lat) # ¡ ensure double precision for calculations
    req, rpol = 6378137., 6356752.  # ℹ equatorial and polar radius in meters
    √(((req^2*cosd(lat))^2 + (rpol^2*sind(lat))^2) / ((req*cosd(lat))^2 + (rpol*sind(lat))^2))
end #function earthradius


"""Convert feet to meters"""
function ft2m(ft::T)::T  where T<:Union{Missing,AbstractFloat}
    0.3048ft
end #function ft2m


"""Convert feet/min to m/s"""
function ftpmin2mps(ftpmin::T)::T  where T<:Union{Missing,AbstractFloat}
    0.00508ftpmin
end #function ftpmin2mps


"""Convert knots to m/s"""
function knot2mps(knot::T)::T  where T<:Union{Missing,AbstractFloat}
    1.852knot/3.6
end #function knot2mps
