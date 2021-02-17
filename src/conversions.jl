### Helper functions for format and unit conversions

## Time format conversion of satellite data
"""
    convertUTC(t::AbstractFloat) -> DateTime

Convert the CALIOP Profile UTC time (`t`) to a `DateTime`.
"""
function convertUTC(t::AbstractFloat)
  # Extract date from Float before decimal point and convert to Date
  date = floor(Int, t)
  d = Date("20"*string(date), "yyyymmdd")
  # Calculate overall time in seconds
  utc = (t - date) * 86400
  # From overall time, calculate hours, minutes, and seconds
  h = floor(Int, utc/3600)
  m = floor(Int, utc - 3600h)÷60
  s = floor(Int, utc - 3600h - 60m)
  ms = round(Int, 1000(utc - 3600h - 60m - s))

  # Return a DateTime from date and time (h/m/s) with timezone UTC
  return DateTime(Dates.yearmonthday(d)..., h, m, s, ms)
end


## Data/unit conversion

"""
    earthradius(lat::T) -> R::T

Calculate the Earth's radius `R` in dependence of the current `lat`itude
condidering the ellipsoidal shape of the Earth due to the rotational flattening.
"""
function earthradius(lat::T)::T where T<:AbstractFloat
  req, rpol = 6378137, 6356752
  √(((req^2*cosd(lat))^2 + (rpol^2*sind(lat))^2) / ((req*cosd(lat))^2 + (rpol*sind(lat))^2))
end #function earthradius


"""Convert feet to kilometers"""
function ft2km(ft::T)::T  where T<:Union{Missing,AbstractFloat}
    0.0003048ft
end #function ft2km


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


## Type conversions

"""
    Float64(track::FlightData)

Overload `Float64` from Base to convert `FlightData{T}` to `FlightData{Float64}`.
"""
# function Float64(track::FlightData)
#   FlightData{Float64}(DataFrame(time=track.time, lat=Float64.(track.lat, lon=Float64.(track.lon),
#     alt=Float64.(track.alt), heading=track.heading), climb=Float64.(track.climb),
#     speed=Float64.(track.speed)),
#     FlightMetadata{Float64}([getproperty(track.metadata, f) for f in fieldnames(FlightMetadata)]...))
# end #Float64 for FlightData
