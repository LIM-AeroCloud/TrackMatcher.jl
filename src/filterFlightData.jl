#=
a = 1
b = 2
c = 0

@info begin
  global c = a + b
end

c

@debug begin
  global c = a * b
end

c

@info "useless"
@debug "useless"
=#


files = joinpath.(pwd(), "data/flightinventory/", readdir("data/flightinventory/")[.!startswith.(readdir("data/flightinventory/"), ".")])
# files = [joinpath(pwd(), "data/flightinventory/1_1_2012_SEGMENT.csv")]
tflights = @timed inventory = loadInventory(files)
flights = tflights[1]

sat = CLay(joinpath(pwd(), "data/CALIOP/"))

import PyPlot; const plt = PyPlot
plt.scatter(tflights[1][2].lon, tflights[1][2].lat, marker="s", color=:black)
plt.grid(ls=":")
plt.gcf()


ctypes = [Int64, Int64, Float64, Float64, Float64, Int64, Int64, Int64, Int64,
  Int64, Int64, Int64, String, String, String, Float64, String, Float64, Float64,
  Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64,
  Float64, Float64]
flights = CSV.read(files[1], datarow=3, footerskip=2,
  ignoreemptylines=true, silencewarnings=true, threaded=true, dateformat="HH:MM:SS.sssm")
lat = flights.LATITUDE[flights.FLIGHT_ID.==2]
lon = flights.LONGITUDE[flights.FLIGHT_ID.==2]
speed = flights.SPEED[flights.FLIGHT_ID.==2]
alt = flights.ALTITUDE[flights.FLIGHT_ID.==2]
t = flights.time[flights.FLIGHT_ID.==2]

lat = lat[alt.≥15000]
lon = lon[alt.≥15000]
lp = any(lon .> 0) ? maximum(lon[lon.≥0]) - minimum(lon[lon.≥0]) : 0
ln = any(lon .< 0) ? maximum(lon[lon.<0]) - minimum(lon[lon.<0]) : 0
useLON = maximum(lat) - minimum(lat) ≤ (lp + ln) * cosd(stats.mean(lat)) ? true : false
useLON ? (x = lon; y = lat) : (x = lat; y = lon)
x
flex = findFlex(x)
x[20:30]
y[20:30]

t = flights.SEGMENT_TIME
Time(t, "HH:MM:SS.sssm")
for t in flights.SEGMENT_TIME
  println(round(Time(t)))
end
println.(names(flights))
flights.time = [ZonedDateTime(flights.SEGMENT_YEAR[i], flights.SEGMENT_MONTH[i],
  flights.SEGMENT_DAY[i], flights.SEGMENT_HOUR[i], flights.SEGMENT_MIN[i],
  flights.SEGMENT_SEC[i], tz.tz"UTC") for i = 1:length(flights.FLIGHT_ID)]

findfirst(flights.time.==ZonedDateTime(2012,1,1,23,6,57, tz.tz"UTC"))
findall(flights.FLIGHT_ID.==3)
a,b,c,d,e = remdup(flights.LONGITUDE[r], flights.LATITUDE[r],
  flights.ALTITUDE[r], flights.SPEED[r], flights.time[r])

r = 278:411

flights.time[r][flights.ALTITUDE[r].≥15000]
flights.time[1:161][flights.ALTITUDE[1:161].≥15000]

dbID = [inventory[i].metadata.dbID for i = 1:length(inventory)]

import PyPlot; const plt = PyPlot

i = 1447
i = findfirst(isequal(i), dbID)
flex = findFlex(inventory[i].lat)
flex = findFlex(inventory[1].lon)

i = 1
f = 43
plt.clf()
# plt.scatter(inventory[i].lat, inventory[i].lon, marker="s", color="black")
# plt.scatter(inventory[i].lat[f], inventory[i].lon[f], marker="D", color="red")
plt.scatter(lat, lon, marker="D", color="red")
plt.scatter(inventory[i].lon[f], inventory[i].lat[f], marker="D", color="red")
plt.minorticks_on()
plt.grid(ls=":")
# plt.xlabel("latitude / °")
# plt.ylabel("longitude / °")
plt.ylabel("longitude / °")
plt.xlabel("latitude / °")
plt.gcf()



i = 1
i = 108
i = 1348
i = 1447

ms = mat.MSession()

lon = inventory[i].lon; lat = inventory[i].lat

lp = any(lon .> 0) ? maximum(lon[lon.≥0]) - minimum(lon[lon.≥0]) : 0
ln = any(lon .< 0) ? maximum(lon[lon.<0]) - minimum(lon[lon.<0]) : 0
useLON = maximum(lat) - minimum(lat) ≤ (lp + ln) * cosd(stats.mean(lat)) ? true : false

useLON ? (x = lon; y = lat) : (x = lat; y = lon)

f = findFlex(x)
itpfunctions = Function[]; itprange = NamedTuple{(:min,:max),Tuple{Real,Real}}[]
for (i, r) in enumerate(f)
  p = "p$i"
  mat.put_variable(ms, :x, x[r])
  mat.put_variable(ms, :y, y[r])
  mat.eval_string(ms, "$p = pchip(x, y);");
  pp = mat.get_mvariable(ms, Symbol("$p"))
  push!(itpfunctions, Minterpolate(ms, pp))
  push!(itprange, (min=r[1], max=r[end]))
end

xr = minimum(x):0.01:maximum(x)
xr2 = x[end]:0.01:maximum(x)
plt.clf()
# plt.scatter(inventory[i].lon, inventory[i].lat, marker="s", color="black")
# plt.plot(xr, itpfunctions[1](xr), color="red")
plt.scatter(inventory[i].lat, inventory[i].lon, marker="s", color="black")
plt.plot(xr, itpfunctions[1](xr), color="red")
plt.plot(xr2, itpfunctions[2](xr2), color="orange")
plt.minorticks_on()
plt.grid(ls=":")
plt.xlabel("latitude / °")
plt.ylabel("longitude / °")
# plt.xlabel("longitude / °")
# plt.ylabel("latitude / °")
plt.gcf()


x1 = float.(x)
mat.put_variable(ms, :x, mat.mxarray(x1[r]))
mat.eval_string(ms, "disp(x)")

x[f[1]] == sort(unique(x[f[1]]))
findfirst(x[f[1]] .== -180)
r = 1:43
p = "p1"
mat.put_variable(ms, :x, x[r])
mat.put_variable(ms, :y, y[r])
mat.eval_string(ms, "$p = pchip(x, y);");
pp = mat.get_mvariable(ms, Symbol("$p"))
push!(itpfunctions, Minterpolate(ms, pp))
push!(itprange, (min=r[1], max=r[end]))


r[1]
r[2]
r[end-1]
r[end]

mat.close(ms)
