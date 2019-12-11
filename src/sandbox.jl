################################################################################

### load FlightAware online data

altmin = 15_000
dir = "data/test/"
files = String[]; files = findFiles(files, dir, ".dat", ".txt")
file = files[4]
file = files[7]
onlineData = loadOnlineData(files)


################################################################################

### load FlightAware archive

altmin = 15_000
dir = "data/flightaware/archive/"
files = String[]; files = findFiles(files, dir, ".csv")


archive = loadArchive(files, altmin=altmin)



flights = loadFlightDB("a", "data/flightaware/archive/")

################################################################################


flights = loadFlightDB("ao", "data/flightaware/archive/", "data/flightaware/online/")

flights = CSV.read(file, datarow=2, normalizenames=true, dateformat="m/d/y H:M:S",
  types = Dict(:feet => Float64, :Groundspeed_knots_ => Float64))
names(flights)

flights = loadFlightDB("iao", "data/flightinventory/", "data/flightaware/archive/",
  "data/flightaware/online/")

sat = SatDB("data/CALIOP/")

f1 = flight = flights.inventory[1]
f2 = flights.inventory[2]
f3 = flights.inventory[3]
f4 = flights.inventory[4]


ms = mat.MSession()
mat.put_variable(ms, :x, flight.lon)
mat.put_variable(ms, :y, flight.lat)
mat.eval_string(ms, "p = pchip(x, y);");
pp = mat.get_mvariable(ms, :p)
fint=Minterpolate(ms, pp)
mat.eval_string("whos")

t1 = findfirst(sat.CLay.time .≥ flight.time[1] - Dates.Minute(30))
t2 = findlast(sat.CLay.time .≤ flight.time[end] + Dates.Minute(30))
c = (flight.metadata.area.latmin .≤ sat.CLay.lat[t1:t2] .≤ flight.metadata.area.latmax) .&
  ((flight.metadata.area.plonmin .≤ sat.CLay.lon[t1:t2] .≤ flight.metadata.area.plonmax) .|
  (flight.metadata.area.nlonmin .≤ sat.CLay.lon[t1:t2] .≤ flight.metadata.area.nlonmax))
cr = UnitRange[]
r = false
for i = 1:length(c)
  if c[i] && !r
    global r = true
    global ind = i
  elseif r && !c[i]
    r = false
    push!(cr, ind:i-1)
  end
end
cr


mat.put_variable(ms, :xs1, sat.CLay.lat[t1:t2][cr[1]])
mat.put_variable(ms, :ys1, sat.CLay.lon[t1:t2][cr[1]])
mat.eval_string(ms, "ps1 = pchip(xs1, ys1);");
pp1 = mat.get_mvariable(ms, :ps1)
s1int=Minterpolate(ms, pp1)

mat.put_variable(ms, :xs2, sat.CLay.lat[t1:t2][cr[2]])
mat.put_variable(ms, :ys2, sat.CLay.lon[t1:t2][cr[2]])
mat.eval_string(ms, "ps2 = pchip(xs2, ys2);");
pp2 = mat.get_mvariable(ms, :ps2)
s2int=Minterpolate(ms, pp2)

xi1 = sat.CLay.lat[t1:t2][cr[1][1]]:0.0001:sat.CLay.lat[t1:t2][cr[1][end]]
yi1 = s1int(xi1)
xi2 = sat.CLay.lat[t1:t2][cr[2][1]]:0.0001:sat.CLay.lat[t1:t2][cr[2][end]]
yi2 = s2int(xi2)

plt.plot(yi1,xi1,label="1")
plt.plot(yi2,xi2,label="2")
plt.gcf()
m1=argmin(abs.(xi1 .- fint(yi1)))
m2=argmin(abs.(xi2 .- fint(yi2)))

close(ms)


import PyPlot; plt = PyPlot


f1 = flights.inventory[1]
t1 = findfirst(sat.CLay.time .≥ f1.time[1] - Dates.Minute(30))
t2 = findlast(sat.CLay.time .≤ f1.time[end] + Dates.Minute(30))
c = (f1.metadata.area.latmin .≤ sat.CLay.lat[t1:t2] .≤ f1.metadata.area.latmax) .&
  ((f1.metadata.area.plonmin .≤ sat.CLay.lon[t1:t2] .≤ f1.metadata.area.plonmax) .|
  (f1.metadata.area.nlonmin .≤ sat.CLay.lon[t1:t2] .≤ f1.metadata.area.nlonmax))
plt.clf()
plt.scatter(sat.CLay.lon[t1:t2], sat.CLay.lat[t1:t2], marker=".", edgecolors=nothing,
  color="royalblue",label="CALIPSO track")
plt.scatter(sat.CLay.lon[t1:t2][c], sat.CLay.lat[t1:t2][c], marker=".", edgecolors=nothing,
  color="m",label="ε")
plt.scatter(f1.lon, f1.lat, marker=".", edgecolors=nothing, color="red",
  label="flight 1")
plt.scatter(yi1[m1], xi1[m1], label="intersection", marker = "s", color=:gold)
plt.scatter(yi2[m2], xi2[m2], marker = "s", color=:gold)
plt.xticks(-180:30:180)
plt.yticks(-90:30:90)
plt.xlabel("longitude / °")
plt.ylabel("latitude / °")
plt.grid(ls=":")
plt.minorticks_on()
ax = plt.axes()
ax.xaxis.set_minor_locator(plt.matplotlib.ticker.MultipleLocator(10))
ax.yaxis.set_minor_locator(plt.matplotlib.ticker.MultipleLocator(10))
plt.legend()
plt.xlim(-180,180); plt.ylim(-90, 90)
plt.gcf()
plt.savefig("data/FIG/f1.pdf")


f2 = flights.inventory[2]
t1 = findfirst(sat.CLay.time .≥ f2.time[1] - Dates.Minute(30))
t2 = findlast(sat.CLay.time .≤ f2.time[end] + Dates.Minute(30))
c = (f2.metadata.area.latmin .≤ sat.CLay.lat[t1:t2] .≤ f2.metadata.area.latmax) .&
  ((f2.metadata.area.plonmin .≤ sat.CLay.lon[t1:t2] .≤ f2.metadata.area.plonmax) .|
  (f2.metadata.area.nlonmin .≤ sat.CLay.lon[t1:t2] .≤ f2.metadata.area.nlonmax))
plt.clf()
plt.scatter(sat.CLay.lon[t1:t2], sat.CLay.lat[t1:t2], marker=".", edgecolors=nothing,
  color="royalblue",label="CALIPSO track")
plt.scatter(sat.CLay.lon[t1:t2][c], sat.CLay.lat[t1:t2][c], marker=".", edgecolors=nothing,
  color="m",label="ε")
plt.scatter(f2.lon, f2.lat, marker=".", edgecolors=nothing, color="red",
  label="flight 2")
plt.scatter(sint(xi[m]), xi[m], label="intersection", marker = "s", color=:gold)
plt.xticks(-180:30:180)
plt.yticks(-90:30:90)
plt.grid(ls=":")
plt.minorticks_on()
plt.legend()
plt.xlim(-180,180); plt.ylim(-90, 90)
plt.gcf()
plt.savefig("data/FIG/f2.pdf")


f3 = flights.inventory[3]
t1 = findfirst(sat.CLay.time .≥ f3.time[1] - Dates.Minute(30))
t2 = findlast(sat.CLay.time .≤ f3.time[end] + Dates.Minute(30))
c = (f3.metadata.area.latmin .≤ sat.CLay.lat[t1:t2] .≤ f3.metadata.area.latmax) .&
  ((f3.metadata.area.plonmin .≤ sat.CLay.lon[t1:t2] .≤ f3.metadata.area.plonmax) .|
  (f3.metadata.area.nlonmin .≤ sat.CLay.lon[t1:t2] .≤ f3.metadata.area.nlonmax))


f3.metadata.useLON
f3.metadata.flex

# Start MATLAB session
ms = mat.MSession()

# Fit satellite data
xs = sat.CLay.lat[t1:t2][c]
ys = sat.CLay.lon[t1:t2][c]
mat.put_variable(ms, :xs, xs)
mat.put_variable(ms, :ys, ys)
mat.eval_string(ms, "ps = pchip(xs, ys);");
ps = mat.get_mvariable(ms, :ps)
sint=Minterpolate(ms, ps)

f3int = NamedTuple{(:int,:min,:max),Tuple{Function,Float64,Float64}}[]
for seg in f3.metadata.flex
  mat.put_variable(ms, :x, f3.lat[seg.range])
  mat.put_variable(ms, :y, f3.lon[seg.range])
  mat.eval_string(ms, "ps = pchip(x, y);");
  ps = mat.get_mvariable(ms, :ps)
  push!(f3int, (int=Minterpolate(ms, ps),min=seg.min,max=seg.max))
end
f3int

xi = sat.CLay.lat[t1:t2][findfirst(c)]:0.001:sat.CLay.lat[t1:t2][findlast(c)]
yi = sint(xi)
mat.put_variable(ms, :i, mat.mxarray(xi)); mat.put_variable(ms, :p, ps)
mat.eval_string(ms, "pp = ppval(ps,i);")
mat.jvalue(mat.get_mvariable(ms, :pp))

plt.clf()
plt.scatter(yi,xi,label="int")
plt.legend()
plt.gcf()

z1 = abs.(abs.(sint(xi[f3int[1].min.≤ xi .≤f3int[1].max])) .- abs.(f3int[1].int(xi[f3int[1].min.≤xi.≤f3int[1].max])))
z2 = abs.(abs.(sint(xi[f3int[2].min.≤ xi .≤f3int[2].max])) .- abs.(f3int[2].int(xi[f3int[2].min.≤xi.≤f3int[2].max])))
z3 = abs.(abs.(sint(xi[f3int[3].min.≤ xi .≤f3int[3].max])) .- abs.(f3int[3].int(xi[f3int[3].min.≤xi.≤f3int[3].max])))
z4 = abs.(abs.(sint(xi[f3int[4].min.≤ xi .≤f3int[4].max])) .- abs.(f3int[4].int(xi[f3int[4].min.≤xi.≤f3int[4].max])))
z5 = abs.(abs.(sint(xi[f3int[5].min.≤ xi .≤f3int[5].max])) .- abs.(f3int[5].int(xi[f3int[5].min.≤xi.≤f3int[5].max])))
z = argmin(minimum.([isempty(z) ? Inf : z  for z in [z1,z2,z3,z4,z5]]))

m=argmin(abs.(sint(xi[f3int[z].min.≤ xi .≤f3int[z].max]) .- f3int[z].int(xi[f3int[z].min.≤ xi .≤f3int[z].max])))

close(ms)


plt.clf()
plt.scatter(sat.CLay.lon[t1:t2], sat.CLay.lat[t1:t2], marker=".", edgecolors=nothing,
  color="royalblue",label="CALIPSO track")
plt.scatter(sat.CLay.lon[t1:t2][c], sat.CLay.lat[t1:t2][c], marker=".", edgecolors=nothing,
  color="m",label="ε")
plt.scatter(f3.lon, f3.lat, marker=".", edgecolors=nothing, color="red",
  label="flight 3")
plt.xticks(-180:30:180)
plt.yticks(-90:30:90)
plt.grid(ls=":")
plt.minorticks_on()
ax = plt.axes()
ax.xaxis.set_minor_locator(plt.matplotlib.ticker.MultipleLocator(10))
ax.yaxis.set_minor_locator(plt.matplotlib.ticker.MultipleLocator(10))
plt.legend()
plt.xlim(-180,180); plt.ylim(-90, 90)
plt.gcf()
plt.savefig("data/FIG/f3.pdf")


f4 = flights.inventory[4]
t1 = findfirst(sat.CLay.time .≥ f4.time[1] - Dates.Minute(30))
t2 = findlast(sat.CLay.time .≤ f4.time[end] + Dates.Minute(30))
c = (f4.metadata.area.latmin .≤ sat.CLay.lat[t1:t2] .≤ f4.metadata.area.latmax) .&
  ((f4.metadata.area.plonmin .≤ sat.CLay.lon[t1:t2] .≤ f4.metadata.area.plonmax) .|
  (f4.metadata.area.nlonmin .≤ sat.CLay.lon[t1:t2] .≤ f4.metadata.area.nlonmax))
plt.clf()
plt.scatter(sat.CLay.lon[t1:t2], sat.CLay.lat[t1:t2], marker=".", edgecolors=nothing,
  color="royalblue",label="CALIPSO track")
plt.scatter(sat.CLay.lon[t1:t2][c], sat.CLay.lat[t1:t2][c], marker=".", edgecolors=nothing,
  color="m",label="ε")
plt.scatter(f4.lon, f4.lat, marker=".", edgecolors=nothing, color="red",
  label="flight 4")
plt.xticks(-180:30:180)
plt.yticks(-90:30:90)
plt.grid(ls=":")
plt.minorticks_on()
plt.legend()
plt.xlim(-180,180); plt.ylim(-90, 90)
plt.gcf()
plt.savefig("data/FIG/f4.pdf")


sat.CLay.time[1:20000]
flight = flights.inventory[1]

flights.inventory[1].time[1]
flights.inventory[end].time[end]

t0 = sat.CLay.time[1]
te = sat.CLay.time[end]

hits = FlightData[]
@pm.showprogress 1 for flight in flights.inventory
  t = (flight.time[1] - Dates.Minute(30) .≤ sat.CLay.time
    .≤ flight.time[end] + Dates.Minute(30))
  if !isempty(t)
    if any((flight.metadata.area.latmin .≤ sat.CLay.lat[t] .≤ flight.metadata.area.latmax) .&
      (flight.metadata.area.plonmin .≤ sat.CLay.lon[t] .≤ flight.metadata.area.plonmax) .|
      (flight.metadata.area.nlonmin .≤ sat.CLay.lon[t] .≤ flight.metadata.area.nlonmax))
      push!(hits, flight)
    end
  end
end


[f3.metadata.area[i] for i = 1:length(f3.metadata.area)]
f3.metadata.area
a = [(a=1, b=2), (a=3,b=4)]
Tuple(a)

flight.time
flight.metadata.flex[1].range

counter = 0
@pm.showprogress 1 for flight in flights.inventory
  if any((flight.metadata.area.latmin .≤ sat.CLay.lat[t] .≤ flight.metadata.area.latmax) .&
    (flight.metadata.area.plonmin .≤ sat.CLay.lon[t] .≤ flight.metadata.area.plonmax) .&
    (flight.metadata.area.nlonmin .≤ sat.CLay.lon[t] .≤ flight.metadata.area.nlonmax))
    push!(hits, flight)
  end
  global counter += 1
  if counter > 10000  break end
end

hits

t = collect(DateTime(2012,1,1,3):Dates.Minute(1):DateTime(2012,1,1,8))
t = ZonedDateTime.(t, tz.tz"UTC")

t = (ZonedDateTime(2012,1,1,3,0,0,tz.tz"UTC") .≤ sat.CLay.time .≤ ZonedDateTime(2012,1,1,4,0,0,tz.tz"UTC"))

late=extrema(sat.CLay.lat[t])
lpe=extrema(sat.CLay.lon[t][sat.CLay.lon[t].≥0])
lne=extrema(sat.CLay.lon[t][sat.CLay.lon[t].≤0])

ZonedDateTime(2012,1,1,3,0,0,tz.tz"UTC") .≤ s .≤ ZonedDateTime(2012,1,1,8,0,0,tz.tz"UTC")

flight = flights.inventory[1]
t = (flight.time[1] - Dates.Minute(30) .≤ sat.CLay.time
  .≤ flight.time[end] + Dates.Minute(30))



fi = findall(length.(timeintersects) .≤ 1)
ti = [timeintersects[i] for i in fi]
timeintersects = []
@pm.showprogress 1 for flight in flights.inventory
  t = findall(flight.time[1] - Dates.Minute(30) .≤ sat.CLay.time
    .≤ flight.time[end] + Dates.Minute(30))
  if !isnothing(t)
    findall
    push!(timeintersects, t)
  end
end
findall(isempty.(timeintersects))

@timed begin
  t = findall(te - Dates.Minute(30) .≤ sat.CLay.time .≤ te + Dates.Minute(30))
  tt = (t[1], t[end])
end

@timed (findfirst(te - Dates.Minute(30) .≤ sat.CLay.time .≤ te + Dates.Minute(30)),
  findlast(te - Dates.Minute(30) .≤ sat.CLay.time .≤ te + Dates.Minute(30)))

findfirst(tm.ZonedDateTime(2012,1,1,5,0,0,tm.tz.tz"UTC") .≤ t .≤ tm.ZonedDateTime(2012,1,1,6,0,0,tm.tz.tz"UTC"))

t[23280]

t0
t0 - tm.tz.Minute(30)
