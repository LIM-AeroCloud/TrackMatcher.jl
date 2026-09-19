# Globally used test data
flight = FlightSet(volpe=joinpath(@__DIR__, "data", "volpe"))
flight64 = FlightSet{Float64}(volpe=joinpath(@__DIR__, "data", "volpe"))
t0 = now()
flight_empty = FlightSet()
cloud = CloudSet(joinpath(@__DIR__, "data", "cloud"))

cpro_src = joinpath(@__DIR__, "data", "caliop", "CPro")
clay_src = joinpath(@__DIR__, "data", "caliop", "CLay")
sat_cpro = SatSet(cpro_src, type=:CPro)
sat_clay = SatSet(clay_src, type=:CLay)

timeindex = [2035:2065]
lidarprofile = TrackMatcher.get_lidarheights((15_000, -Inf), Float32)
cpro = CPro([joinpath(@__DIR__, "data", "caliop", "CPro", "CPro_4.h5")], timeindex, lidarprofile)
clay = CLay([joinpath(@__DIR__, "data", "caliop", "CLay", "CLay_4.h5")], timeindex, (15_000, -Inf))
