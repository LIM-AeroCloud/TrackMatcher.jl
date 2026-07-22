using Test, TrackMatcher, Logging
using Dates, TimeZones, DataFrames, StructArrays, HDF5, CSV

global_logger(ConsoleLogger(stderr, Error))

haskey(ENV, "TRACKMATCHER_PROGRESS") || (ENV["TRACKMATCHER_PROGRESS"] = "false")

include("test_flightdata.jl")
include("test_clouddata.jl")
include("test_satdata.jl")
include("test_lidar.jl")
include("test_dataprocessing.jl")
