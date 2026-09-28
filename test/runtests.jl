using Test, TrackMatcher, Logging
using Dates, TimeZones, DataFrames, StructArrays, HDF5, CSV
import IntervalArithmetic.Symbols: (..)

global_logger(ConsoleLogger(stderr, Error))
# Always print @info messages issued directly in this file, regardless of the global logger level
testing_logger = ConsoleLogger(stderr, Info)

haskey(ENV, "TRACKMATCHER_PROGRESS") || (ENV["TRACKMATCHER_PROGRESS"] = "false")

with_logger(testing_logger) do
    @info "Setting up test data"
end
include("init.jl")
include("setup.jl")

with_logger(testing_logger) do
    @info "Running tests on data imports"
end
include("test_flightdata.jl")
include("test_clouddata.jl")
include("test_satdata.jl")

with_logger(testing_logger) do
    @info "Running tests on data processing"
end
include("test_lidar.jl")
include("test_dataprocessing.jl")

with_logger(testing_logger) do
    @info "Running tests on intercept finding"
end
include("test_intercepts.jl")
