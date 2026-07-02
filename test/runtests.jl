using Test, TrackMatcher
using Dates, DataFrames, StructArrays, HDF5, CSV

include("test_flightdata.jl")
include("test_satdata.jl")
include("test_lidar.jl")
include("test_dataprocessing.jl")
