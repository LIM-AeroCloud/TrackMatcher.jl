@testset "cloud data" begin
    # Test data
    t = [
        DateTime(2012, 2, 6, 0, 15, 0), DateTime(2012, 2, 6, 0, 30, 0),
        DateTime(2012, 2, 6, 0, 45, 0), DateTime(2012, 2, 6, 1, 0, 0)
    ]
    lat = [6.2779856, 5.648815, 4.928818, 4.211443]
    lon = [23.687, 23.552141, 23.398176, 23.249731]
    # Instantiate test sets
    cloud = CloudSet(joinpath(@__DIR__, "data", "cloud"))
    cloud64 = CloudSet{Float64}(joinpath(@__DIR__, "data", "cloud"))
    cloud16 = CloudSet{Float16}(cloud64)
    t0 = now()
    cloud_empty = CloudSet()
    cloud_nothing = CloudSet(joinpath(@__DIR__, "data", "archive"))
    # Test data
    @testset "data integrity" begin
        @test_throws Base.IOError CloudSet(joinpath(@__DIR__, "data", "foo"))
        @test length(cloud.tracks) == 3
        @test cloud.tracks.time[1] isa Vector{<:DateTime}
        @test cloud.tracks.lat[1] isa Vector{<:Float32}
        @test cloud.tracks.lon[1] isa Vector{<:Float32}
        @test cloud.tracks.time[1] == t
        @test cloud.tracks.lat[1] ≈ lat
        @test cloud.tracks.lon[1] ≈ lon
        @test cloud64.tracks.time[1] isa Vector{<:DateTime}
        @test cloud64.tracks.lat[1] isa Vector{<:Float64}
        @test cloud64.tracks.lon[1] isa Vector{<:Float64}
        @test cloud64.tracks.time[1] == t &&
            cloud64.tracks.lat[1] ≈ lat && cloud64.tracks.lon[1] ≈ lon
        @test cloud16.tracks.time[1] isa Vector{<:DateTime}
        @test cloud16.tracks.lat[1] isa Vector{<:Float16}
        @test cloud16.tracks.lon[1] isa Vector{<:Float16}
        @test cloud16.tracks.time[1] == t &&
            cloud16.tracks.lat[1] ≈ lat && cloud16.tracks.lon[1] ≈ lon
        @test cloud_empty.tracks isa StructArray{CloudData{Float32}}
        @test cloud_nothing.tracks isa StructArray{CloudData{Float32}}
        @test isempty(cloud_empty.tracks)
        @test isempty(cloud_nothing.tracks)
        @test t0 ≤ cloud_empty.metadata.date.start == cloud_empty.metadata.date.stop ≤
            DateTime(cloud_empty.metadata.created)
        @test t0 ≤ cloud_nothing.metadata.date.start == cloud_nothing.metadata.date.stop ≤
            DateTime(cloud_nothing.metadata.created)
    end
    @testset "constructors" begin
        # Define datasets
        primary = PrimarySet(joinpath(@__DIR__, "data", "cloud"))
        primary64 = PrimarySet{Float64}(joinpath(@__DIR__, "data", "cloud"))
        cloudprimary = PrimarySet{Float64}(primary)

        fields, fields64 = [], []
        for field in fieldnames(CloudData)
            push!(fields, getproperty(cloud.tracks[1], field))
            push!(fields64, getproperty(cloud64.tracks[1], field))
        end
        track = CloudTrack(fields...)
        track64 = CloudTrack{Float64}(fields64...)
        track_converted = CloudTrack{Float64}(track)
        track_empty = CloudTrack()
        cloudtrack = CloudData(fields...)
        primmeta = PrimaryMetadata()
        cloudmeta = CloudMetadata()

        # Test constructors and datasets
        @test primary isa CloudSet{Float32}
        @test primary64 isa CloudSet{Float64}
        @test cloudprimary isa CloudSet{Float64}
        @test track.lat == Float32.(lat)
        @test track.lon == Float32.(lon)
        @test track64.time == t
        @test track64.lat ≈ lat
        @test track64.lon ≈ lon
        @test track_converted.time == t
        @test track_converted.lat ≈ lat atol = 1e-5
        @test track_converted.lon ≈ lon atol = 1e-5
        @test track_empty isa CloudData{Float32}
        @test isempty(track_empty.time) && isempty(track_empty.lat) &&
            isempty(track_empty.lon)
        @test cloudtrack isa CloudData{Float32}
        @test cloudtrack.lat ≈ lat
        @test cloudtrack.lon ≈ lon
        @test cloudtrack.lat isa Vector{Float32} && cloudtrack.lon isa Vector{Float32}
        @test primmeta isa PrimaryMetadata{Float32}
    end
end
