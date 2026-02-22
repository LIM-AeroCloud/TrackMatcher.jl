## Define expected results

cpro_files = [joinpath("Level1", "CPro1.h5"), joinpath("Level1", "Level2", "CPro2.h5")]
clay_files = [joinpath("Level1", "CLay.h5")]
mixed_files = [joinpath("Level1", "CLay.h5"), joinpath("Level1", "CPro1.h5"),
    joinpath("Level1", "Level2", "CPro2.h5")]
empty_files = String[]


## Helper functions to evaluate tests

"""
    test_sat_datafiles(files, type, satfiles, expected_type=type) -> Bool

Test function `TrackMatcher.sat_datafiles` to return the `expected_type` from the given `type`
and clean `files` to match `satfiles`.
Returns `true` if both tests pass, `false` otherwise.
"""
function test_sat_datafiles(files, type, satfiles, expected_type=type)::Bool
    success = false
    input_files = copy(files)
    sattype = TrackMatcher.sat_datafiles!(input_files, type)
    success |= sattype == expected_type
    success |= input_files == satfiles
    return success
end


## Testsets

@testset "read satellite data" begin
    @test TrackMatcher.scandir(joinpath("data", "caliop", "correct"), ".h5") == mixed_files
    @test TrackMatcher.scandir(joinpath("data", "caliop", "correct"), [".h5"]) == mixed_files
    @test TrackMatcher.scandir(joinpath("data", "caliop", "correct"), [".h5", ".hdf"]) == mixed_files
    @test isempty(TrackMatcher.scandir(joinpath("data", "caliop", "no_satdata"), ".h5"))
    @test @test_logs(
        (:info, "satellite data type auto-detected as CPro"),
        (:warn, "the given folder contains 1 non-CPro files which will be ignored"),
        test_sat_datafiles(mixed_files, :undef, cpro_files, :CPro)
    )
    @test @test_logs(
        (:warn, "the given folder contains 1 non-CPro files which will be ignored"),
        test_sat_datafiles(mixed_files, :CPro, cpro_files)
    )
    @test @test_logs(
        (:warn, "the given folder contains 2 non-CLay files which will be ignored"),
        test_sat_datafiles(mixed_files, :CLay, clay_files)
    )
    @test test_sat_datafiles(cpro_files, :CPro, cpro_files)
    @test @test_logs(
        (:warn, "the given folder contains 2 non-CLay files which will be ignored"),
        (:warn, "no CLay data files found, TrackMatcher requires CALIOP data"),
        test_sat_datafiles(cpro_files, :CLay, empty_files)
    )
    @test @test_logs(
        (:warn, "no satellite data files found, TrackMatcher requires CALIOP data"),
        test_sat_datafiles(empty_files, :undef, empty_files, :undef)
    )
end

@testset "I/O checks" begin
    # TODO @test checkcols!
end
