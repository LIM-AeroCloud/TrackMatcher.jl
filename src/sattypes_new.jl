rootdir, type = "../ICAREdata/", :undef

files = scandir(rootdir, ".h5")
type = sat_datafiles!(files, type)
files

folders = ["../ICAREdata/05kmCPro.v4.51/2006/", "../ICAREdata/05kmCPro.v4.51/2019/"]
paths = []
for folder in folders
    files = scandir(folder, ".h5")
    type = sat_datafiles!(files, type)
    push!(paths, (;root = folder, files))
end
paths

## structs for satellite data

struct SatData{T} <: SatTrack{T}
    data::DataFrame

    """ Unmodified constructor for `SatData` with basic checks for correct `data`"""
    function SatData{T}(data::DataFrame) where T
        # Ensure floats of correct precision
        convert_floats!(data, T)
        # Check for correct column names and data types
        standardnames = ["time", "lat", "lon"]
        standardtypes = [Vector{DateTime}, Vector{T}, Vector{T}]
        bounds = (:time => (DateTime(2006), Dates.now()), :lat => (-90,90),
            :lon => (-180,180))
        checkcols!(data, standardnames, standardtypes, bounds, "CLay")
        new{T}(data)
    end #constructor 1 SatData
end

struct SecondaryMetadata{T} <: SecondarySet{T}
  granules::DataFrame
  type::Symbol
  date::NamedTuple{(:start,:stop),Tuple{DateTime,DateTime}}
  created::Union{DateTime,ZonedDateTime}
  loadtime::Dates.CompoundPeriod
  remarks
end #struct SecondaryMetadata

struct SatSet{T} <: SecondarySet{T}
    granules::Vector{SatData{T}}
    metadata::SecondaryMetadata{T}

    """ Unmodified constructor for `SatData` with basic checks for correct `data`"""
    function SatSet{T}(granules::Vector{SatData{T}}, metadata::SecondaryMetadata{T}) where T
        # Ensure floats of correct precision
        convert_floats!.(granules, T)
        # Instantiate struct
        new{T}(granules, metadata)
    end #constructor 1 SatData
end

function SatSet{T}(
    folders::String...;
    type::Symbol = :undef,
    attachments = nothing
)#::SatSet{T} where T<:AbstractFloat
    # Scan folders for satellite data files
    paths = []
    for folder in folders
        files = scandir(folder, ".h5")
        type = sat_datafiles!(files, type)
        push!(paths, (root = folder; files))
    end
    for file in files
        println(file)
    end
    # Return SatSet constructor
    return #SatSet{T}()
end
