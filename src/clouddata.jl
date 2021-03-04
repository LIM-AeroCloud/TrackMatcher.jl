### Routines related to loading cloud track data

"""
    loadCloudTracks(files::String...; Float::DataType=Float32)

Load cloud track data from mat `files` and store in Julia format as `CloudTrack` structs.
Set the floating point precision to `Float` (default: `Float32`).
"""
function loadCloudTracks(files::String...; Float::DataType=Float32)
  # Open MATLAB session
  ms = mat.MSession(0)
  # Initialise
  tracks = CloudData[]
  # Loop over all mat files
  for (i, file) in enumerate(files)
    # Read data from mat files
    t, lonlat = readMAT(ms, file)
    # Store data in Julia format as Vector of CloudTrack structs
    storeMAT!(tracks, t, lonlat, i, file; Float=Float)
  end #loop over files
  mat.close(ms)

  return tracks
end #function loadCloudTracks


"""
    readMAT(ms::mat.MSession, file::String)

Use MATLAB session `ms` to read cloud track data from a mat `file`.
"""
function readMAT(ms::mat.MSession, file::String)
  # Send file name to MATLAB
  mat.put_variable(ms, :file, file)
  # Load Data from mat file using MATLAB
  mat.eval_string(ms, "data = load('-mat', file)")
  mat.eval_string(ms,
  """
  t = cell(length(data.cloud), 1)
  coords = cell(length(data.cloud), 1)
  for i = 1:length(data.cloud)
    t{i} = data.cloud(i).timestamp
    coords{i} = data.cloud(i).centrLonLat
  end
  """)
  # Load data from MATLAB into TrackMatcher/Julia
  t = mat.jarray(mat.get_mvariable(ms, :t))
  lonlat = mat.jarray(mat.get_mvariable(ms, :coords))
  return t, lonlat
end


"""
    storeMAT!(
      tracks::Vector{CloudTrack},
      t::Array,
      lonlat::Array;
      Float::DataType=Float32
    )

Append the vector with cloud `tracks` by timestamps `t` and coordinates `lonlat`
using the floating point precision set by `Float` (default: `Float32`) for the
positional data.
"""
function storeMAT!(
  tracks::Vector{CloudData},
  t::Array,
  lonlat::Array,
  fileID::Int,
  filename::String;
  Float::DataType=Float32
)
  # Transform MATLAB data into Julia Format and store as CloudTrack struct in a vector
  for i = 1:length(t)
    data = DataFrame(time = vec(DateTime.(t[i], "yyyymmddHHMM")),
      lat = Float.(lonlat[i][:,2]), lon = Float.(lonlat[i][:,1]))
    # Determine predominant trajectory direction, inflection points, and remove duplicate entries
    flex, useLON = preptrack!(data)
    isempty(flex) && continue
    push!(tracks, CloudData{Float}(data, CloudMetadata{Float}(string(fileID, ".", i), data,
      flex, useLON, filename)))
  end
end #function storeMAT!
