"""
    findtextfiles(folder::String)

Load inventory files saved in `folder` to the `inventory` holding the file names
and locations as a vector of strings.
"""
function findtextfiles(inventory::Vector{String}, folder::String)
  # Scan directory for files and folders and save directory
  dir = readdir(folder); path = abspath(folder)
  for file in dir
    cwd = joinpath(path, file)
    if endswith(file, ".txt") || endswith(file, ".dat")
      push!(inventory, cwd)
    elseif isdir(cwd)
      inventory = findcsv(inventory, cwd)
    end
  end

  return inventory
end # function findcsv
