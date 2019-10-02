function checklength(vect, ref)
  len = length(ref) - length(vect)
  if len > 0
    @warn string("$(len) entries missing in vector compared to reference. ",
      "Missing entries are filled with `missing` at the end of the vector.")
    vect = [vect; [missing for i = 1:len]]
  elseif len < 0
    @warn string("$(-len) additional entries found in vector compared to reference. ",
      "Additional entries at the end of the vector are ignored.")
    vect = vect[1:length(ref)]
  end

  return vect
end


function convertUTC(t::AbstractFloat)
  date = floor(Int, t)
  d = Date("20"*string(date), "yyyymmdd")
  utc = (t - date) * 86400
  h = floor(Int, utc/3600)
  m = floor(Int, utc - 3600h)÷60
  s = floor(Int, utc - 3600h - 60m)

  return ZonedDateTime(DateTime(d, Time(h,m,s)), tz.tz"UTC")
end


"""
    findFiles(inventory::Vector{String}, folder::String, filetypes::String...) -> inventory

Load inventory files saved in `folder` to the `inventory` holding the file names
and locations ending with `filetypes` as a vector of strings.
"""
function findFiles(inventory::Vector{String}, folder::String, filetypes::String...)
  # Scan directory for files and folders and save directory
  fileendings = Regex(join(filetypes,'|'))
  dir = readdir(folder); path = abspath(folder)
  for file in dir
    cwd = joinpath(path, file)
    if endswith(file, fileendings)
      push!(inventory, cwd)
    elseif isdir(cwd)
      inventory = findFiles(inventory, cwd, filetypes...)
    end
  end

  return inventory
end # function findcsv
