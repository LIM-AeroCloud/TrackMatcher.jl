"""
    checklength(vect, ref) -> vect

Check that length of vector `vect` is the same as reference vector `ref`,
otherwise warn and either fill missing entries with `missing` at the end or
delete additional entries at the end of `vect` and return `vect`.
"""
function checklength(vect, ref)
  # Compare lengths of vectors
  len = length(ref) - length(vect)
  if len > 0
    # Warn of vector shorter than reference and fill with missing
    @warn string("$(len) entries missing in vector compared to reference. ",
      "Missing entries are filled with `missing` at the end of the vector.")
    vect = [vect; [missing for i = 1:len]]
  elseif len < 0
    # Warn of vector longer than reference and delete additonal entries
    @warn string("$(-len) additional entries found in vector compared to reference. ",
      "Additional entries at the end of the vector are ignored.")
    vect = vect[1:length(ref)]
  end

  # Return (modified) vector
  return vect
end


"""
    convertUTC(t::AbstractFloat) -> ZonedDateTime

Convert the CALIOP Profile UTC time (`t`) to a `ZonedDateTime` with `TimeZone` `UTC`.
"""
function convertUTC(t::AbstractFloat)
  # Extract date from Float before decimal point and convert to Date
  date = floor(Int, t)
  d = Date("20"*string(date), "yyyymmdd")
  # Calculate overall time in seconds
  utc = (t - date) * 86400
  # From overall time, calculate hours, minutes, and seconds
  h = floor(Int, utc/3600)
  m = floor(Int, utc - 3600h)÷60
  s = floor(Int, utc - 3600h - 60m)

  # Return a ZonedDateTime from date and time (h/m/s) with timezone UTC
  return ZonedDateTime(DateTime(d, Time(h,m,s)), tz.tz"UTC")
end


"""
    findFiles(inventory::Vector{String}, folder::String, filetypes::String...) -> inventory

Load inventory files saved in `folder` to the `inventory` holding the file names
and locations ending with `filetypes` as a vector of strings.
"""
function findFiles(inventory::Vector{String}, folder::String, filetypes::String...)
  # Construct Regex of file endings from filetypes
  fileendings = Regex(join(filetypes,'|'))
  # Scan directory for files and folders and save directory
  dir = readdir(folder); path = abspath(folder)
  for file in dir
    # Save current directory/file
    cwd = joinpath(path, file)
    if endswith(file, fileendings)
      # Save files of correct type
      push!(inventory, cwd)
    elseif isdir(cwd)
      # Step into subdirectories and scan them, too
      inventory = findFiles(inventory, cwd, filetypes...)
    end
  end

  return inventory
end # function findcsv
