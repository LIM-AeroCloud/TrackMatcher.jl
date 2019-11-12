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
    convertUTC(t::Float64) -> ZonedDateTime

Convert the CALIOP Profile UTC time (`t`) to a `ZonedDateTime` with `TimeZone` `UTC`.
"""
function convertUTC(t::Float64)
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


function remdup(x::Vector{<:Float64}, y::Vector{<:Float64},
  alt::Vector{<:Float64}, speed::Vector{<:Float64}, t::Vector{<:ZonedDateTime})
  i = 1
  iEnd = length(x)
  while i < iEnd
    j = i + 1
    while j ≤ iEnd && x[i] == x[j]
      δ = eps(x[i]); Δ = 0
      if y[i] == y[j]
        deleteat!(x, i); deleteat!(y, i); deleteat!(alt, i); deleteat!(speed, i)
        deleteat!(t, i)
        iEnd -= 1
      else
        Δ += δ
        x[j] += Δ
        j += 1
      end
    end
    i += 1
  end

  return x, y, alt, speed, t
end


function findFlex(x::Vector{<:Real})
  flex = Int[1]
  for i = 2:length(x)-1
    if x[i-1] > x[i] < x[i+1] || x[i-1] < x[i] > x[i+1]
      if count(isequal(-180), x[i-1:i+1]) == 1 &&
        !(sign(x[i-1]) == sign(x[i+1]) && x[i] == -180)
        continue
      else
        push!(flex, i)
      end
    end
  end
  push!(flex, length(x))
  ranges = UnitRange[]
  for i = 2:length(flex)
    push!(ranges, flex[i-1]:flex[i])
  end

  return ranges
end


function Minterpolate(ms::mat.MSession, p::mat.MxArray)
  function (i::Union{Real,Vector{<:AbstractFloat},StepRangeLen})
    mat.put_variable(ms, :i, mat.mxarray(i)); mat.put_variable(ms, :p, p)
    mat.eval_string(ms, "pp = ppval(p,i);")
    mat.jvalue(mat.get_mvariable(ms, :pp))
  end
end
