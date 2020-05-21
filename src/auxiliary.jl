### Helper functions

"""
    convertUTC(t::AbstractFloat) -> DateTime

Convert the CALIOP Profile UTC time (`t`) to a `DateTime`.
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
  ms = round(Int, 1000(utc - 3600h - 60m - s))

  # Return a DateTime from date and time (h/m/s) with timezone UTC
  return DateTime(Dates.yearmonthday(d)..., h, m, s, ms)
end


"""
    findfiles!(inventory::Vector{String}, folder::String, filetypes::String...) -> inventory

Scan `folder` recursively for files of `filetype` and add to the `inventory`.
"""
function findfiles!(inventory::Vector{String}, folder::String, filetypes::String...)
  # Scan directory for files and folders and save directory
  dir = readdir(folder); path = abspath(folder)
  for file in dir
    # Save current directory/file
    cwd = joinpath(path, file)
    if isdir(cwd)
      # Step into subdirectories and scan them, too
      findfiles!(inventory, cwd, filetypes...)
    elseif any(endswith.(file, filetypes)) && !startswith(file, ".")
      # Save files of correct type
      push!(inventory, cwd)
    end
  end

  return inventory
end # function findfiles!


"""
    remdup!(data::DataFrame, useLON::Bool)

Remove entries with duplicate x and y (`lat`/`lon` or `lon`/`lat`) values from
`data` or increase x by an infinitessimal number if x data is identical, but y data
is not.
"""
function remdup!(data::DataFrame, useLON::Bool)
  # Define x and y data
  x, y = useLON ? (:lon, :lat) : (:lat, :lon)
  # Initialise
  i = 1
  iEnd = length(data[!,x])
  # Loop over entries in vector
  while i < iEnd
    j = i + 1 # index for next consecutive line
    while j ≤ iEnd && data[i, x] == data[j, x]
      # Define a infinitessimal small value δ which can be repeatedly applied via Δ
      δ = eps(data[i, x]); Δ = 0
      if data[i, y] == data[j, y]
        # Delete entries from all arrays with equal x and y data
        df.deleterows!(data, i)
        # Decrease the counter for the end of the arrays
        iEnd -= 1
      else
        # If x values are equal, but y values differ add multiples of δ repeatedly
        # by adding δ to Δ, and Δ to x.
        Δ += δ
        data[j, x] += Δ
        j += 1 # set to next index
      end
    end
    i += 1 # increase line counter
  end

  # Return revised data
  return data
end #function remdup!


"""
    findflex(x::Vector{<:Real}) -> Vector{ NamedTuple{(:range, :min, :max), Tuple{UnitRange, AbstractFloat, AbstractFloat}}}

Find inflection points in `x` and return a vector of named tuples with index ranges
between inflection points and corresponding extrema.
"""
function findflex(x::Vector{<:Real})
  # Save starting index
  flex = Int[1]
  # Loop over data points
  for i = 2:length(x)-1
    # Find and save inflection points
    if x[i-1] > x[i] < x[i+1] || x[i-1] < x[i] > x[i+1]
      push!(flex, i)
    end
  end
  # Save end index
  push!(flex, length(x))
  # Initialise a vector for ranges between flex points, and min/max values
  ranges = NamedTuple{(:range, :min, :max), Tuple{UnitRange, AbstractFloat, AbstractFloat}}[]
  # Loop over saved inflection points
  for i = 2:length(flex)
    # Evaluate ascending/descending order of data and save extrema and range,
    # if it range contains more than one data point
    xmin, xmax = x[flex[i-1]] < x[flex[i]] ? (x[flex[i-1]], x[flex[i]]) :
      (x[flex[i]], x[flex[i-1]])
    flex[i] - flex[i-1] > 1 && push!(ranges, (range=flex[i-1]:flex[i], min=xmin, max=xmax))
  end

  # Return all ranges between flex points together with the corresponding extrema
  return Tuple(ranges)
end #function findflex


"""
    Minterpolate(ms::mat.MSession, p::mat.MxArray) -> λ(i::Union{Real,Vector{AbstractFloat},StepRangeLen})

From a MATLAB piecewise polynomial structure `p` and the corresponding MATLAB
session `ms` return a λ function taking a Union{Real,Vector{AbstractFloat},StepRangeLen}
as input to return the interpolated data of `p`.
"""
function Minterpolate(ms::mat.MSession, p::mat.MxArray)
  """ Return the interpolated data of `i` using MATLAB's pchip algorithm"""
  function (i::Union{Real,Vector{<:AbstractFloat},StepRangeLen})
    # Send data to MATLAB
    mat.put_variable(ms, :i, mat.mxarray(i)); mat.put_variable(ms, :p, p)
    # Interpolate data with MATLAB
    mat.eval_string(ms, "pp = ppval(p,i);")
    # Retrieve interpolated data from MATLAB
    mat.jvalue(mat.get_mvariable(ms, :pp))
  end
end


function checkcols!(
  data::DataFrame,
  standardnames::Vector{Symbol},
  standardtypes::Vector{<:Type},
  bounds::Tuple{Vararg{Pair{<:Union{Int,Symbol},<:Tuple}}},
  dataset::T where T<:Union{Nothing,<:AbstractString}=nothing,
  id::Union{Nothing,Int,AbstractString}=nothing;
  essentialcols::Vector{Int}=[1,2,3]
)
  # Warn of non-standardised data
  names(data) == standardnames ? (return data) :
    (@warn "Non-standard names and/or order used for data columns. Trying to correct..." dataset id)

  # Bring column bounds into the right format
  colbounds = definebounds(bounds, standardnames)
  # Setup vector holding checked and correct data column indices
  correctcols = zeros(Int, length(standardnames))
  # Find columns by name
  findbyname!(correctcols, data, standardnames, standardtypes, colbounds)
  # Find unchecked columns
  opencols = setdiff(collect(1:min(length(names(data)), length(standardnames))),
    filter(x -> x .≠ 0, correctcols))
  # Find columns by position
  findbyposition!(correctcols, opencols, data, standardtypes, colbounds)
  # Find columns in data not associated with a correct column
  remainingcols = setdiff(collect(1:length(names(data))),
    filter(x -> x .≠ 0, correctcols))
  # Find unchecked columns
  opencols = findall(isequal(0), correctcols)
  # Find columns by type and bounds
  findbytype!(correctcols, opencols, remainingcols, standardtypes, colbounds)
  # Correct DataFrame to the right format
  correctDF!(data, correctcols, standardnames, essentialcols)
  return data
end #function checkcols!


function definebounds(
  bounds::Tuple{Vararg{Pair{<:Union{Int,Symbol},<:Tuple}}},
  colnames::Vector{Symbol}
)
  colbounds = Dict(bounds)
  bounds = Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}}(undef,length(colnames))
  for (i, name) in enumerate(colnames)
    lower, upper = if haskey(colbounds, name)
      colbounds[name]
    elseif haskey(colbounds, i)
      colbounds[i]
    else
      (-Inf, Inf)
    end
    bounds[i] = (lower, upper)
  end

  return bounds
end #function definebounds!


function checkbounds!(
  correctcols::Vector{Int},
  bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}},
  data::DataFrame,
  col::T where T <: Union{Int,Symbol},
  pos::Int,
  val::Int
)
  if typeof(bounds[pos]) .== Tuple{DateTime,DateTime} &&
    all(bounds[pos][1] .< skipmissing(data[!,col]) .< bounds[pos][2])
    correctcols[pos] = val
  elseif typeof(bounds[pos]) .== Tuple{DateTime,DateTime}
    return
  elseif all(isinf.(bounds[pos])) ||
    all(bounds[pos][1] .< skipmissing(data[!,col]) .< bounds[pos][2])
    correctcols[pos] = val
  end
end


function findbyname!(
  correctcols::Vector{Int},
  data::DataFrame,
  standardnames::Vector{Symbol},
  standardtypes::Vector{<:Type},
  bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}}
)
  for (i, name) in enumerate(standardnames)
    col = findfirst(isequal(name), names(data))
    col ≠ nothing && typeof(data[!,name]) <: standardtypes[i] &&
      checkbounds!(correctcols, bounds, data, name, i, col)
  end
  isempty(findall(isequal(0), correctcols)) &&
    @warn "all columns corrected based on column names"
end #function findbyname!

function findbyposition!(
  correctcols::Vector{Int},
  opencols::Vector{Int},
  data::DataFrame,
  standardtypes::Vector{<:Type},
  bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}}
)
  isempty(findall(isequal(0), correctcols)) && return
  for pos in opencols
    try typeof(data[!,pos]) <: standardtypes[pos] &&
      checkbounds!(correctcols, bounds, data, pos, pos, pos)
    catch
      continue
    end
  end

  isempty(findall(isequal(0), correctcols)) &&
    @warn "all columns corrected based on\n- column names\n- column positions"
end #function findbyposition!


function findbytype!(
  correctcols::Vector{Int},
  opencols::Vector{Int},
  remianingcols::Vector{Int},
  standardtypes::Vector{<:Type},
  bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}}
)
  isempty(findall(isequal(0), correctcols)) && return
  for i in opencols, j in remianingcols
    try typeof(data[!,j]) <: standardtypes[i] && correctcols[i] == 0 &&
      checkbounds!(correctcols, bounds, data, j, i, j)
    catch
      continue
    end
  end

  isempty(findall(isequal(0), correctcols)) &&
    @warn "all columns corrected based on\n- column names\n- column positions\n- column types and bounds"
end #function findbytype!


function correctDF!(data::DataFrame, correctcols::Vector{Int},
  standardnames::Vector{Symbol}, essentialcols::Vector{Int})
  for (i, col) in enumerate(findall(isequal(0), correctcols))
    data[!,Symbol("missing$i")] = [missing for i = 1:length(data[!,1])]
    correctcols[col] = length(names(data))
    col in essentialcols ? throw(@error "column $(standardnames[col]) not found in data") :
      @warn("column $(standardnames[col]) not found in data; filled with `missing`")
  end
  additionalcols = setdiff(collect(1:length(names(data))), correctcols)
  !isempty(additionalcols) &&
    @warn "additional columns $(join(additionalcols, ", ", " and ")) deleted in data"
  df.select!(data, correctcols)
  df.rename!(data, standardnames)
end #function correctDF!


"""
    checkDBtype(DB::Vector{FlightData}, type::String) -> Vector{FlightData}

Check whether the `DB` vector with `FlightData` contains only entries of the desired
`type` (`VOLPE AEDT inventory`, `FlightAware archive`, or `flightaware.com online data`).
Remove entries of a wrong type and return the revised vector.
"""
function checkDBtype(DB::Vector{FlightData}, type::String)
  removed = count([d.metadata.source.≠type for d in DB])
  removed > 0 && @warn "$removed entries removed of database of type $type."
  DB[[d.metadata.source.==type for d in DB]]
end #function checkDBtype


"""
    find_timespan(data::DataFrame, t::Int, timespan::Int=15) -> DataFrame

From the `data` in a `DataFrame`, extract a subset at index (row) `t` ± `timespan`
(rows).
"""
function find_timespan(sat::DataFrame, t::Int, timespan::Int=15)
  t1 = max(1, min(t-timespan, length(sat.time)))
  t2 = min(length(sat.time), t+timespan)
  # t2 = t-timespan > length(data[:,1]) ? 0 : t2
  return sat.time[t1:t2], unique(sat.fileindex[t1:t2])
end #function extract_timespan


function extract_timespan(sat::Union{CLay,CPro}, timespan::Vector{DateTime})
  timeindex = [findfirst(sat.data.time .== t) for t in timespan
    if findfirst(sat.data.time .== t) ≠ nothing]
  satdata = sat.data[timeindex,:]
  typeof(sat) == CPro ? CPro(satdata) : CLay(satdata)
end


"""
    swap_sattype(sattype::Symbol)::Symbol

Returns the other Symbol from the options

- `:CLay`
- `:CPro`
"""
function swap_sattype(sattype::Symbol)::Symbol
  swap = Dict{Symbol,Symbol}(:CPro => :CLay, :CLay => :CPro)
  swap[sattype]
end


"""
    get_trackdata(flight::FlightData, tflight::DateTime, flightspan::Int)
      -> flightdata::FlightData

From the measured `flight` data and the time of the aircrafat the intersection (`tflight`),
save the closest measured value to the interpolated intersection ±`flightspan` data points
to `flightdata`.
"""
function get_flightdata(flight::FlightData, tflight::DateTime, flightspan::Int)
  # Find the index (DataFrame row) of the intersection in the flight data
  tf = argmin(abs.(flight.data.time .- tflight))
  # Construct FlightData at Intersection
  t1 = max(1, min(tf-flightspan, length(flight.data.time)))
  t2 = min(length(flight.data.time), tf+flightspan)
  # t2 = t-timespan > length(data[:,1]) ? 0 : t2
  flightdata = FlightData(flight.data[t1:t2,:], flight.metadata)

  return flightdata, argmin(abs.(flightdata.data.time .- tflight))
end


function get_satdata(ms::mat.MSession, sat::SatData, tsat::DateTime, satspan::Int,
  flightalt::AbstractFloat, lidar::NamedTuple, lidarrange::Tuple{Real,Real},
  savesecondtype::Bool)
  # Get satellite data used to find the intersection and find DataFrame row of intersection
  ts = argmin(abs.(sat.data.time .- tsat))
  sattype = sat.metadata.type
  # Retrieve DataFrame at Intersection ± 15 time steps
  timespan, fileindex = find_timespan(sat.data, ts, satspan)
  primfiles = map(f -> get(sat.metadata.files, f, 0), fileindex)
  secfiles = if sattype == :CPro && savesecondtype
    replace.(primfiles, "CPro" => "CLay")
  elseif sattype == :CLay && savesecondtype
    replace.(primfiles, "CLay" => "CPro")
  else
    String[]
  end

  # Get CPro/CLay data from near the intersection
  clay = sattype == :CLay ? CLay(ms, primfiles) : CLay(ms, secfiles)
  cpro = sattype == :CPro ? CPro(ms, primfiles, timespan, lidar) :
    CPro(ms, secfiles, timespan, lidar)
  clay = extract_timespan(clay, timespan)
  cpro = extract_timespan(cpro, timespan)

  # Define primary data and index of intersection in primary data
  primdata = sattype == :CPro ? cpro : clay
  ts = argmin(abs.(primdata.data.time .- tsat))

  # Get feature classification
  feature = atmosphericinfo(primdata, lidar.fine, flightalt, ts)

  return cpro, clay, feature, ts
end #function get_satdata


"""
    timesec(t::Real) -> DateTime

From a UNIX DateTime, return a DateTime rounded to the second
"""
timesec(t::Real) = round(Dates.unix2datetime(t), Dates.Second)


"""
    timesec(t::Union{DateTime,ZonedDateTime}) -> DateTime

Round a `DateTime` or `ZonedDateTime` to the second.
"""
timesec(t::Union{DateTime,ZonedDateTime}) = round(t, Dates.Second)


"""
    preptrack(flight::DataFrame) -> flight, flex, useLON

Use the `flight` data to prepare the track for interpolation. Find the predominant
flight direction and inflection (flex) points in the flight track.

Return a tidied flight data from which duplicate entries are removed together
with the `flex` points in the x data for interpolation and a boolean `useLON`,
which is true for longitude values used x data in the track interpolation.
"""
function preptrack(flight::DataFrame)
  # calculate area covered by flight
  lp = any(flight.lon .≥ 0) ? maximum(filter(l -> l ≥ 0, flight.lon)) -
    minimum(filter(l -> l ≥ 0, flight.lon)) : 0
  ln = any(flight.lon .< 0) ? maximum(filter(l -> l < 0, flight.lon)) -
    minimum(filter(l -> l < 0, flight.lon)) : 0
  # Determine main direction of flight (N<>S, E<>W) and use it as x values
  # for flight interpolation (info stored as bool useLON)
  useLON = maximum(flight.lat) - minimum(flight.lat) ≤ (lp + ln) *
    cosd(stats.mean(flight.lat)) ? true : false
  # Adjust for duplicate entries and
  remdup!(flight, useLON)
  # find flex points to cut data in segments needed for the interpolation
  flex = useLON ? findflex(flight.lon) : findflex(flight.lat)

  return flight, flex, useLON
end #function preptrack


"""Convert feet to kilomenters"""
ft2km(ft::Real) = 0.0003048ft
