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

Scan `folder` recursively for files of `filetypes` and add to the `inventory`.
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
end # function findfiles!


"""
    standardisecols!(df::DataFrame)

Convert possible SentinelArrays from `CSV` to `Vector{<:Union{Missing, Type}}`.
"""
function standardisecols!(df::DataFrame)
  for col in names(df)
    df[!,col] = [el for el in df[!,col]]
    df[!,col] isa Vector{Bool} && (df[!,col] = BitVector(df[!,col]))
  end #loop over DataFrame columns
end #function standardisecols


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
    while j ≤ iEnd && data[i, x] ≈ data[j, x]
      if data[j-1, y] ≈ data[j, y]
        # Delete datarow, if x and y are identical
        delete!(data, j-1)
        # Decrease the counter for the end of the arrays
        iEnd -= 1
      else
        # If x values are equal but y values differ,
        # add infinitessimal small increment to x value
        data[j, x] = nextfloat(data[j-1, x])
        j += 1 # set to next index
      end
    end
    i += 1 # increase line counter
  end
end #function remdup!


"""
    findflex(x::Vector{<:Real}) -> Vector{ NamedTuple{(:range, :min, :max), Tuple{UnitRange, AbstractFloat, AbstractFloat}}}

Find inflection points in `x` and return a vector of named tuples with index ranges
between inflection points and corresponding extrema.
"""
function findflex(x::Union{Vector{<:Real},CSV.SentinelArrays.SentinelArray{<:Real}})
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
    checkcols!(
      data::DataFrame,
      standardnames::Vector{String},
      standardtypes::Vector{<:Type},
      bounds::Tuple{Vararg{Pair{<:Union{Int,Symbol},<:Tuple}}},
      dataset::T where T<:Union{Nothing,<:AbstractString}=nothing,
      id::Union{Nothing,Int,AbstractString}=nothing;
      essentialcols::Vector{Int}=[1,2,3]
    ) -> data::DataFrame

Check the columns of the `data` `DataFrame` for the use of correct `standardnames`,
`standardtypes`, and `bounds`. Throw an error for deviations in `essentialcols`,
otherwise issue a warning stating the `dataset` and `id` and try to correct the `data`.
"""
function checkcols!(
  data::DataFrame,
  standardnames::Vector{String},
  standardtypes::Vector{<:Type},
  bounds::Tuple{Vararg{Pair{<:Union{Int,Symbol},<:Tuple}}},
  dataset::T where T<:Union{Nothing,<:AbstractString}=nothing,
  id::Union{Nothing,Int,AbstractString}=nothing;
  essentialcols::Vector{Int}=[1,2,3]
)
  # Warn of non-standardised data
  check = (names(data) == standardnames && all([typeof(d) <: t for (d, t) in zip(df.eachcol(data), standardtypes)])) ||
    @warn "Non-standard names and/or order used for data columns. Trying to correct..." dataset id

  # Bring column bounds into the right format
  colbounds = definebounds(bounds, standardnames)
  # Setup vector holding checked and correct data column indices
  correctcols = zeros(Int, length(standardnames))
  # Find columns by name
  findbyname!(correctcols, data, standardnames, standardtypes, colbounds, check)
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


"""
    definebounds(
      bounds::Tuple{Vararg{Pair{<:Union{Int,Symbol},<:Tuple}}},
      colnames::Vector{Symbol}
    ) -> bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}}

From a Tuple of Pairs with column names or indices pointing to tuples with min/max
bounds, return a vector `bounds` of tuples with the min/max bounds corresponding
to each of the `colnames`. If `bounds` are not defined, they are set to `-Inf`/`Inf`
in the output vector.
"""
function definebounds(
  bounds::Tuple{Vararg{Pair{<:Union{Int,Symbol},<:Tuple}}},
  colnames::Vector{String}
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


"""
    checkbounds!(
      correctcols::Vector{Int},
      bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}},
      data::DataFrame,
      col::T where T <: Union{Int,Symbol},
      pos::Int,
      val::Int
    )
Check the `bounds` of a `col`umn in `data` and add the index `val`, if correct,
to `correctcols`. The index is added in `correctcols` at the `pos`ition, the column
will have in the final corrected DataFrame.
"""
function checkbounds!(
  correctcols::Vector{Int},
  bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}},
  data::DataFrame,
  col::T where T <: Union{Int,String},
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


"""
    findbyname!(
      correctcols::Vector{Int},
      data::DataFrame,
      standardnames::Vector{String},
      standardtypes::Vector{<:Type},
      bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}},
      check::T where T<:Union{Nothing,Bool}
    )

Find columns by `standardnames` in `data` and add the current position as `Int`
at the correct position in `correctcols`, if the type corresponds to the
`standardtypes` and `bounds` are correct. Issue a warning, if the previous `check`
found a mismatch of the column order, but all `data` columns could be identified by
`standardnames`.
"""
function findbyname!(
  correctcols::Vector{Int},
  data::DataFrame,
  standardnames::Vector{String},
  standardtypes::Vector{<:Type},
  bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}},
  check::T where T<:Union{Nothing,Bool}
)
  for (i, name) in enumerate(standardnames)
    col = findfirst(isequal(name), names(data))
    col ≠ nothing && typeof(data[!,name]) <: standardtypes[i] &&
      checkbounds!(correctcols, bounds, data, name, i, col)
  end
  isempty(findall(isequal(0), correctcols)) && check == nothing &&
    @warn "all columns corrected based on column names"
end #function findbyname!


"""
    findbyposition!(
      correctcols::Vector{Int},
      opencols::Vector{Int},
      data::DataFrame,
      standardtypes::Vector{<:Type},
      bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}}
    )

Find columns by the `standardtypes` and `bounds` at each column index in `data`
and add the current in `correctcols` for positive matches. Issue a warning,
if all data could be corrected by type (and previously name).
"""
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


"""
    findbytype!(
      correctcols::Vector{Int},
      opencols::Vector{Int},
      remianingcols::Vector{Int},
      standardtypes::Vector{<:Type},
      bounds::Vector{Tuple{<:Union{Real,DateTime},<:Union{Real,DateTime}}}
    )

Find the missing columns `opencols` in `data` using the `remainingcols` in `data`
which have not yet been assigned and add the current position as `Int`
at the correct position in `correctcols`, if `standardtypes` and `bounds` are
correct in a column. Issue a warning, if all data could be corrected by type
(and previously name and position).
"""
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


"""
    correctDF!(
      data::DataFrame,
      correctcols::Vector{Int},
      standardnames::Vector{String},
      essentialcols::Vector{Int}
    )

Correct the column order and/or names in `data` based on the order of indices in
`correctcols` and the `standardnames`. Throw an error for missing columns in
`essentialcols` otherwise issue warnings for missing columns filled with `missing`
values or additional columns being deleted.
"""
function correctDF!(data::DataFrame, correctcols::Vector{Int},
  standardnames::Vector{String}, essentialcols::Vector{Int})
  for (i, col) in enumerate(findall(isequal(0), correctcols))
    data[!,Symbol("missing$i")] = [missing for i = 1:length(data[!,1])]
    correctcols[col] = length(names(data))
    col in essentialcols ? throw(@error "column $(standardnames[col]) not found in data") :
      @warn("column $(standardnames[col]) not found in data; filled with `missing`")
  end
  additionalcols = setdiff(collect(1:length(names(data))), correctcols)
  !isempty(additionalcols) &&
    @warn "additional column $(join(additionalcols, ", ", " and ")) deleted in data"
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
    find_timespan(sat::DataFrame, X::Tuple{<:AbstractFloat, <:AbstractFloat}, dataspan::Int=15)
      -> DataFrame, Vector{Int}

From the `sat` data in a `DataFrame` and the intersection `X` (as lat/lon pair),
find the time indices `t` ± `dataspan` for which the distance at `t` is minimal to `X`.

Return a `DataFrame` with `sat` data in the time span together with a `Vector{Int}`
holding the file indices of the corresponding granule file(s).
The `sat` data may be smaller than the `dataspan` at the edges of the `sat` `DataFrame`.
"""
function find_timespan(sat::DataFrame, X::Tuple{<:AbstractFloat, <:AbstractFloat},
  dataspan::Int=15)
  # Find index in sat data array with minimum distance to analytic intersection solutin
  Xlatlon = geo.LatLon(X...)
  coords = geo.LatLon.(sat.lat, sat.lon)
  imin = argmin([geo.distance(coord, Xlatlon) for coord in coords])
  # Find first/last index of span acknowledging bounds of the data array
  t1 = max(1, min(imin-dataspan, length(sat.time)))
  t2 = min(length(sat.time), imin+dataspan)

  return sat.time[t1:t2], unique(sat.fileindex[t1:t2])
end #function find_timespan


"""
    extract_timespan(sat::Union{CLay,CPro}, timespan::Vector{DateTime}) -> T where T<:Union{CLay,CPro}

From the `sat` data of type `CLay` or `CPro`, extract a subset within `timespan`
and return the reduced struct.
"""
function extract_timespan(sat::Union{CLay,CPro}, timespan::Vector{DateTime})
  timeindex = [findfirst(sat.data.time .== t) for t in timespan
    if findfirst(sat.data.time .== t) ≠ nothing]
  satdata = sat.data[timeindex,:]
  typeof(sat) == CPro ? CPro(satdata) : CLay(satdata)
end #function extract_timespan


"""
    get_flightdata(flight::FlightData, X::Tuple{<:AbstractFloat, <:AbstractFloat}, flightspan::Int)
      -> flightdata::FlightData, index::Int

From the measured `flight` data and lat/lon coordinates the intersection `X`,
save the closest measured value to the interpolated intersection ±`flightspan` data points
to `flightdata` and return it together with the `index` in `flightdata` of the time with the coordinates
closest to `X`.
"""
function get_flightdata(flight::FlightData, X::Tuple{<:AbstractFloat, <:AbstractFloat}, flightspan::Int)
  # Convert to LatLon structs
  Xlatlon = geo.LatLon(X...)
  coords = geo.LatLon.(flight.data.lat, flight.data.lon)
  # Find the index (DataFrame row) of the intersection in the flight data
  imin = argmin([geo.distance(coord, Xlatlon) for coord in coords])
  # Construct FlightData at Intersection
  t1 = max(1, min(imin-flightspan, length(flight.data.time)))
  t2 = min(length(flight.data.time), imin+flightspan)
  # t2 = t-timespan > length(data[:,1]) ? 0 : t2
  flightdata = FlightData(flight.data[t1:t2,:], flight.metadata)

  return flightdata, argmin([geo.distance(coord, Xlatlon)
    for coord in geo.LatLon.(flightdata.data.lat, flightdata.data.lon)])
end #function get_flightdata


"""
    get_satdata(
      ms::mat.MSession,
      sat::SatData,
      X::Tuple{<:AbstractFloat, <:AbstractFloat},
      satspan::Int,
      flightalt::Real,
      flightid::Union{Int,String},
      lidarprofile::NamedTuple,
      lidarrange::Tuple{Real,Real},
      savesecondtype::Bool
    ) -> cpro::CPro, clay::CLay, feature::Symbol, ts::Int

Using the `sat` data measurements within the overlap region and the MATLAB session
`ms`, extract CALIOP cloud profile (`cpro`) and/or layer data (`clay`) together with
the atmospheric `feature` at flight level (`flightalt`) for the data point closest
to the calculated intersection `X` ± `satspan` timesteps. In addition, return the
index `ts` within `cpro`/`clay` of the data point closest to `X`.
When `savesecondtype` is set to `false`, only the data type (`CLay`/`CPro`) in `sat`
is saved; if set to `true`, the corresponding data type is saved if available.
The lidar column data is saved for the height levels givin in the `lidarprofile` data
for the `lidarrange`.
"""
function get_satdata(
  ms::mat.MSession,
  sat::SatData,
  X::Tuple{<:AbstractFloat, <:AbstractFloat},
  satspan::Int,
  flightalt::Real,
  flightid::Union{Int,String},
  lidarprofile::NamedTuple,
  lidarrange::Tuple{Real,Real},
  savesecondtype::Bool
)
  # Retrieve DataFrame at Intersection ± 15 time steps
  timespan, fileindex = find_timespan(sat.data, X, satspan)
  primfiles = map(f -> get(sat.metadata.files, f, 0), fileindex)
  secfiles = if sat.metadata.type == :CPro && savesecondtype
    replace.(primfiles, "CPro" => "CLay")
  elseif sat.metadata.type == :CLay && savesecondtype
    replace.(primfiles, "CLay" => "CPro")
  else
    String[]
  end

  # Get CPro/CLay data from near the intersection
  clay = sat.metadata.type == :CLay ? CLay(ms, primfiles, lidarrange, flightalt) :
    CLay(ms, secfiles, lidarrange, flightalt)
  cpro = sat.metadata.type == :CPro ? CPro(ms, primfiles, timespan, lidarprofile) :
    CPro(ms, secfiles, timespan, lidarprofile)
  clay = extract_timespan(clay, timespan)
  cpro = extract_timespan(cpro, timespan)

  # Define primary data and index of intersection in primary data
  primdata = sat.metadata.type == :CPro ? cpro : clay
  ts = argmin([geo.distance(coord, geo.LatLon(X...))
    for coord in geo.LatLon.(primdata.data.lat, primdata.data.lon)])

  # Get feature classification
  feature = sat.metadata.type == :CPro ?
    atmosphericinfo(primdata, lidarprofile.fine, ts, flightalt, flightid) :
    atmosphericinfo(primdata, flightalt, ts)

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


"""Convert feet to kilometers"""
ft2km(ft::Union{Missing,Real}) = 0.0003048ft

"""Convert feet to meters"""
ft2m(ft::Union{Missing,Real}) = 0.3048ft

"""Convert feet/min to m/s"""
ftpmin2mps(ftpmin::Union{Missing,Real}) = 0.00508ftpmin

"""Convert knots to m/s"""
knot2mps(knot::Union{Missing,Real}) = 1.852knot/3.6
