### Helper functions

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
end #function checklength


"""
    convertUTC(t::Float64) -> DateTime

Convert the CALIOP Profile UTC time (`t`) to a `DateTime`.
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

  # Return a DateTime from date and time (h/m/s) with timezone UTC
  return DateTime(Dates.yearmonthday(d)..., h, m, s)
end


"""
    findfiles!(inventory::Vector{String}, folder::String, filetypes::String...) -> inventory

Scan `folder` recursively for files of `filetype` and add to the `inventory`.
"""
function findfiles!(inventory::Vector{String}, folder::String, filetypes::String...)
  # Construct Regex of file endings from filetypes
  fileextensions = Regex(join(filetypes,'|'))
  # Scan directory for files and folders and save directory
  dir = readdir(folder); path = abspath(folder)
  for file in dir
    # Save current directory/file
    cwd = joinpath(path, file)
    if isdir(cwd)
      # Step into subdirectories and scan them, too
      findfiles!(inventory, cwd, filetypes...)
    elseif endswith(file, fileextensions) && !startswith(file, ".")
      # Save files of correct type
      push!(inventory, cwd)
    end
  end

  return inventory
end # function findfiles!


"""
    remdup!(x::Vector{<:Float64}, y::Vector{<:Float64},
      alt::Vector{<:Float64}, speed::Vector{<:Float64}, t::Vector{<:ZonedDateTime})

Remove entries with duplicate `x` and `y` (`lat`/`lon` or `lon`/`lat`) values from
these arrays as well as `alt`, `speed`, and `t`. Increase x by an infinitessimal
number if x data is identical, but y data is not.
"""
function remdup!(x::Vector{Float64}, y::Vector{Float64}, t::Vector{ZonedDateTime},
  alt::Vector{<:Union{Missing,Float64}}, speed::Vector{<:Union{Missing,Float64}},
  head::Vector{<:Union{Missing,Int}}, climb::Vector{<:Union{Missing,Int}})
  # Initialise
  i = 1
  iEnd = length(x)
  # Loop over entries in vector
  while i < iEnd
    j = i + 1 # index for next consecutive line
    while j ≤ iEnd && x[i] == x[j]
      # Define a infinitessimal small value δ which can be repeatedly applied via Δ
      δ = eps(x[i]); Δ = 0
      if y[i] == y[j]
        # Delete entries from all arrays with equal x and y data
        deleteat!(x, i); deleteat!(y, i); deleteat!(t, i)
        deleteat!(alt, i); deleteat!(speed, i); deleteat!(head, i); deleteat!(climb, i)

        # Decrease the counter for the end of the arrays
        iEnd -= 1
      else
        # If x values are equal, but y values differ add multiples of δ repeatedly
        # by adding δ to Δ, and Δ to x.
        Δ += δ
        x[j] += Δ
        j += 1 # set to next index
      end
    end
    i += 1 # increase line counter
  end

  # Return revised data
  return x, y, t, alt, speed, head, climb
end #function remdup!



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
    findflex(x::Vector{<:Real}) -> Vector{ NamedTuple{(:range, :min, :max), Tuple{UnitRange, Float64, Float64}}}

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
  ranges = NamedTuple{(:range, :min, :max), Tuple{UnitRange, Float64, Float64}}[]
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
    Minterpolate(ms::mat.MSession, p::mat.MxArray) -> λ(i::Union{Real,Vector{<:Float64},StepRangeLen})

From a MATLAB piecewise polynomial structure `p` and the corresponding MATLAB
session `ms` return a λ function taking a Union{Real,Vector{<:Float64},StepRangeLen}
as input to return the interpolated data of `p`.
"""
function Minterpolate(ms::mat.MSession, p::mat.MxArray)
  """ Return the interpolated data of `i` using MATLAB's pchip algorithm"""
  function (i::Union{Real,Vector{<:Float64},StepRangeLen})
    # Send data to MATLAB
    mat.put_variable(ms, :i, mat.mxarray(i)); mat.put_variable(ms, :p, p)
    # Interpolate data with MATLAB
    mat.eval_string(ms, "pp = ppval(p,i);")
    # Retrieve interpolated data from MATLAB
    mat.jvalue(mat.get_mvariable(ms, :pp))
  end
end


"""
    checkcols(data::DataFrame, standardnames::Vector{Symbol},
      standardtypes::Vector{<:Union{Union,DataType}}, bounds::Vector{Tuple{Real,Real}},
      dataset::T where T<:AbstractString, id::Union{Nothing,Int,AbstractString})
      -> DataFrame

Check `DataFrame` `data` for correct `standardnames`, order, and `standardtypes`
together with correct `bounds` of each column. Issue warnings on mismatches for
data giving the `dataset` and `id` for clarification.
"""
function checkcols(data::DataFrame, standardnames::Vector{Symbol},
  standardtypes::Vector{<:Union{Union,DataType}}, bounds::Vector{Tuple{Real,Real}},
  dataset::T where T<:AbstractString, id::Union{Nothing,Int,AbstractString})

  # Warn of non-standardised data
  if df.names(data) ≠ standardnames
    @warn "Non-standard names and/or order used for data columns. Trying to correct..."
  end

  ### Check column types, for correctly ordered DataFrames
  drev = DataFrame() # init DataFrame for revised data
  # Sample column numbers or return an empty DataFrame with the correct structure
  # for empty input data
  unchecked = try collect(1:length(data[1,:]))
  catch
    d=DataFrame()
    for (n, t) in zip(standardnames, standardtypes)
      d[!,n] = t[]
    end
    return d
  end
   # init vector with column numbers to check
  for i = 1:length(standardnames)
    try checktype(data, standardnames[i], standardtypes[i])
      # Remove column from unchecked list, if column was identified by name and
      # tests passed. Save data to revised data and remove column from unchecked list.
      drev[!,standardnames[i]] = data[!, standardnames[i]]
      unchecked = unchecked[unchecked.≠i]
    catch
      # If column couldn't be identified by name, take the first column not yet checked
      # and run tests.
      try checktype(data, unchecked[1], standardtypes[i])
        # If tests pass, save revised data and remove column from unchecked list.
        drev[!,standardnames[i]] = data[!, unchecked[1]]
        unchecked = unchecked[unchecked.≠i]
      catch e
        # If tests fail throw error for essential time, lat/lon columns
        if i ≤ 3
          rethrow(e)
        elseif !isempty(unchecked)
          # or return column with missing values for non-essential data
          # and issue warnings for wrong data types
          @warn string("Mismatch in column type for column $(standardnames[i]). ",
            "Column filled with `missing` in: "), dataset, id
          drev[!,standardnames[i]] = [missing for j = 1:length(drev[!,1])]
        else
          # or warn of missing columns
          @warn string("Missing data for column $(standardnames[i]). ",
            "Column filled with `missing`: "), dataset, id
          drev[!,standardnames[i]] = [missing for j = 1:length(drev[!,1])]
        end
        # Remove column from unchecked list
        unchecked = unchecked[unchecked.≠i]
      end
    end
  end

  ### Check data bounds in various columns
  # Throw error for lat/lon (essential)
  checkrange(drev.lat, (-90, 90))
  checkrange(drev.lon, (-180, 180))

  # Replace column with missing column for non-essential data
  for (i, r) in enumerate(bounds)
    try checkrange(drev[!,i+3], r)
    catch
      @warn string("Mismatch in data range for column $(standardnames[i+3]) ",
        "in flight $id of $dataset dataset. Column filled with `missing`.")
      drev[!,i+3] = [missing for j = 1:length(drev[!,1])]
    end
  end

  # Warn of additonal columns
  if length(unchecked) > 0
    @warn "More columns than expected in data. Extra columns ignored."
  end

  # Return revised data
  return drev
end #function checkcols


""" Check the type of a `DataFrame` column, throw AssertionError on failure """
checktype(df::DataFrame, col::Union{Symbol,Int}, coltype::Union{Union,DataType}) =
  @assert typeof(df[!,col]) <: coltype "column $col is not of type $coltype"


""" Check the data range of all values in a `DataFrame` column, throw AssertionError on failure """
checkrange(v::Vector, bounds::Tuple{Real,Real}=(-Inf,Inf)) =
  @assert(all(bounds[1] .≤ v[.!ismissing.(v)] .≤ bounds[2]),
    "Vector out of range. Expected [$(bounds[1])...$(bounds[2])], got [$(minimum(v))...$(maximum(v))].")


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
    extract_timespan(data::DataFrame, t::Int, timespan::Int=15) -> DataFrame

From the `data` in a `DataFrame`, extract a subset at index (row) `t` ± `timespan`
(rows).
"""
function extract_timespan(data::DataFrame, t::Int, timespan::Int=15)
  t1 = max(1, min(t-timespan, length(data[:,1])))
  t2 = min(length(data[:,1]), t+timespan)
  t2 = t-timespan > length(data[:,1]) ? 0 : t2
  return data[t1:t2,:]
end #function extract_timespan


"""
    swap_sattype(sattype::Symbol)::Symbol

Returns the other Symbol from the options

- `:CLay`
- `:CPro`

when `sattype` is given to switch between `SatDB` datasets.
"""
function swap_sattype(sattype::Symbol)::Symbol
  swap = Dict(:CPro => :CLay, :CLay => :CPro)
  swap[sattype]
end


"""
    get_trackdata(flight::FlightData, sat::SatDB, sattype::Symbol,
      tflight::DateTime, tsat::DateTime, flightspan::Int, satspan::Int)
      -> flightdata::FlightData, satdata::SatDB

From the measured `flight` and `sat` data and the times of the aircraf (`tflight`)
and the satellite (`tsat`) at the intersection, save the closest measured value
to the interpolated intersection of each dataset ±`flightspan`/`satspan` data points
to `flightdata` and `satdata`, respectively, and return the datasets.
"""
function get_trackdata(flight::FlightData, sat::SatDB, sattype::Symbol,
  tflight::DateTime, tsat::DateTime, flightspan::Int, satspan::Int)

  # Find the index (DataFrame row) of the intersection in the flight data
  tf = argmin(abs.(flight.data.time .- tflight))

  # Construct FlightData at Intersection
  flightdata =
    FlightData(extract_timespan(flight.data, tf, flightspan), flight.metadata)

  # Get satellite data used to find the intersection and find DataFrame row of intersection
  satprim = getfield(sat, sattype).data
  tsp = argmin(abs.(satprim.time .- tsat))
  # Retrieve DataFrame at Intersection ± 15 time steps
  primdata = extract_timespan(satprim, tsp, satspan)

  # Switch to other satellite data (CLay/CPro) and check whether data is available
  # at intersection and retrieve ±15 time steps as well
  sattype = swap_sattype(sattype)
  secdata = try satsec = getfield(sat, sattype).data
    tss = argmin(abs.(satsec.time .- tsat))
    extract_timespan(satsec, tss, satspan)
  catch
    # Return an empty DataFrame, if no data is available
    DataFrame(time=DateTime[], lat=Float64[], lon=Float64[])
  end
  # Save satellite data in SatDB
  satdata = sattype == :CPro ? SatDB(CLay(primdata), CPro(secdata), sat.metadata) :
    SatDB(CLay(secdata), CPro(primdata), sat.metadata)

  # Instatiate new Intersection
  return flightdata, satdata
end


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
