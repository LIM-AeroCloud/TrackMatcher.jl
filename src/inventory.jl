"""
    loadInventory(files)

documentation
"""
function loadInventory(files)
  # Initialise inventory file array
  inventory = flightData[]
  for file in files
    # Load data
    flights = DataFrame(csv.load(file, skiplines_begin=3,
      header_exists=false, colnames=[:id, :seg, :lat, :lon, :alt, :month, :day, :year,
      :hour, :min, :sec, :emiss, :temp, :press, :rH, :speed, :segtime, :segdist, :thrust,
      :weight, :fuel, :co, :hc, :nox, :pmnv, :pmso, :pmfo, :co2, :h2o, :sox],
      colparsers=Dict(:id => Int32, :lat => Float32, :lon =>Float32, :alt => Float32, :speed => Float32)))
    df.deleterows!(flights, length(flights.id)) #delete additional line at end of file
    # Calculate time from individual columns and add as DateTime to DataFrame
    flights.time = ZonedDateTime.(flights.year, flights.month, flights.day,
      flights.hour, flights.min, flights.sec, tz.tz"UTC")

    # Initialise loop over file
    global iStart = 1
    # Loop over all data points
    for i = 1:length(flights.time)-1
      if flights.id[i] ≠ flights.id[i+1]
        # When flight ID changes save data as flightData and set start index to next dataset
        push!(inventory, flightData(flights.time[iStart:i], flights.lat[iStart:i],
          flights.lon[iStart:i], flights.alt[iStart:i], [missing for j=iStart:i],
          [missing for j=iStart:i], flights.speed[iStart:i]))
        global iStart = i+1
      end
    end
    # Save last flight of the dataset
    i = length(flights.time)
    push!(inventory, flightData(flights.time[iStart:i], flights.lat[iStart:i],
      flights.lon[iStart:i], flights.alt[iStart:i], [missing for j=iStart:i],
      [missing for j=iStart:i], flights.speed[iStart:i]))
  end

  return inventory
end #function loadInventory


"""
    findcsv(folder::String)

Load inventory files saved in `folder` to the `inventory` holding the file names
and locations as a vector of strings.
"""
function findcsv(inventory::Vector{String}, folder::String)
  # Scan directory for files and folders and save directory
  dir = readdir(folder); path = abspath(folder)
  for file in dir
    cwd = joinpath(path, file)
    if endswith(file, ".csv")
      push!(inventory, cwd)
    elseif isdir(cwd)
      inventory = findcsv(inventory, cwd)
    end
  end

  return inventory
end # function findcsv
