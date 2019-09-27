"""
    loadInventory(files)

documentation
"""
function loadInventory(files)
  # Initialise inventory file array
  inventory = flightData[]
  for file in files
    # Load data
    flights = CSV.read(file, datarow=3, copycols=true)
    # Delete additional line at end of file
    df.deleterows!(flights, length(flights.FLIGHT_ID)-1:length(flights.FLIGHT_ID))
    # Calculate time from individual columns and add as DateTime to DataFrame
    flights.time = ZonedDateTime.(flights.SEGMENT_YEAR, flights.SEGMENT_MONTH,
      flights.SEGMENT_DAY, flights.SEGMENT_HOUR, flights.SEGMENT_MIN,
      flights.SEGMENT_SEC, tz.tz"UTC")

    # Initialise loop over file
    global iStart = 1
    # Loop over all data points
    for i = 1:length(flights.time)-1
      if flights.FLIGHT_ID[i] ≠ flights.FLIGHT_ID[i+1]
        # When flight ID changes save data as flightData and set start index to next dataset
        push!(inventory, flightData(flights.time[iStart:i], flights.LATITUDE[iStart:i],
          flights.LONGITUDE[iStart:i], flights.ALTITUDE[iStart:i], [missing for j=iStart:i],
          [missing for j=iStart:i], flights.SPEED[iStart:i]))
        global iStart = i+1
      end
    end
    # Save last flight of the dataset
    i = length(flights.time)
    push!(inventory, flightData(flights.time[iStart:i], flights.LATITUDE[iStart:i],
      flights.LONGITUDE[iStart:i], flights.ALTITUDE[iStart:i], [missing for j=iStart:i],
      [missing for j=iStart:i], flights.SPEED[iStart:i]))
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
