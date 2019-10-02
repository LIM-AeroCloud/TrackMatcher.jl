### Select individual flights from Database files
# Find all flight IDs
id = unique(flights.id)
# Initialise
global iStart = 1 # Start index of row range in of current flight in database file
FlightDB = []; # Database with individual flights
# Loop over all flights
for i in id
  # Find end index of current flight in database file
  iEnd = findlast(flights.id.==i)
  println(i, ", ", iStart, ", ", iEnd)
  # Add current flight to database
  push!(FlightDB, FlightData(flights.time[iStart:iEnd], flights.lat[iStart:iEnd],
    flights.lon[iStart:iEnd], flights.alt[iStart:iEnd], [missing for j=iStart:iEnd],
    [missing for j=iStart:iEnd], flights.speed[iStart:iEnd]))
  # Set start index in database file to next flight
  global iStart = iEnd + 1
end
