library(ncdf4)

# Spatial
EARTH_RADIUS <- 6371

# We only focus on the North Atlantic
MIN_LON <- -110
MAX_LON <- 60
MIN_LAT <- -30
MAX_LAT <- 90

# Longitudes and latitudes
# In ERA-Interim, the latitudes start at the north pole
LONS <- seq(-110, 60, 0.125)
ERA_LATS <- seq(90, -30, -0.125)
LATS <- seq(-30, 90, 0.125)

# Pre-defined netCDF dimensions for latitudes and longitudes
LON_DIM <- ncdim_def("longitude", "degrees_east", LONS)
LAT_DIM <- ncdim_def("latitude", "degrees_north", LATS)

# Convert a string lat/lon value to a numeric one
# negative_hemisphere indicates which ending character
# results in a negative number (west or south)
extract_coord <- function(coord_string, negative_hemisphere) {
    last_pos <- nchar(coord_string)

    coord_number <- as.numeric(substring(coord_string, 1, last_pos - 1))
    coord_hemisphere <- substring(coord_string, last_pos, last_pos)
    if (coord_hemisphere == negative_hemisphere) {
        coord_number <- coord_number * -1
    }

    return (coord_number)
}

# Extract a numerical latitude value from a string
extract_lat <- function(lat_string) {
    return (extract_coord(lat_string, "S"))
}

# Extract a numerical longitude value from a string
extract_lon <- function(lon_string) {
    return (extract_coord(lon_string, "W"))
}

# Get the index of a given longitude
get_lon_index <- function(lon) {
    return (which(LONS >= lon)[1])
}

# Get the index of a given latitude
get_lat_index <- function(lat) {
    return (which(LATS >= lat)[1])
}

# Get the position of a point that is a given
# direction and distance from a starting point
# The result is a vector of (lon, lat)
get_position_at_distance <- function(start_lon, start_lat, direction, distance) {

    start_lon_radians <- to_radians(start_lon)
    start_lat_radians <- to_radians(start_lat)
    dir_radians <- direction * (pi / 180)

    end_lat <- asin(sin(start_lat_radians) * cos(distance / EARTH_RADIUS) +
        cos(start_lat_radians) * sin(distance / EARTH_RADIUS) * cos(dir_radians))

    end_lon <- start_lon_radians + atan2(sin(dir_radians) * sin(distance / EARTH_RADIUS) * cos(start_lat_radians),
        cos(distance / EARTH_RADIUS) - sin(start_lat_radians) * sin(end_lat))

    return (c(to_degrees(end_lon), to_degrees(end_lat)))
}

# Get the distance between two points
get_distance_between <- function(lon1, lat1, lon2, lat2) {

    lat1_radians <- to_radians(lat1)
    lat2_radians <- to_radians(lat2)

    delta_lat <- lat2_radians - lat1_radians
    delta_lon <- to_radians(lon2 - lon1)

    a <- sin(delta_lat / 2) * sin(delta_lat / 2) +
         cos(lat1_radians) * cos(lat2_radians) *
         sin(delta_lon / 2) * sin(delta_lon / 2)

    c <- 2 * atan2(sqrt(a), sqrt(1 - a))

    return (EARTH_RADIUS * c)

}

# Convert a degrees value to radians
to_radians <- function(degrees) {
    return (degrees * (pi / 180))
}

# Convert a radians value to degrees
to_degrees <- function(radians) {
    return (radians * (180 / pi))
}

# Temporal

# The number of days in each month of the year
MONTH_DAYS <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# Calculate the netCDF time dimension for a given month
# The dimension contains the number of days * 4 to give 6-hourly values
get_time_dim <- function(month) {
    return(ncdim_def("time", "day of month", seq(1, MONTH_DAYS[month] + 0.75, 0.25)))
}

# Get the index of a given day and time within the
# month. The time must be 0000, 0600, 1200 or 1800.
# Any other values return a negative index indicating
# that it is invalid
get_time_index <- function(day, time) {
    index <- (day - 1) * 4

    if (time == 0) {
        index <- index + 1
    } else if (time == 600) {
        index <- index + 2
    } else if (time == 1200) {
        index <- index + 3
    } else if (time == 1800) {
        index <- index + 4
    } else {
        index <- -1
    }

    return (index)
}

get_surrounding_cells <- function(lon_index, lat_index, step, distance_limit) {

    cell_count <- step * 8

    cells <- vector(mode="numeric", length=cell_count * 2)
    cells[cells == 0] <- NA
    dim(cells) <- c(cell_count,2)

    cell_count <- 0
    for (y in (step * -1):step) {

        cell_y <- get_surround_y_cell(lat_index, y)
        if (!is.na(cell_y)) {

            # For the top and bottom rows of the step grid, add all horizontal cells
            if (abs(y) == step) {
                for (x in (step * -1):step) {
                    cell_count <- cell_count + 1
                    cells[cell_count,1] <- get_surround_x_cell(lon_index, x)
                    cells[cell_count,2] <- cell_y
                }
            } else {
                # For all other rows, just add the left and right edges
                cell_count <- cell_count + 1
                cells[cell_count,1] <- get_surround_x_cell(lon_index, (step * -1))
                cells[cell_count,2] <- cell_y

                cell_count <- cell_count + 1
                cells[cell_count,1] <- get_surround_x_cell(lon_index, step)
                cells[cell_count,2] <- cell_y
            }
        }
    }

    return (cells)
}

get_surround_x_cell <- function(lon, x) {
    new_x <- lon + x
    if (new_x > length(LONS) || new_x < 1) {
        new_x <- NA
    }
    return (new_x)
}

get_surround_y_cell <- function(lat, y) {
    new_y <- lat + y
    if (new_y > length(LATS) || new_y < 1) {
        new_y <- NA
    }
    return (new_y)
}
