library(ncdf4)

# Spatial
EARTH_RADIUS <- 6371

MIN_LON <- -110
MAX_LON <- 60
MIN_LAT <- -30
MAX_LAT <- 90

LONS <- seq(-110, 60, 0.125)
ERA_LATS <- seq(90, -30, -0.125)
LATS <- seq(-30, 90, 0.125)

LON_DIM <- ncdim_def("longitude", "degrees_east", LONS)
LAT_DIM <- ncdim_def("latitude", "degrees_north", LATS)

extract_coord <- function(coord_string, negative_hemisphere) {
    last_pos <- nchar(coord_string)

    coord_number <- as.numeric(substring(coord_string, 1, last_pos - 1))
    coord_hemisphere <- substring(coord_string, last_pos, last_pos)
    if (coord_hemisphere == negative_hemisphere) {
        coord_number <- coord_number * -1
    }

    return (coord_number)
}

extract_lat <- function(lat_string) {
    return (extract_coord(lat_string, "S"))
}

extract_lon <- function(lon_string) {
    return (extract_coord(lon_string, "W"))
}

get_lon_index <- function(lon) {
    return (which(LONS >= lon)[1])
}

get_lat_index <- function(lat) {
    return (which(LATS >= lat)[1])
}

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

to_radians <- function(degrees) {
    return (degrees * (pi / 180))
}

to_degrees <- function(radians) {
    return (radians * (180 / pi))
}

# Temporal
MONTH_DAYS <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

get_time_dim <- function(month) {
    return(ncdim_def("time", "day of month", seq(1, MONTH_DAYS[month] + 0.75, 0.25)))
}

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
