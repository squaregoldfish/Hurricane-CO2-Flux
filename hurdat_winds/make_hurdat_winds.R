library(ncdf4)
source("../common/coords.R")
source("../common/conversion.R")

HURDAT_FILE <- "../data/HurDAT2/hurdat2.txt"
OUTPUT_DIR <- "output"

write_netcdf <- function(year, month, winds) {

    mask <- winds
    mask[!is.na(mask)] <- 1
    mask[is.na(mask)] <- 0

    if (month < 10) {
        month_string <- paste("0", month, sep="")
    } else {
        month_string <- month
    }
    nc_file <- paste(OUTPUT_DIR, "/hurdat_wind-", year, "_", month_string, ".nc", sep="")
    wind_var <- ncvar_def("hurdat_wind", "ms-1", list(LON_DIM, LAT_DIM, get_time_dim(as.numeric(month))), -9999)
    nc <- nc_create(nc_file, list(wind_var))
    ncvar_put(nc, wind_var, winds)
    nc_close(nc)
}

add_winds <- function(winds, time_index, start_lon, start_lat, lon_direction, lat_direction, radius, speed) {

    lon_start_index <- get_lon_index(start_lon)
    if (lon_direction == "e") {
        lon_degrees <- 90
    } else {
        lon_degrees <- 270
    }
    
    end_lon <- get_position_at_distance(start_lon, start_lat, lon_degrees, radius)[1]
    lon_end_index <- get_lon_index(end_lon)

    lat_start_index <- get_lat_index(start_lat)
    if (lat_direction == "n") {
        lat_degrees <- 0
    } else {
        lat_degrees <- 180
    }

    end_lat <- get_position_at_distance(start_lon, start_lat, lat_degrees, radius)[2]
    lat_end_index <- get_lat_index(end_lat)

    for (lon_loop in lon_start_index:lon_end_index) {
        real_lon <- LONS[lon_loop]

        for (lat_loop in lat_start_index:lat_end_index) {

            real_lat <- LATS[lat_loop]

            if (get_distance_between(start_lon, start_lat, real_lon, real_lat) <= radius) {
                if (is.na(winds[lon_loop, lat_loop, time_index]) || winds[lon_loop, lat_loop, time_index] < speed) {
                    winds[lon_loop, lat_loop, time_index] <- speed
                }
            }
        }
    }

    return (winds)
}

#################################################################

# Command line parameters
for (arg in commandArgs()) {
    argset <- strsplit(arg, "=", fixed=TRUE)
    if (!is.na(argset[[1]][2])) {
        if (argset[[1]][1] == "year") {
            assign("year",argset[[1]][2])
        } else if (argset[[1]][1] == "month") {
            assign("month",argset[[1]][2])
        }
    }
}

year <- as.numeric(year)
month <- as.numeric(month)

month_winds <- vector(mode="numeric", length=(length(LONS) * length(LATS) * (MONTH_DAYS[as.numeric(month)] * 4)))
dim(month_winds) <- c(length(LONS), length(LATS), (MONTH_DAYS[as.numeric(month)] * 4))
month_winds[month_winds == 0] <- NA
wind_found <- FALSE

row_count <- 0
con <- file(HURDAT_FILE, open="r")
while (length(line <- readLines(con, n=1, warn=F)) > 0) {
    row_count <- row_count + 1

    if (substring(line, 1, 2) != "AL") {
        row <- strsplit(gsub(" ", "", line), ",")[[1]]
        row_date <- row[1]
        row_year <- as.numeric(substring(row_date, 1, 4))
        row_month <- as.numeric(substring(row_date, 5, 6))

        if (row_year == year && row_month == month) {
            cat("\r", row_count)
            row_day <- as.numeric(substring(row_date, 7, 8))
            time <- as.numeric(row[2])
            time_index <- get_time_index(row_day, time)

            centre_lon <- extract_lon(row[6])
            centre_lat <- extract_lat(row[5])
            centre_wind <- knots_to_ms(row[7])

            month_winds[get_lon_index(centre_lon), get_lat_index(centre_lat), time_index] <- centre_wind
            wind_found <- TRUE

            ms64 <- knots_to_ms(64)
            k64ne_radius <- nm_to_km(row[17])
            if (k64ne_radius > 0) {
                month_winds <- add_winds(month_winds, time_index, centre_lon, centre_lat, "e", "n", k64ne_radius, ms64)
            }
            
            k64se_radius <- nm_to_km(row[18])
            if (k64se_radius > 0) {
                month_winds <- add_winds(month_winds, time_index, centre_lon, centre_lat, "e", "s", k64se_radius, ms64)
            }
            
            k64sw_radius <- nm_to_km(row[19])
            if (k64sw_radius > 0) {
                month_winds <- add_winds(month_winds, time_index, centre_lon, centre_lat, "w", "s", k64sw_radius, ms64)
            }
            
            k64nw_radius <- nm_to_km(row[20])
            if (k64nw_radius > 0) {
                month_winds <- add_winds(month_winds, time_index, centre_lon, centre_lat, "w", "n", k64nw_radius, ms64)
            }
            
            ms50 <- knots_to_ms(50)
            k50ne_radius <- nm_to_km(row[13])
            if (k50ne_radius > 0) {
                month_winds <- add_winds(month_winds, time_index, centre_lon, centre_lat, "e", "n", k50ne_radius, ms50)
            }
            
            k50se_radius <- nm_to_km(row[14])
            if (k50se_radius > 0) {
                month_winds <- add_winds(month_winds, time_index, centre_lon, centre_lat, "e", "s", k50se_radius, ms50)
            }
            
            k50sw_radius <- nm_to_km(row[15])
            if (k50sw_radius > 0) {
                month_winds <- add_winds(month_winds, time_index, centre_lon, centre_lat, "w", "s", k50sw_radius, ms50)
            }
            
            k50nw_radius <- nm_to_km(row[16])
            if (k50nw_radius > 0) {
                month_winds <- add_winds(month_winds, time_index, centre_lon, centre_lat, "w", "n", k50nw_radius, ms50)
            }
            
            ms35 <- knots_to_ms(35)
            k35ne_radius <- nm_to_km(row[9])
            if (k35ne_radius > 0) {
                month_winds <- add_winds(month_winds, time_index, centre_lon, centre_lat, "e", "n", k35ne_radius, ms35)
            }
            
            k35se_radius <- nm_to_km(row[10])
            if (k35se_radius > 0) {
                month_winds <- add_winds(month_winds, time_index, centre_lon, centre_lat, "e", "s", k35se_radius, ms35)
            }
            
            k35sw_radius <- nm_to_km(row[11])
            if (k35sw_radius > 0) {
                month_winds <- add_winds(month_winds, time_index, centre_lon, centre_lat, "w", "s", k35sw_radius, ms35)
            }
            
            k35nw_radius <- nm_to_km(row[12])
            if (k64ne_radius > 0) {
                month_winds <- add_winds(month_winds, time_index, centre_lon, centre_lat, "w", "n", k35nw_radius, ms35)
            }
        }
    }
}

if (wind_found) {
    write_netcdf(year, month, month_winds)
}

cat("\n")
