#
# Generate monthly netCDF files of wind speeds in the HurDAT2 data set.
# Files are 0.125 degree resoltion and 6 hourly time steps
#

# Imports
library(ncdf4)
source("../common/coords.R")
source("../common/conversion.R")

# Paths
HURDAT_FILE <- "../data/HurDAT2/hurdat2.txt"
ERA_DIR <- "ERA-Interim"
OUTPUT_DIR <- "output"

# Write the netCDF file for a given month
write_netcdf <- function(year, month, winds) {

    mask <- winds
    mask[!is.na(mask)] <- 1
    mask[is.na(mask)] <- 0

    if (month < 10) {
        month_string <- paste("0", month, sep="")
    } else {
        month_string <- month
    }
    nc_file <- paste(OUTPUT_DIR, "/era_hurricane_wind-", year, "_", month_string, ".nc", sep="")
    wind_var <- ncvar_def("era_hurricane_wind", "ms-1", list(LON_DIM, LAT_DIM, get_time_dim(as.numeric(month))), -9999)
    nc <- nc_create(nc_file, list(wind_var))
    ncvar_put(nc, wind_var, winds)
    nc_close(nc)
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


# Open the ERA winds file
if (month < 10) {
    era_file <- paste(ERA_DIR, "/wind_", year, "-0", month, ".nc", sep="")
} else {
    era_file <- paste(ERA_DIR, "/wind_", year, "-", month, ".nc", sep="")
}

nc <- nc_open(era_file)
lons <- ncvar_get(nc, "longitude")
lats <- ncvar_get(nc, "latitude")
times <- ncvar_get(nc, "time")
era_wind <- ncvar_get(nc, "wind")
nc_close(nc)

# Set up the output data set
hurricane_winds <- vector(mode="numeric", length=(length(lons) * length(lats) * length(times)))
dim(hurricane_winds) <- c(length(lons), length(lats), length(times))
hurricane_winds[hurricane_winds == 0] <- NA
wind_found <- FALSE

# Open the HurDAT file and read one line at a time
row_count <- 0
con <- file(HURDAT_FILE, open="r")
while (length(line <- readLines(con, n=1, warn=F)) > 0) {
    row_count <- row_count + 1

    # Only use actual data lines
    if (substring(line, 1, 2) != "AL") {

        # Extract the year and month to see if
        # we need to process this line
        row <- strsplit(gsub(" ", "", line), ",")[[1]]
        row_date <- row[1]
        row_year <- as.numeric(substring(row_date, 1, 4))
        row_month <- as.numeric(substring(row_date, 5, 6))

        if (row_year == year && row_month == month) {
            
            # Calculate the time index
            row_day <- as.numeric(substring(row_date, 7, 8))
            time <- as.numeric(row[2])
            time_index <- get_time_index(row_day, time)

            if (time_index != -1) {

                # Get the central point and wind speed
                centre_lon <- extract_lon(row[6])
                centre_lat <- extract_lat(row[5])

                ne_radius <- nm_to_km(row[9])
                se_radius <- nm_to_km(row[10])
                sw_radius <- nm_to_km(row[11])
                nw_radius <- nm_to_km(row[12])


                max_radius <- ne_radius
                if (se_radius > max_radius) {
                    max_radius <- se_radius
                }
                if (sw_radius > max_radius) {
                    max_radius <- sw_radius
                }
                if (nw_radius > max_radius) {
                    max_radius <- nw_radius
                }

                print(max_radius)
            }
        }
    }
}

# Write the netCDF file
if (wind_found) {
    write_netcdf(year, month, month_winds)
}

cat("\n")
