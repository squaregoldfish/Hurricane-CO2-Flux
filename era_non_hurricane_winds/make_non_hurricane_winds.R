library(ncdf4)

INPUT_DIR <- "ERA-Interim"
HURRICANE_DIR <- "hurricane_winds"
OUTPUT_DIR <- "output"
YEAR_RANGE <- seq(2004, 2015)


lons <- NULL
lats <- NULL

for (month in 1:12) {

	month_string <- month
	if (month < 10) {
		month_string <- paste("0", month, sep="")
	}

	time_steps <- 0
	wind_totals <- NULL
	wind_counts <- NULL

	for (year in YEAR_RANGE) {
		wind_file <- paste(INPUT_DIR, "/wind_", year, "-", month_string, ".nc", sep="")

		if (file.exists(wind_file)) {
			cat("\r", month, year, "   ")
			
			nc <- nc_open(wind_file)
			if (is.null(lats)) {
				lons <- ncvar_get(nc, "longitude")
				lats <- ncvar_get(nc, "latitude")
			}

			
			in_wind <- ncvar_get(nc, "wind")

			if (is.null(wind_totals)) {
				time_steps <- dim(in_wind)[3]
				wind_totals <- vector(mode="numeric", length=(length(lons) * length(lats) * time_steps))
				dim(wind_totals) <- c(length(lons), length(lats), time_steps)
			}


			hurricane_file <- paste(HURRICANE_DIR, "/era_hurricane_wind-", year, "_", month_string, ".nc", sep="")
			print(hurricane_file)
			

			quit("no")
		}
	}

}

cat("\n")
