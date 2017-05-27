library(ncdf4)

#in_files <- list.files(pattern="uvwind*")
in_files <- c("uvwind_2005-09.nc", "uvwind_2008-09.nc", "uvwind_2010-09.nc", "uvwind_2012-09.nc", "uvwind_2014-09.nc", "uvwind_2015-09.nc")

for (file_loop in 1:length(in_files)) {

	in_file <- in_files[file_loop]
	print(in_file)
	out_file <- paste("wind_", substring(in_file, 8, 14), ".nc", sep="")

	nc <- nc_open(in_file)
	lons <- ncvar_get(nc, "longitude")
	lats <- ncvar_get(nc, "latitude")
	times <- ncvar_get(nc, "time")

	u10 <- ncvar_get(nc, "u10")
	v10 <- ncvar_get(nc, "v10")
	nc_close(nc)

	wind <- sqrt(u10^2 + v10^2)

	lon_dim <- ncdim_def("longitude", "degrees_east", lons)
	lat_dim <- ncdim_def("latitude", "degrees_north", lats)
	time_dim <- ncdim_def("time", "hours since 1900-01-01 00:00:0.0", times)

	wind_var <- ncvar_def("wind", "m s**-1", list(lon_dim, lat_dim, time_dim), prec="float", compression=9)

	nc <- nc_create(out_file, list(wind_var))
	ncvar_put(nc, wind_var, wind)
	nc_close(nc)
}