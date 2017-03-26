Creates netCDF files of winds specified in the HurDAT data,
and a binary map of hurricane/not hurricane wind locations.

Run make_jobs.tcl to create a Makefile, then

make -jX

where X is the number of parallel threads to run

The output is put in ./output (make sure you create it!)