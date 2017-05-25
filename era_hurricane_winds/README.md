This script extracts hurricanes from the ERA wind fields.
It takes the hurricane radii from the HurDAT data and uses that
to select winds from the ERA files.

The script assumes that hurricanes are circular, as the HurDAT data
frequently only registers winds in one quadrant.

Run make_jobs.tcl to create a Makefile, then

make -jX

where X is the number of parallel threads to run

The output is put in ./output (make sure you create it!)

