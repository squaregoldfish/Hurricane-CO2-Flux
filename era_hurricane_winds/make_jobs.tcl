#!/usr/bin/tclsh

set jobList ""
set jobText ""

set wind_file_list [glob -nocomplain "ERA-Interim/wind*nc"]

foreach wind $wind_file_list {
	regexp "wind_(\[0-9\]*)-(\[0-9\]*)" $wind dummy year month
    append jobList "${year}_${month} "

    append jobText "${year}_${month} :\n"
    append jobText "\tR --no-save --slave year=\"${year}\" month=\"${month}\" < make_era_hurricane_winds.R\n\n"
}

set outChan [open "Makefile" w]
puts $outChan "all : $jobList\n"
puts $outChan $jobText
close $outChan
