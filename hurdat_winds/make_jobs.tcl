#!/usr/bin/tclsh


set jobList ""
set jobText ""

for {set year 2004} {$year <= 2015} {incr year} {
    for {set month 1} {$month <= 12} {incr month} {
        append jobList "${year}_${month} "

        append jobText "${year}_${month} :\n"
        append jobText "\tR --no-save --slave year=\"${year}\" month=\"${month}\" < make_hurdat_winds.R\n\n"
    }
}

set outChan [open "Makefile" w]
puts $outChan "all : $jobList\n"
puts $outChan $jobText
close $outChan
