
puts "Usage: data_avg <logfile> \[<first timestep> <last timestep>]"
puts "   <first timestep> and <last timestep> may be entered as numbers or"
puts "   <first timestep> = 'first' will start at the beginning of the simulation"
puts "   <last timestep> = 'last' will go to the end of the simulation"
puts "Usage: data_time <data stream> <logfile> \[<first timestep> <last timestep>]"
puts "   <data stream> = BOND, ANGLE, DIHED, IMPRP, ELECT, VDW, BOUNDARY, MISC, KINETIC, TOTAL, TEMP, TOTAL2, TOTAL3, TEMPAVG"


proc data_avg {logfile {first 0} {last -1}} {

set file [open $logfile r]
while { [gets $file line] != -1 } {
  if [regexp "ETITLE:" $line] {set etitles $line}
  }
close $file
puts "Calculating averages..."

set file [open $logfile r]
set ener " "
while { [gets $file line] != -1 } {
  if [regexp "ENERGY:  " $line] {set ener "$ener $line"}
  }
close $file

set l [llength $etitles]
set lc [expr [llength $ener]/$l]

for {set k 0} {$k < $l} {incr k} {
  set tstest [lindex $etitles $k]
  if {$tstest == "TS"} {set tsnum $k}
  }

for {set k 0} {$k < $lc} {incr k} {
  set ts([expr $k+1]) [lindex $ener [expr $l*$k+$tsnum]]
  }
set lastts $ts($k)
set lastline $lc
set firstts $ts(1)
set firstline 1

if {$first != "first"} {
  set k 1
  while {$first > $ts($k)} {
    set firstts $ts([expr $k+1])
    set firstline [expr $k+1]
    incr k
    }
  }

if {$last != "last"} { 
  set k 1
  while {$last >= $ts($k)} {
    set lastts $ts($k)
    set lastline $k
    incr k
    }
  }
  
puts "CALCULATING DATA FROM TIMESTEP $firstts TO $lastts:"

set k 0
for {set i 1} {$i <= $lc} {incr i} {
  for {set j 1} {$j <= [llength $etitles]} {incr j} {
    set e($i,$j) [lindex $ener $k]
    incr k
    }
  }

for {set j 3} {$j <= [llength $etitles]} {incr j} {
  set total 0
  for {set i $firstline} {$i <= $lastline} {incr i} {
    set total [expr $total + $e($i,$j)]
    }
  set avg [expr $total/(1+$lastline-$firstline)]
  puts "[lindex $etitles [expr $j-1]]: $avg"
  }

unset ener
unset ts
}


proc data_time {data logfile {first 0} {last -1}} {

set file [open $logfile r]
while { [gets $file line] != -1 } {
  if [regexp "ETITLE:" $line] {set etitles $line}
  }
close $file
puts "Getting $data data..."

set file [open $logfile r]
set ener " "
while { [gets $file line] != -1 } {
  if [regexp "ENERGY:  " $line] {set ener "$ener $line"}
  }
close $file

set lc [expr [llength $ener]/[llength $etitles]]

set k 0
for {set i 1} {$i <= $lc} {incr i} {
  for {set j 1} {$j <= [llength $etitles]} {incr j} {
    set e($i,$j) [lindex $ener $k]
    incr k
    }
  }

set lastts $e($lc,2)
set lastline $lc
set firstts $e(1,2)
set firstline 1

if {$first != "first"} {
  set k 1
  while {$first > $e($k,2)} {
    set firstts $e([expr $k+1],2)
    set firstline [expr $k+1]
    incr k
    }
  }

if {$last != "last"} { 
  set k 1
  while {$last >= $e($k,2)} {
    set lastts $e($k,2)
    set lastline $k
    incr k
    }
  }
  
set j 0

switch $data {
  BOND {set j 3}
  ANGLE {set j 4}
  DIHED {set j 5}
  IMPRP {set j 6}
  ELECT {set j 7}
  VDW {set j 8}
  BOUNDARY {set j 9}
  MISC {set j 10}
  KINETIC {set j 11}
  TOTAL {set j 12}
  TEMP {set j 13}
  TOTAL2 {set j 14}
  TOTAL3 {set j 15}
  TEMPAVG {set j 16}
  default {
    puts "Invalid data stream selection"
    }
  }
if {$j != 0} {
  set file [open $data.dat w]
  for {set i $firstline} {$i <= $lastline} {incr i} {
    puts $file "$e($i,2)	$e($i,$j)"
    }
  close $file
  }
unset ener
puts "Done"
}
