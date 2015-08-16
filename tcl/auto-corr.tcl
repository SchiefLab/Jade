
set file [open TEMP.dat r]
while { [gets $file line] != -1 } {
  lappend templist $line
  }
close $file
set lc [llength $templist]
unset templist

set data " "
set file [open TEMP.dat r]
while { [gets $file line] != -1 } {
  set data "$data $line"
  }
close $file

set endlag 25

proc avg {dataarray} {
global $dataarray
set numdata [array size $dataarray]
set total 0
for {set k 0} {$k < $numdata} {incr k} {
  set total [expr "$[list $dataarray]($k)" + $total]
  }
return [expr $total/$numdata]
}

for {set k 0} {$k < $lc} {incr k} {
#  set timestep($k) [lindex $data [expr 2*$k]]
  set temper($k) [lindex $data [expr 2*$k+1]]
  }
unset data

set avg_temp [avg temper]

for {set k 0} {$k < $lc} {incr k} {
  set temper_adj($k) [expr $temper($k) - $avg_temp]
  }

set file [open auto-corr.dat w]

for {set lag 0} {$lag <= $endlag} {incr lag} {
  for {set k 0} {$k < [expr $lc-$lag]} {incr k} {
    set data1($k) $temper_adj($k)
    set data2($k) $temper_adj([expr $k+$lag])
    set data2sq($k) [expr $data2($k) * $data2($k)]
    set dataprod($k) [expr $data1($k) * $data2($k)]
    }
  set mean1 [avg data1]
  set mean2 [avg data2sq]
  set meanprod [avg dataprod]
  
# Calculate the Autocorrelation Function
  puts $file "$lag\t [expr ($meanprod - $mean1*$mean1)/($mean2 - $mean1*$mean1)]"
  }
  
close $file
