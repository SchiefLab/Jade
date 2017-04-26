set outfile myenergy.dat
set energy Total

set file [open $outfile r]
while { [gets $file line] != -1 } {
  lappend count $line
  }
set lc [llength $count]
unset count
close $file

set file [open $outfile r]
set avg_count 0
set avg_square 0

for {set i 1} {$i < $lc} {incr i} {
  
  if {$i == 1} {
    set one [gets $file]
    for {set j 0} {$j < [llength $one]} {incr j} {
      set lab [lindex $one $j]
      if {$lab == $energy} {
        set index $j
        }
      }
    }
    
  set ene [lindex [gets $file] $index]
  set avg_count [expr $avg_count + $ene]
  set avg_square [expr ($avg_square + $ene*$ene)]
  }
  
puts "Average $energy: [expr $avg_count/($lc-1)]"
puts "Squared Average $energy: [expr $avg_square/($lc-1)]"
close $file
