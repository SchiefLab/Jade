set mut_fraction 0.3
if {$mut_fraction <= 0.0 || $mut_fraction > 1.0} {
 puts "Fraction of residues for mutation (mut_fraction)"
 puts "should satisfy 0 < mut_fraction < 1;"
 puts "exiting."
 return
}

set head [atomselect top "name PCH"]
set nl [$head num]
$head delete


set ns [expr int($nl*$mut_fraction)]
set selected {}
for {set i 0} {$i<$nl} {incr i} {
   set a($i) [expr $i +1]
}
for {set k 1} {$k<=$ns} {incr k} {
   set rang [expr $nl-$k+1]
   set n [expr {int(rand()*$rang)}]
   lappend selected $a($n)
   for {set j 0} {$j< [expr $nl - $k]} {incr j} {
      if {$j<$n} {
         set a($j) $a($j)
      } else {
         set jj [expr $j+1]
         set a($j) $a($jj)
      }
   }
}
set fid [open "selectedresidues.txt" w+]
puts $fid $selected
close $fid


for {set i 0} {$i<$ns} {incr i} {
   set resi [lindex $selected $i]
   set ps [atomselect top "residue $resi and name PCH"]
   $ps set name "PSH"
   set dops [atomselect top "residue $resi"]
   $dops set resname "DOPS"
   $dops delete
   $ps delete
}
set all [atomselect top all]
$all writepdb mixture.pdb
