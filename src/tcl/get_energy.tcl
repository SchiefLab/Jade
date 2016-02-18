
mol new ../common/ubq_wb.psf 
mol addfile ubq_wb_eq_1fs.restart.vel type namdbin waitfor all

set all [atomselect top all]
set fil [open energy.dat w]
foreach m [$all get mass] v [$all get {x y z}] {
puts $fil [expr 0.5* $m * [vecdot $v $v]]
}

close $fil
