for {set i 0} {$i<=$nf} {incr i 10} {[atomselect top all frame $i] writepdb BeeLH_1ns_${i}ps.pdb}
