#!/usr/local/bin/wish


#----------------------TK----------------------

tk::text .result
entry .min -width 25 -relief sunken -textvariable min
.result insert 1.0 "Boundary Conditions\n"
button .calculate -text "Calculate BC" \

grid .min .calculate .max -padx 1m -pady 1m

#---------------------Output------------------

.result insert end [format "cellBasisVector1\t%.2f 0 0\n" [lindex $result 0]]
