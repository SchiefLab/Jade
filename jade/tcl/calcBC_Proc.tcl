proc calcBC {min max} {
	set mins [split $min]
	set maxs [split $max]
	set result {}
	lappend result [expr ([lindex $maxs 0] - [lindex $mins 0])]
	lappend result [expr ([lindex $maxs 1] - [lindex $mins 1])]
	lappend result [expr ([lindex $maxs 2] - [lindex $mins 2])]
	return $result
}
