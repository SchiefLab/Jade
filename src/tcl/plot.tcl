set log		push.log
set var		"Center of mass"

set varcol	[expr 2 + [llength [split $var]]]

package require multiplot

set inStream [open $log r]

set x {}
set y {}
set i 0

foreach line [split [read $inStream] \n] {
    set columns [split $line]
    if { [lindex $columns 0] == "TCL:" } {
	if { [string first $var $line] != -1 } {
	    lappend x $i
	    lappend y [lindex $columns $varcol]
	    incr i
	}
    }
}

set plothandle [multiplot -x $x -y $y -ylabel $var -lines -marker point -plot]
