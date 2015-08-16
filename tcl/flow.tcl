# This script counts the total net water flow through the
# nanutube layer in the trajectory

# Before running it, first load the dcd file in VMD, and make sure
# it's the "top" molecule.

# Specify the upper and lower boundaries of the nanotube layer
set upperEnd 6.7
set lowerEnd -6.7


# The following function sets the status for each water molecule
# status 0: Inside the nanotube layer
# status 1: Above the nanotube layer
# status -1: Below the nanotube layer
proc set_status {} {
  global wat statusList upperEnd lowerEnd
  set statusList {}
  foreach z [$wat get z] {
    if {$z < $lowerEnd} {
      lappend statusList -1
    } elseif {$z > $upperEnd} {
      lappend statusList 1
    } else {
      lappend statusList 0
    }
  }
}


set wat [atomselect top "name OH2"]
set numFrame [molinfo top get numframes]
set total 0

molinfo top set frame 0
set_status

# For every frame, the status of each water molecule is
# calculated, and compared with its status in the previous frame.
# If the status changes from 0 to +-1 or vice versa, it means that
# this water molecule has crossed one of the boundaries of the
# nanotube layer. The variable "total" records the total number
# of such crossing events (for each event, either +1 or -1 is
# added to "total", according to the crossing direction). However,
# due to the periodic boundary condition, a change of the status
# from +1 to -1 or vice versa doesn't mean the water molecule has
# crossed the channel.

for {set fr 1} {$fr < $numFrame} {incr fr} {
  molinfo top set frame $fr
  set oldList $statusList
  set_status
  foreach oldSt $oldList newSt $statusList {
    if {$oldSt!=$newSt && $oldSt+$newSt!=0} {
      incr total [expr $newSt - $oldSt]
    }
  }
}


# The net flow is taken as the average of the numbers of the
# crossing events for the two boundaries, i.e., one half of "total"

if {$total > 0} {
  puts "The net flow is [expr $total/2.0] water molecules along +z"
} elseif {$total < 0} {
  puts "The net flow is [expr -$total/2.0] water molecules along -z"
} else {
  puts "The net flow is 0"
}
