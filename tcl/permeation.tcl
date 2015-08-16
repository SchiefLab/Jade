# This script identifies water molecules that permeate through the
# nanotube layer

# Before running it, first load the dcd file in VMD, and make sure
# it's the "top" molecule.

# Specify the upper and lower boundaries of the nanotube layer
set upperEnd 6.7
set lowerEnd -6.7

# Start counting permeation events after this number of frames
set skipFrame 500


puts "Computing permeation events... (please wait)"

set wat [atomselect top "name OH2"]
set segList [$wat get segname]
set ridList [$wat get resid]
set labelList {}
foreach foo $segList {
  lappend labelList 0
}

set num1 0
set num2 0
set numFrame [molinfo top get numframes]


# Each water molecule has a label, which has 5 possible values
#  2: Above the nanotube layer
# -2: Below the nanotube layer
#  1: Inside the nanotube layer, entering from upper surface
# -1: Inside the nanotube layer, entering from lower surface
#  0: Inside the nanotube layer from the beginning

# For every frame, the label of each water molecule is
# determined, and compared with its label in the previous frame.
# If the new label is +2 (or -2), while the old label is -1 (or +1),
# it means the water molecule has traversed the nanotube, thus a
# permeation event is reported and counted. If a water molecule
# is inside the nanotube layer in the current frame, its label
# will be determined by its old label.

for {set fr 0} {$fr < $numFrame} {incr fr} {
  molinfo top set frame $fr
  set oldList $labelList
  set labelList {}
  foreach z [$wat get z] oldLab $oldList segname $segList resid $ridList {
    if {$z > $upperEnd} {
      set newLab 2
      if {$oldLab == -1} {
        puts "$segname:$resid permeated through the nanotubes along +z direction at frame $fr"
        if {$fr >= $skipFrame} {
          incr num1
        }
      }
    } elseif {$z < $lowerEnd} {
      set newLab -2
      if {$oldLab == 1} {
        puts "$segname:$resid permeated through the nanotubes along -z direction at frame $fr"
        if {$fr >= $skipFrame} {
          incr num2
        }
      }
    } elseif {abs($oldLab) > 1} {
      set newLab [expr $oldLab / 2]
    } else {
      set newLab $oldLab
    }
    lappend labelList $newLab
  }
}

puts ""
set nf [expr $numFrame-$skipFrame]
if {$nf >= 0} {
  puts "The total number of permeation events during $nf frames in +z direction is: $num1"
  puts "The total number of permeation events during $nf frames in -z direction is: $num2"
} else {
  puts "The specified first frame ($skipFrame) is larger than the total number of frames ($numFrame)"
}
