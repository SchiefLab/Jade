# Nx and Ny are the number of "molecules" in x and y directions, respectively;
# dx and dy are the spacings between "molecules" in x and y directions;
# each "molecule" is two beads connected by a spring;
# hname is the atom name for the head bead;
# tname is the atom name for the tail bead;
# rname is the residue name;
# the procedure generates Nx*Ny "molecules" in the upper leaflet
# and Nx*Ny "molecules" in the lower leaflet.
# outPDB is the name for the output pdb file.
proc generate_patch { Nx Ny dx dy hname tname rname outPDB } {
   set output [open $outPDB "w"]
   
   set resid 0
   set atomid 0
   # Upper leaflet
   for {set kx 1} {$kx <= $Nx} {incr kx} {
      set x [expr ($kx - 1)*$dx]
      for {set ky 1} {$ky <= $Ny} {incr ky} {
         set y [expr ($ky - 1)*$dy]
         incr resid
         incr atomid
         puts $output [format \
         "ATOM  %5d %4s %4s %4d   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s" \
         $atomid $tname $rname $resid $x $y 6.25 1.00 0.0 "L" ]
         incr atomid
         puts $output [format \
         "ATOM  %5d %4s %4s %4d   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s" \
         $atomid $hname $rname $resid $x $y 18.75 1.00 0.0 "L" ]
      }
   }
   # Lower leaflet
   for {set kx 1} {$kx <= $Nx} {incr kx} {
      set x [expr ($kx - 1)*$dx]
      for {set ky 1} {$ky <= $Ny} {incr ky} {
         set y [expr ($ky - 1)*$dy]
         incr resid
         incr atomid
         puts $output [format \
         "ATOM  %5d %4s %4s %4d   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s" \
         $atomid $tname $rname $resid $x $y -6.25 1.00 0.0 "L" ]
         incr atomid
         puts $output [format \
         "ATOM  %5d %4s %4s %4d   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s" \
         $atomid $hname $rname $resid $x $y -18.75 1.00 0.0 "L" ]
      }
   }
   
   puts $output "END"

   close $output
}