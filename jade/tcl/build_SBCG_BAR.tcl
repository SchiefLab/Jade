# Load the two monomers of a BAR domain dimer and make two selections containing the two monomers respectively.
mol load pdb cg_monomer.pdb
set p1 [atomselect top all]
mol load pdb cg_monomer-2.pdb
set p2 [atomselect top all]
# Move the two dimers to the location so that the BAR domain lattice starts on an appropriate place on the membrane. The vector {-130.0 80.0 0.0} is used to position the lattice start apporpriately; the choice of the vector depends on the location of the membrane.  
$p1 moveby {-130.0 80.0 0.0}
$p2 moveby {-130.0 80.0 0.0}


# Move the two monomers by a vector that can result in a lattice shown in Figure 10, and record the locations of the two dimers as separate pdb files. Please note that two rows, each containing three BAR domain dimers (six BAR domain dimers in total), will be generated. 
set k 1
for {set irow 0} {$irow < 2} {incr irow} {
  for {set icopy 0} {$icopy < 3} {incr icopy} {
    $p1 writepdb p1-$k.pdb
    $p2 writepdb p2-$k.pdb
    incr k
    $p1 moveby {125.0 0.0 0.0}
    $p2 moveby {125.0 0.0 0.0}
  }
  $p1 moveby {-330.0 -50.0 0.0}
  $p2 moveby {-330.0 -50.0 0.0}
}
set N_k $k

# Delete the selections and molecules loaded in VMD. 
$p1 delete
$p2 delete
mol delete top
mol delete top


# Use psfgen to generate pdb and psf of the combined system. 
resetpsf
package require psfgen
#Load SBCG topology files for BAR domains, lipids and ions. 
topology lipid-ion.top
topology cg_monomer.top

#Input the coordinates of the six BAR domains for psfgen.
for {set k 1} {$k < $N_k} {incr k} {
  segment P1$k { pdb p1-$k.pdb
  }
  coordpdb p1-$k.pdb P1$k

  segment P2$k { pdb p2-$k.pdb
  }
  coordpdb p2-$k.pdb P2$k
}

#Input the coordinates of membrane for psfgen.
segment L1 { pdb mixture.pdb
}
coordpdb mixture.pdb L1

#Input the coordinates of ions for psfgen.
segment I { pdb ions.pdb
}
coordpdb ions.pdb I

# Guess coordinates (if any are missing) according to the topology files.
guesscoord

# Write pdb and psf of the combined system.
writepdb 6bar.pdb
writepsf 6bar.psf


for {set k 1} {$k < $N_k} {incr k} {
  file delete p1-$k.pdb
}
