# Build CG structure.
package require psfgen
resetpsf

topology cg_monomer.top

segment P1 { pdb cg_monomer-2.pdb
}
coordpdb cg_monomer-2.pdb P1

guesscoord

writepdb cg_monomer-2-psfgen.pdb
writepsf cg_monomer-2-psfgen.psf
