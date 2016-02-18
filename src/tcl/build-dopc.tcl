resetpsf
package require psfgen

topology lipid-ion.top

segment L1 { pdb dopc.pdb
}
coordpdb dopc.pdb L1

guesscoord

writepdb dopc.pdb
writepsf dopc.psf
