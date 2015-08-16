resetpsf
package require psfgen

topology lipid-ion.top

segment L1 { pdb mixture.pdb
}
coordpdb mixture.pdb L1

guesscoord

writepdb mixture.pdb
writepsf mixture.psf
