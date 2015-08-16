#!/usr/local/bin/wish
source /Users/Madsci/Desktop/Lab_Primary/Dynamics/NAMD/My_Scripts/calcBC_Proc.tcl
#MinMax = {-11.192999839782715 67.40699768066406 76.11000061035156} {81.64700317382813 129.11599731445313 138.61000061035156}
#Center = 35.2236328125 98.40187072753906 107.4239501953125


entry .min -width 25 -relief sunken -textvariable min
entry .max -width 25 -relief sunken -textvariable max
entry .cen -width 50 -relief sunken -textvariable cen
text .info -width 60 -height 10
text .result -width 70 -height 20
text .setup -width 100 -height 60
#-----------Useful Info Print -----------
.info insert 1.0 "Useful Info:\n\n"
.info insert end "NAMD Start Mac:  ../../../Program/NAMD_2.7_MacOSX-x86/namd2 +p 4\n"
.info insert end "NAMD Out   Mac:  ../../../Data/Trajectories\n\n"
.info insert end "Cluster Run:     /usr/local/bin/qsub -d /common/madsci/Dynamics/NAMD/Run/temp -V /common/madsci/Dynamics/NAMD/Run/Run_Namd\n\n"
.result insert end "Boundary Conditions:\n\n"
#----------------------------------------
button .calculate -text "Calculate BC" \
	-command {
		set result [calcBC $min $max]
		set center [split $cen]
		#---------------------CellBasisVector------------------
		.result insert end [format "cellBasisVector1\t%.2f 0 0\n" [lindex $result 0]]
		.result insert end [format "cellBasisVector2\t0 %.2f 0\n" [lindex $result 1]]
		.result insert end [format "cellBasisVector3\t0 0 %.2f\n" [lindex $result 2]]
		#-------------------------CENTER-----------------------
		.result insert end [format "cellOrigin\t%.2f %.2f %.2f\n" [lindex $cen 0] [lindex $cen 1] [lindex $cen 2]] 
		}
	

#---------Config Generator -----------
label .pIn -text "Enter PSF and PDB filename:"
entry .protein -width 25 -relief sunken -textvariable protein

label .pathInText -text "Enter the PDB/PSF Folder (Within /PDBs):"
entry .pathIn -width 25 -relief sunken -textvariable pathIn

set input "../../../PDBs"
set output "../../../Data/Trajectories"

button .rigid -text "NPT Rigid" \
	-command {
		.setup insert end "NPT Rigid:\n"
		.setup insert end "structure          ${input}${pathIn}/${protein}.pdb\n"
		.setup insert end "coordinates        ${input}${pathIn}/${protein}.psf\n"
		.setup insert end "set output         ${output}${pathIn}/${protein}_Rigid\n"
}
button .flex -text "NPT Flex" \
	-command {
		.setup insert end "\nNPT Flex:\n"
		.setup insert end "structure          ${input}${pathIn}/${protein}.pdb\n"
		.setup insert end "coordinates        ${input}${pathIn}/${protein}.psf\n"
		.setup insert end "set output         ${output}${pathIn}/${protein}_Flex\n"
		.setup insert end "set inputname      ${output}${pathIn}/${protein}_Rigid\n"
}
button .nve -text "NVE" \
	-command {
		.setup insert end "\nNVE Full:\n"
		.setup insert end "structure          ${input}${pathIn}/${protein}.pdb\n"
		.setup insert end "coordinates        ${input}${pathIn}/${protein}.psf\n"
		.setup insert end "set output         ${output}${pathIn}/${protein}_NVE\n"
		.setup insert end "set inputname      ${output}${pathIn}/${protein}_Flex\n"
}

.setup insert 1.0 "Config File Builder:\n\n"


grid .min .calculate .max -padx 1m -pady 1m
grid .cen - - -padx 1m -pady 1m
grid .pIn .protein -padx 1m -pady 1m
grid .pathInText .pathIn -padx 1m -pady 1m
grid .rigid .flex .nve -padx 1m -pady 1m

grid .info - - -padx 1m -pady 1m
grid .result - - -padx 1m -pady 1m
grid .setup - - -padx 1m -pady 1m


