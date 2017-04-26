### Open the log file for reading and the output .dat file for writing
set file [open ubq_ww_pcv.log r]
set output [open analysis/ft.dat w]

### Gather input from user.
puts "Enter a value for n_x:"
set nx [gets stdin]
puts "Enter a value for n_y:"
set ny [gets stdin]
puts "Enter a value for n_z:"
set nz [gets stdin]

### Loop over all lines of the log file
set file [open ubq_ww_pcv.log r]
while { [gets $file line] != -1 } {

### Determine if a line contains SMD output. If so, write the 
### timestep followed by f(dot)n to the output file
  if {[lindex $line 0] == "SMD"} {
    puts $output "[lindex $line 1] [expr $nx*[lindex $line 5] + $ny*[lindex $line 6] + $nz*[lindex $line 7]]"
    }
  }

### Close the log file and the output .dat file
close $file
close $output
