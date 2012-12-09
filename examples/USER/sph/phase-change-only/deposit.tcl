# read lammps files with "fix deposit"
log vmd.tcl

# read variable number of atoms
topo readvarxyz [lindex $argv 0]/data.xyz
mol modstyle 0 0 Points 16

set sel [atomselect top all]
$sel set radius 0.03

mol modstyle 0 0 VDW 0.600000 12.000000
