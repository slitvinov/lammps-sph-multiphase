# read lammps files with "fix deposit"
log vmd.tcl

# read variable number of atoms
topo readvarxyz [lindex $argv 0]/data.xyz
mol modstyle 0 0 Points 16
pbc set {1.0 1.0 1.0} -all

set sel [atomselect top all]
$sel set radius 0.018

mol modstyle 0 0 VDW 0.600000 15.000000

color Display Background white
color Display FPS black
color Axes Labels black
axes location off

# show only gas
#mol modselect 0 0 (all) and user > 0 and name B
