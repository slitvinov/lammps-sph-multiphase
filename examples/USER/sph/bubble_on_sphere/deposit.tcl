# read lammps files with "fix deposit"
log vmd.tcl
user add key q exit

# read variable number of atoms
set filename [lindex $argv 0]
topo readvarxyz $filename
mol modstyle 0 0 Points 16

pbc set [list ${xmax} ${ymax} ${zmax}] -all

set sel [atomselect top all]
set nx  30
set dx [expr {1.0/${nx}}]
set V  [expr {pow($dx,3)}]
set pi 3.1415926
set r  [expr pow($V/(4.0/3.0*$pi), 1.0/3.0)]
$sel set radius $r

#mol modstyle 0 0 VDW 0.600000 12.000000
color Name A blue
color Name B purple

mol modstyle 0 0 VDW 1.000000 12.000000

mol modselect 0 0 (all) and user > 0

