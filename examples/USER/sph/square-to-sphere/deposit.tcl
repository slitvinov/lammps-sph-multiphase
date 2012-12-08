# read lammps files with "fix deposit"
log vmd.tcl

# read variable number of atoms
topo readvarxyz [lindex $argv 0]/data.xyz
pbc set {1.0 1.0 1.0} -all
pbc box

set sel [atomselect top all]
$sel set radius 0.05

mol modselect 0 0 (all) and user > 0 and name B
mol modstyle 0 0 VDW 0.600000 12.000000

# show only gas
