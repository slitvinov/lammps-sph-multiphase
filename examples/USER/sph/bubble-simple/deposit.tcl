# read lammps files with "fix deposit"
log vmd.tcl
user add key q exit

# read variable number of atoms
topo readvarxyz [lindex $argv 0]/data.xyz
mol modstyle 0 0 Points 16
pbc set {1.0 1.0 1.0} -all

set sel [atomselect top all]
$sel set radius 0.03

mol modselect 0 0 (all) and user > 0 and z>0.42 and z<0.58
#mol modstyle 0 0 VDW 0.600000 12.000000

pbc box

# show only gas
#mol modselect 0 0 (all) and user > 0 and name B
