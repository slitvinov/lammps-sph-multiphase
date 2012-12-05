# read lammps files with "fix deposit"
log vmd.tcl

# read variable number of atoms
topo readvarxyz data-wall/data.xyz
mol modstyle 0 0 Points 16
pbc set {1.0 1.0 1.0} -all

set sel [atomselect top all]
$sel set radius 0.03

mol modselect 0 0 (all) and user > 0 and name B
mol modstyle 0 0 VDW 0.600000 12.000000

pbc box

# show only gas
#mol modselect 0 0 (all) and user > 0 and name B
