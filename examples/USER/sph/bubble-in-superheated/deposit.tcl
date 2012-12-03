# read lammps files with "fix deposit"
log vmd.tcl

# read variable number of atoms
topo readvarxyz data-wall/data.xyz
mol modstyle 0 0 Points 16

#mol modselect 0 0 (all) and user > 0 and name B