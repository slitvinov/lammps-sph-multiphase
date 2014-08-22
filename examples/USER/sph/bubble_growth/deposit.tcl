# read lammps files with "fix deposit"
log vmd.tcl
user add key q exit

# read variable number of atoms
topo readvarxyz [lindex $argv 0]
pbc set {1.0 1.0 1.0} -all

set sel [atomselect top all]
$sel set radius 0.05

#  show region of interest in the center
proc roi {dw} {
    set x1 [expr {0.5-$dw}]
    set x2 [expr {0.5+$dw}]
    set y1 [expr {0.5-$dw}]
    set y2 [expr {0.5+$dw}]
    set z1 [expr {0.5-$dw}]
    set z2 [expr {0.5+$dw}]
    mol modselect 0 0 (all) and user > 0 and x>$x1 and x<$x2 and y>$y1 and y<$y2 and z>$z1 and z<$z2
}

proc slice {dw {dim x}} {
    set z1 [expr {0.5-$dw}]
    set z2 [expr {0.5+$dw}]
    mol modselect 0 0 (all) and user > 0 and ${dim}>$z1 and ${dim}<$z2
}

#mol modstyle 0 0 VDW 0.300000 12.000000
mol modstyle 0 0 Points 16
mol modselect 0 0 (all) and user > 0
#pbc box

# show only bubble
#
