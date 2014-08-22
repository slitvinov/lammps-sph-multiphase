# read lammps files with "fix deposit"
log vmd.tcl
user add key q exit

# read variable number of atoms
set filename [lindex $argv 0]
topo readvarxyz $filename
mol modstyle 0 0 Points 16

set xmax 2.0
set ymax 2.0
set zmax 2.0
pbc set [list ${xmax} ${ymax} ${zmax}] -all

set sel [atomselect top all]
set nx 30
set dx [expr {1.0/${nx}}]
set V  [expr {pow($dx,3)}]
set pi 3.1415926
set r  [expr pow($V/(4.0/3.0*$pi), 1.0/3.0)]
$sel set radius $r


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

proc getn {} {
    global sel
    global ypos
    set n 0
    set vals [$sel get user]
    set ycoords [$sel get y]
    foreach u $vals y $ycoords {
	if {$u>0 && $y>${ypos}} {
	    incr n
	}
    }
    puts $n
}

set linelist {}
proc putline {} {
    global xmax ymax zmax
    global linelist
    global ypos
    foreach id $linelist {
	graphics top delete $id
    }
    set linelist {}
    lappend linelist [graphics top line [list 0 ${ypos} ${zmax}] [list ${xmax} ${ypos} ${zmax}]]
    lappend linelist [graphics top line [list 0 ${ypos} 0] [list ${xmax} ${ypos} 0]]
    lappend linelist [graphics top line [list ${xmax} ${ypos} 0] [list ${xmax} ${ypos} ${zmax}]]
    lappend linelist [graphics top line [list 0 ${ypos} 0] [list 0 ${ypos} ${zmax}]]
}

proc incrypos {} {
    global ypos
    set ypos [expr {$ypos + 0.1}]
    putline
}

proc decrypos {} {
    global ypos
    set ypos [expr {$ypos - 0.1}]
    putline
}

set ypos 0.0
user add key s getn
user add key n decrypos
user add key p incrypos

#mol modstyle 0 0 VDW 0.600000 12.000000
color Name B purple
mol modstyle 0 0 VDW 1.000000 12.000000

mol modselect 0 0 (all) and user > 0 and y<1.5
proc makemov {}  {
set frame 0
while {${frame} < 69} {
    set filename snap.[format "%04d" $frame].tga
    animate goto ${frame}
    incr frame
    render snapshot ${filename}
    if {${frame} == 69} {
	exit
    }
}
}


