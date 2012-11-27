log vmd.tcl

mol modstyle 0 0 Points 21.0
set sel_solvent [atomselect top "resname 1"]
$sel_solvent set name water

set sel_solvent [atomselect top "resname 2"]
$sel_solvent set name droplet
