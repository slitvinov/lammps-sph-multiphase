log vmd.tcl

mol modstyle 0 0 Points 21.0
set sel_solvent [atomselect top "resname 1"]
$sel_solvent set name water

user add key q exit

set sel_solvent [atomselect top "resname 2"]
$sel_solvent set name droplet


set sel_solvent [atomselect top "resname 3"]
$sel_solvent set name boundary
