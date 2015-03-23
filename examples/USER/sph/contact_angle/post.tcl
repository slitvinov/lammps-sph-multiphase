log vmd.tcl

mol modstyle 0 0 Points 5.0 0.0 0.0
user add key q exit

set filename [lindex $argv 0]
mol new ${filename} filebonds off autobonds off 
