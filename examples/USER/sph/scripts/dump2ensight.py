#!/usr/bin/python
"""
  function: translate lammps dump file into ensight format
  output: ensight formated data
  usage: dump2ensight [dump_file]
  NOTE set $PYTHONPATH to .../pizza/src
  Author: Jones, Reese
  http://lammps.sandia.gov/threads/msg20083.html
"""
import sys
import re
from dump import dump
from ensight import ensight
if not globals().has_key("argv"): argv = sys.argv

# main script
if len(argv) < 1:
  raise StandardError, "Syntax: dump2ensight.py dump.1 ..."

# read column:name dictionary and massage names
atoms_pattern = re.compile('ITEM: ATOMS')
file1 = open(argv[1],"r");
line = file1.readline()
while line and not atoms_pattern.match(line) :
  line = file1.readline()
labels = line.split()
labels.remove("ITEM:")
labels.remove("ATOMS")
print "found: ", ( ' '.join(labels) )
i = 1
pairs = []
names = []
for label in labels:
  label  = label.replace('[','_')
  label  = label.replace(']','')
  label  = label.replace('f_','')
  label  = label.replace('c_','')
  label  = label.replace('v_','')
  pairs.append(i)
  pairs.append(label)
  names.append(label)
  names.append(label)
  i=i+1
arg_pairs = tuple(pairs)
arg_names = tuple(names)

# read files
files = ' '.join(argv[1:])
prefix = argv[1]
d = dump(files)
# map columns to names
d.map(*arg_pairs)
e = ensight(d)
e.change = 1 # moving coord
e.one(prefix,*arg_names)

