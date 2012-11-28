# For emacs set
# (setenv "PYTHONPATH" "/scratch/work/Pizza.py/src/")
# (setenv "PYTHONPATH" "/home/vital303/work/Pizza.py/src/")

import dump
import numpy
d = dump.dump("data-wall/dump*.dat")

d.tselect.all()
t = d.time()
d.scale()
d.aselect.test("$x > 0.5 and $y > 0.5 and $type == 2")

com = numpy.zeros(numpy.size(t))
i = 0
for tc in t:
    y = d.vecs(tc, "y")
    com[i] = numpy.mean(y)
    i = i + 1

numpy.savetxt("hm", com)
