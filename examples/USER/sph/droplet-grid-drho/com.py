# For emacs set
# (setenv "PYTHONPATH" "/scratch/work/Pizza.py/src/")
# (setenv "PYTHONPATH" "/home/vital303/work/Pizza.py/src/")

import dump
import numpy as np
d = dump.dump("data-wall/dump*.dat")

d.tselect.all()
t = d.time()
#d.scale()
# upper right corner of the babel
#d.aselect.test("$x >= 0.5417 and $y >= 0.5417 and $type == 2")

com_droplet = np.zeros([np.size(t), 2])
d.aselect.test("$type == 2")
i = 0
for tc in t:
    x = d.vecs(tc, "x")
    com_droplet[i, 0] = np.mean(x)
    y = d.vecs(tc, "y")
    com_droplet[i, 1] = np.mean(y)
    i = i + 1

com = np.zeros([np.size(t), 2])
i = 0
for tc in t:
    idx = (x>com_droplet[i, 0]) & (y>com_droplet[i, 1])
    x = np.array(d.vecs(tc, "x"))
    y = np.array(d.vecs(tc, "y"))
    com[i, 0] = np.mean(x[idx])
    com[i, 1] = np.mean(y[idx])
    i = i + 1

# com = np.zeros([np.size(t), 2])
# i = 0
# for tc in t:
#     x = d.vecs(tc, "x")
#     com[i, 0] = np.mean(x) - 0.5
#     y = d.vecs(tc, "y")
#     com[i, 1] = np.mean(y) - 0.5
#     i = i + 1

np.savetxt("com.dat", np.column_stack( (t, com) ))
