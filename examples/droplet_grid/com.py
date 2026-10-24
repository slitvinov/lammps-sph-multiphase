# For emacs set
# (setenv "PYTHONPATH" "<Pizza>/src/")

import dump
import numpy as np
d = dump.dump("data-wall-nx40/dump*.dat")

d.tselect.all()
t = d.time()

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
vom = np.zeros([np.size(t), 2])
i = 0
for tc in t:
    idx = (x>com_droplet[i, 0]) & (y>com_droplet[i, 1])
    x = np.array(d.vecs(tc, "x"))
    y = np.array(d.vecs(tc, "y"))
    vx = np.array(d.vecs(tc, "vx"))
    vy = np.array(d.vecs(tc, "vy"))
    com[i, 0] = np.mean(x[idx]) - com_droplet[i, 0]
    com[i, 1] = np.mean(y[idx]) - com_droplet[i, 0]
    vom[i, 0] = np.mean(vx[idx])
    vom[i, 1] = np.mean(vy[idx])

    i = i + 1

np.savetxt("com.dat", np.column_stack( (t, com) ))
np.savetxt("vom.dat", np.column_stack( (t, vom) ))
