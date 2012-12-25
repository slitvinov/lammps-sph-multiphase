# For emacs set
# (setenv "PYTHONPATH" "/scratch/work/Pizza.py/src/")
# (setenv "PYTHONPATH" "/home/vital303/work/Pizza.py/src/")

import os.path
import glob
import dump
import numpy as np
import math

os.system("./aux.sh")
dname = "last1"
d = dump.dump(os.path.join(dname, "dump*"))

d.tselect.all()
# list of times
t = d.time()

def interp(xi, x, y, h):
    def kernel(r):
        return math.exp(-r*r/ (h*h) )
        #return (math.fabs(r)<h)+0.0
    yc = np.zeros(np.size(xi))
    wc = np.zeros(np.size(xi))
    for i, xc in enumerate(xi):
        for xp, yp in zip(x, y):
            w = kernel(xp-xc)
            yc[i] += w*yp
            wc[i] += w
    yc[wc>0] = yc[wc>0]/wc[wc>0]
    return yc

d.aselect.test("$type == 2")
for tc in t:
    d.aselect.test("$type == 2")
    xd = d.vecs(tc, "x")
    xcm = np.mean(xd)
    yd = d.vecs(tc, "y")
    ycm = np.mean(yd)
    zd = d.vecs(tc, "z")
    zcm = np.mean(zd)
    d.aselect.all()
    x = d.vecs(tc, "x")
    y = d.vecs(tc, "y")
    z = d.vecs(tc, "z")
    r = np.sqrt(np.square(x-xcm) + np.square(y-ycm) + np.square(z-zcm))
    temperature = d.vecs(tc, "c_it_atom")
    ri = np.linspace(min(r), max(r), 100)
    #data = [r, temperature]
    #Z2 = ndimage.gaussian_filter(data, 0.2, mode='constant');
    temperature_smoothed = interp(ri, r, temperature, 0.01)
    np.savetxt(os.path.join(dname, "com.int.%i" % tc),
               np.column_stack( (ri, temperature_smoothed) ))
    np.savetxt(os.path.join(dname, "com.dat.%i" % tc),
               np.column_stack( (r, temperature) ))
    




