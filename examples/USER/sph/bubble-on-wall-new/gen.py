import numpy as np
mi=2.0
ma=10.0
g = 1.0/np.linspace(1/mi, 1/ma, 8)

np.random.shuffle(g)
np.savetxt("g.dat", g, "%.3f")

#parallel -a g.dat echo
