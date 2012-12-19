import numpy as np
mi=0.1
ma=5.0
g = 1.0/np.linspace(1/ma, 1/mi, 100)

np.random.shuffle(g)
np.savetxt("g.dat", g, "%.2f")

#parallel -a g.dat echo
