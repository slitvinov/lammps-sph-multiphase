import numpy as np
mi=0.1
ma=5.0
g = 1.0/np.linspace(0.5, 1.4, 100)

np.random.shuffle(g)
np.savetxt("g.dat", g, "%.3f")

#parallel -a g.dat echo
