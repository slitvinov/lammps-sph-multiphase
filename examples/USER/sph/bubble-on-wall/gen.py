import numpy as np
g = 1.0/np.linspace(10.0, 20.0, 20)

np.random.shuffle(g)
np.savetxt("g.dat", g, "%.3f")

#parallel -a g.dat echo
