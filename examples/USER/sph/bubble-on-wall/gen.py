import numpy as np
g = 1.0/np.linspace(1.0/10.0, 1.0/4.0, 7)

np.random.shuffle(g)
np.savetxt("g.dat", g, "%.2f")

#parallel -a g.dat echo
