import numpy as np
g = 1.0/np.linspace(0.15, 0.3, 20)

np.random.shuffle(g)
np.savetxt("g.dat", g, "%.2f")

parallel -a g.dat echo
