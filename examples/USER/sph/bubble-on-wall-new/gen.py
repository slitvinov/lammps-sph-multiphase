import numpy as np
g = 1.0/np.linspace(1.0/16.0, 1.0/50.0, 30)

np.random.shuffle(g)
np.savetxt("g.dat", g, "%.2f")

#parallel -a g.dat echo
