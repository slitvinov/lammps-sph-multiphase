import numpy as np
g = 1.0/np.linspace(1.0/3.0, 1.0/200.0, 20)

np.random.shuffle(g)
np.savetxt("g.dat", g, "%.2f")

#parallel -a g.dat echo
