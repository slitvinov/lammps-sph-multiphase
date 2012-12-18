import numpy as np
min=5.0
max=20.0
g = 1.0/np.linspace(1.0/max, 1.0/min, 10)

np.random.shuffle(g)
np.savetxt("g.dat", g, "%.2f")

#parallel -a g.dat echo
