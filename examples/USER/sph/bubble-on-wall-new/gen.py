import numpy as np
min=0.1
max=100.0
g = 1.0/np.linspace(1/100.0, 1/1.22, 100)

np.random.shuffle(g)
np.savetxt("g.dat", g, "%.2f")

#parallel -a g.dat echo
