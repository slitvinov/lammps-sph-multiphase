import numpy as np
mi=2e5
ma=9e5
g = 1.0/np.linspace(1/mi, 1/ma, 8)

#np.random.shuffle(g)
np.savetxt("g.dat", g, "%6.2e")

#parallel -a g.dat echo
