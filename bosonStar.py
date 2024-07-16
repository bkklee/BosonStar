import numpy as np
from scipy.integrate import simpson, odeint
import matplotlib.pyplot as plt

from diffEqSolver import *

mBosonCU = 2e-77
aCU = 1e-73
hBarCU = 1.1977151493389159e-76
targetMass = 0.26843

mu = findMu(targetMass, mBosonCU, aCU)
print("MU=", mu)
allThings = getProfile(mu, mBosonCU, aCU)
print("MASS=", allThings["MCU"])
x = allThings["x"]
ans = allThings["profile"]

np.savetxt("a.txt", ans)
plt.plot(x, ans, 'b', label='f(x)')
plt.legend(loc='best')
plt.xlabel('x')
plt.grid()
plt.show()