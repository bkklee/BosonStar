import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def diffEq(y, x, mu):
    f, f1, f2, f3 = y
    dydx = [f1, f2, f3, 2.0*f*f - 4.0*f3/x + 10.0*f1*f2/f/x - 6.0*f1*f1*f1/f/f/x + 3.0*f3*f1/f + 2.0*f2*f2/f - 7.0*f1*f1*f2/f/f + 3.0*f1*f1*f1*f1/f/f/f + 2.0*f*mu*(f2+2.0*f1/x)]
    return dydx

mu = 3.0
y0 = [1.0, 0.0, -0.216, 0.0]

x = np.linspace(0.0001, 25.0, 100000)
sol = odeint(diffEq, y0, x, args=(mu,))

plt.plot(x, sol[:, 0], 'b', label='f(x)')
plt.legend(loc='best')
plt.xlabel('x')
plt.grid()
plt.show()