import numpy as np
from scipy.integrate import simpson, odeint
import matplotlib.pyplot as plt

mBosonCU = 2e-77
aCU = 1e-73
hBarCU = 1.1977151493389159e-76

def diffEq(y, x, mu):
    f, f1, f2, f3 = y
    if(f<=0):
        return [0.0, 0.0, 0.0, 0.0]
    else:
        dydx = [f1, f2, f3, 2.0*f*f - 4.0*f3/x + 10.0*f1*f2/f/x - 6.0*f1*f1*f1/f/f/x + 3.0*f3*f1/f + 2.0*f2*f2/f - 7.0*f1*f1*f2/f/f + 3.0*f1*f1*f1*f1/f/f/f + 2.0*f*mu*(f2+2.0*f1/x)]
        return dydx

def solver(mu, A2):
    y0 = [1.0, 0.0, A2, 0.0]
    x = np.linspace(1e-20, 30.0, 100000)
    sol = odeint(diffEq, y0, x, args=(mu,))
    return x, sol[:, 0]

def findA2(mu):
    left = -10.0
    right = -1e-9
    ans = 0.0

    for i in range(100):
        x, ansMid = solver(mu, (left+right)/2)

        for i in range(len(ansMid)-1):
            if(ansMid[i+1] <= 0 or abs(ansMid[i+1]/ansMid[i])>2.0):
                #Too negative
                left = (left+right)/2
                break
            elif(ansMid[i+1] > ansMid[i]):
                #Too positive
                right = (left+right)/2
                break
        
        if((left+right)/2 - ans == 0.0):
            print("A2= ", (left+right)/2, " X=", x[np.argmin(ansMid)])
            break
        else:
            ans = (left+right)/2
    
    return ans

def findMass(mu):
    A2 = findA2(mu)
    x, ans = solver(mu, A2)
    for i in range(len(ans)-1):
        if(ans[i] <= 0 or ans[i+1] > ans[i]):
            ans[i] = 0.0
            ans[i+1] = 0.0
    ans[len(ans)-1] = 0.0

    n0 = pow(1.0/simpson(ans*x*x, x=x), 4.0)
    chi = mu/np.sqrt(n0)
    MCU = np.sqrt(chi*hBarCU*hBarCU/mBosonCU/aCU/4.0)

    return MCU

def findMu(mass):
    x0 = 0.0
    y0 = 0.0 - mass
    x1 = 1.0
    y1 = findMass(x1) - mass

    for _ in range(100):
        newX = x1 - y1*(x1-x0)/(y1-y0)
        newY = findMass(newX) - mass
        x0 = x1
        y0 = y1
        x1 = newX
        y1 = newY
        if(newY < 1e-6):
            return newX
    
    return x1

mu = findMu(0.1)
print("MU=", mu)
mass = findMass(mu)
print("MASS=", mass)
A2 = findA2(mu)
print("A2=", A2)
x, ans = solver(mu, A2)
np.savetxt("a.txt", ans)

plt.plot(x, ans, 'b', label='f(x)')
plt.legend(loc='best')
plt.xlabel('x')
plt.ylim(-0.05,1.05)
plt.grid()
plt.show()