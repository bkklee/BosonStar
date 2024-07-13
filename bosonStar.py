import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def diffEq(y, x, mu):
    f, f1, f2, f3 = y
    dydx = [f1, f2, f3, 2.0*f*f - 4.0*f3/x + 10.0*f1*f2/f/x - 6.0*f1*f1*f1/f/f/x + 3.0*f3*f1/f + 2.0*f2*f2/f - 7.0*f1*f1*f2/f/f + 3.0*f1*f1*f1*f1/f/f/f + 2.0*f*mu*(f2+2.0*f1/x)]
    return dydx

mu = 3.0
left = -10.0
right = -1e-9
ans = 0.0

def solver(mu, A2):
    y0 = [1.0, 0.0, A2, 0.0]
    x = np.linspace(1e-20, 30.0, 200000)
    sol = odeint(diffEq, y0, x, args=(mu,), atol=1e-12, rtol=1e-12)
    return x, sol[:, 0]

for i in range(100):
    x, ansLeft = solver(mu, left)
    x, ansMid = solver(mu, (left+right)/2)
    x, ansRight = solver(mu, right)

    print((left+right)/2, x[np.argmin(ansLeft)], x[np.argmin(ansMid)], x[np.argmin(ansRight)])

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
        break
    ans = (left+right)/2

x, ansMid = solver(mu, ans)
#for i in range(len(ansMid)-1):
#    if(ansMid[i+1] <= 0 or abs(ansMid[i+1]/ansMid[i])>2.0):
#        #Too negative
#        print(i, "Too negative")
#        break
#    elif(ansMid[i+1] > ansMid[i]):
#        #Too positive
#        print(i, "Too positive")
#        break
np.savetxt("a.txt", ansMid)

plt.plot(x, ansMid, 'b', label='f(x)')
plt.legend(loc='best')
plt.xlabel('x')
plt.ylim(-0.05,1.05)
plt.grid()
plt.show()