import numpy as np
from scipy.integrate import simpson, odeint

hBarCU = 1.1977151493389159e-76

def cleanProfile(profile):
    for i in range(len(profile)-1):
        if(profile[i] <= 0 or profile[i+1] > profile[i]):
            profile[i] = 0.0
            profile[i+1] = 0.0
    profile[len(profile)-1] = 0.0

    return profile

def diffEq(y, x, mu):
    f, f1, f2, f3 = y
    if(f<=0):
        return [0.0, 0.0, 0.0, 0.0]
    else:
        dydx = [f1, f2, f3, 2.0*f*f - 4.0*f3/x + 10.0*f1*f2/f/x - 6.0*f1*f1*f1/f/f/x + 3.0*f3*f1/f + 2.0*f2*f2/f - 7.0*f1*f1*f2/f/f + 3.0*f1*f1*f1*f1/f/f/f + 2.0*f*mu*(f2+2.0*f1/x)]
        return dydx
    
def solver(mu, A2, x=np.linspace(1e-20, 30.0, 100000)):
    y0 = [1.0, 0.0, A2, 0.0]
    sol = odeint(diffEq, y0, x, args=(mu,))
    return x, sol[:, 0]

def findA2(mu):
    left = -10.0
    right = -1e-9
    A2 = 0.0

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
        
        if((left+right)/2 - A2 == 0.0):
            print("A2= ", (left+right)/2, " X=", x[np.argmin(ansMid)])
            break
        else:
            A2 = (left+right)/2
    
    return A2

def getProfile(mu, mBosonCU, aCU):
    A2 = findA2(mu)
    x, profile = solver(mu, A2)
    profile = cleanProfile(profile)

    n0 = pow(1.0/simpson(profile*x*x, x=x), 4.0)
    chi = mu/np.sqrt(n0)
    MCU = np.sqrt(chi*hBarCU*hBarCU/mBosonCU/aCU/4.0)
    b = hBarCU*hBarCU/2.0/MCU/mBosonCU/mBosonCU

    return {"x": x*b/pow(n0,1.0/4.0), "profile": profile*n0*MCU/4.0/np.pi/b/b/b, "MCU": MCU}

def findMu(mass, mBosonCU, aCU):
    x0 = 0.0
    y0 = 0.0 - mass
    x1 = 1.0
    y1 = getProfile(x1, mBosonCU, aCU)["MCU"] - mass

    for _ in range(100):
        newX = x1 - y1*(x1-x0)/(y1-y0)
        newY = getProfile(newX, mBosonCU, aCU)["MCU"] - mass
        x0 = x1
        y0 = y1
        x1 = newX
        y1 = newY
        if(abs(newY) < 1e-5):
            return newX
    
    return x1