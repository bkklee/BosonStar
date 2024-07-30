import numpy as np
import matplotlib.pyplot as plt
from diffEqSolver import *

hBarCU = 1.1977151493389159e-76

def outputProfile(mBosonCU, aCU, targetMass, boxSize):
    mu = findMu(targetMass, mBosonCU, aCU)
    allThings = getProfile(mu, mBosonCU, aCU)
    print("MU=", mu, "MASS=", allThings["MCU"])
    x = allThings["x"]
    ans = allThings["profile"]

    #Find farpoint
    farPoint = 0
    for i in range(len(ans)):
        if ans[i] == 0.0:
            farPoint = i
            print("FAR POINT:", farPoint, len(ans))
            break

    simulationBox = boxSize*2

    sli = int((farPoint*2.0)//simulationBox)
    x = x[::sli]
    ans = ans[::sli]
    x = x[:simulationBox]
    ans = ans[:simulationBox]

    plt.plot(x, ans, 'b', label='f(x)')
    plt.legend(loc='best')
    plt.xlabel('x')
    plt.grid()
    plt.show()

    return {"x": x[1::2], "profile": ans[1::2], "MCU": allThings["MCU"], "mBoson": mBosonCU, "aCU": aCU}

#outputProfile(2e-77, 1e-73, 0.1)