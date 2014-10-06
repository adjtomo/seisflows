
import numpy as np


def backtrack2(f0,g0,x1,f1,b1=0.1,b2=0.5):

    # parabolic backtrack
    x2 = -g0*x1**2/(2*(f1-f0-g0*x1))

    # apply constraints
    if x2 > b2*x1:
        x2 = b2*x1
    elif x2 < b1*x1:
        x2 = b1*x1
    return x2

def backtrack3(f0,g0,x1,f1,x2,f2):
    # cubic backtrack
    raise NotImplementedError

def polyfit2(x,f):
    # parabolic fit
    i = np.argmin(f)
    p = np.polyfit(x[i-1:i+2],f[i-1:i+2],2)

    if p[0] > 0:
        return -p[1]/(2*p[0])
    else:
        print -1
        raise Exception()

def lsq2(x,f):
    # parabolic least squares fit
    p = np.polyfit(x,f,2)
    if p[0] > 0:
        return -p[1]/(2*p[0])
    else:
        print -1
        raise Exception()
