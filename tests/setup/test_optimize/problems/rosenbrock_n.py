
import numpy as np

n = 1e8

def model_init():
    return 0.1*np.ones(n)

def model_true():
    return np.ones(n)

def func(x):
    return sum(100*(x[:-1]**2. - x[1:])**2. + 
               (x[:-1] - 1.)**2.)

def grad(x):
    g = np.zeros(n)

    g[1:-1] = -200*(x[:-2]**2. - x[1:-1]) + \
               400.*x[1:-1]*(x[1:-1]**2. - x[2:]) + \
               2.*(x[1:-1]-1.)

    g[0] = 400.*x[0]*(x[0]**2. - x[1]) + \
           2.*(x[0] - 1)

    g[-1] = -200.*(x[-2]**2. - x[-1])

    return g

