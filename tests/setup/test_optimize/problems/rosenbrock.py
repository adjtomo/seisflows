
import numpy as np

def model_init():
    return np.array([-1.2, 1])

def model_true():
    return np.array([1, 1])

def func(x):
    return np.array([((1-x[0])**2 + 100*(-x[0]**2+x[1])**2)])

def grad(x):
    return np.array([-2*(1-x[0]) - 400*x[0]*(-x[0]**2+x[1]),
                     200*(- x[0]**2+x[1])])

