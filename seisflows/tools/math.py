#!/usr/bin/env python3
"""
Mathematical tools for Seisflows
"""
import sys
import numpy as np
from scipy.signal import hilbert as analytic

from seisflows import logger
from seisflows.tools import msg


def angle(x, y):
    """
    Determine the angle between two vectors using dot products

    :type x: np.array
    :param x: vector 1
    :type y: np.array
    :param y: vector 2
    :rtype: float
    :return: the angle in degrees between `x` and `y`
    """
    xy = dot(x, y)
    xx = dot(x, x)
    yy = dot(y, y)
    return np.arccos(xy / (xx * yy) ** 0.5)


def dot(x, y):
    """
    Calculate the dot product between two vectors

    :type x: np.array
    :param x: vector 1
    :type y: np.array
    :param y: vector 2
    :rtype: float
    :return: The dot product between `x` and `y`
    """
    return np.dot(np.squeeze(x), np.squeeze(y))


def hilbert(w):
    """
    Take the Hilbert transform of some function to get the analytic signal

    TODO Change the naming here, it seems confusing to rename a scipy function
    TODO and then overwrite its name with this function.

    :type w: np.array
    :param w: signal data, must be real
    :rtype: float
    :return: imaginary part of the analytic signal
    """
    return np.imag(analytic(w))


def poissons_ratio(vp, vs):
    """
    Calculate Poisson's Ratio based on the definition given in the Specfem3D
    source code

    :type vp: float or np.array
    :param vp: P-wave velocity
    :type vs: float or np.array
    :param vs: S-wave velocity
    :rtype: float or np.array
    :return: Poissons ratio
    """
    return 0.5 * (vp * vp - 2 * vs * vs) / (vp * vp - vs * vs)


def parabolic_backtrack(f0, g0, x1, f1, b1=0.1, b2=0.5):
    """
    Safeguarded parabolic backtracking function
    Equation provided in Nocedal & Wright, 2006 ??

    :type f0: float
    :param f0: initial misfit function value
    :type g0: float
    :param g0: slope
    :type x1: float
    :param x1: step length value
    :type f1: float
    :param f1: current misfit function value (?)
    :type b1: float
    :param b1: constant for safeguard
    :type b2: float
    :param b2: constant for safeguard
    :rtype: float
    :return: trial step length (alpha)
    """
    # Parabolic backtrack
    x2 = -g0 * x1 ** 2 / (2 * (f1 - f0 - g0 * x1))

    # Apply safeguards
    if x2 > b2 * x1:
        x2 = b2 * x1
    elif x2 < b1 * x1:
        x2 = b1 * x1

    return x2


def gaussian(x, y, mu, sigma, normalize=True):
    """
    Evaluates Gaussian over points of X, Y

    :type x: np.ndarray
    :param x: x-axis to evaluate gaussian over
    :type y: np.ndarray
    :param y: y-axis to evaluate gaussian over
    :type mu: ???
    :param mu: expected value
    :type sigma: ???
    :param sigma: standard deviation
    :type normalize: bool
    :param normalize: normalize the results
    """
    # evaluates Gaussian over X,Y
    d = (sigma[0, 0] * sigma[1, 1]) - (sigma[0, 1] * sigma[1, 0])
    b = np.linalg.inv(sigma)

    x = x - mu[0]
    y = y - mu[1]
    z = (b[0, 0] * x ** 2. +
         b[0, 1] * x * y +
         b[1, 0] * x * y +
         b[1, 1] * y ** 2.)
    z = np.exp(-0.5 * z)

    if normalize:
        z *= (2. * np.pi * np.sqrt(d)) ** -1.

    return z


def polynomial_fit(x, f):
    """
    Least squares (polynomial) line fitting used to fit a line to the
    objective function.

    :type x: np.array
    :param x: trial step lengths
    :type f: np.array
    :param f: misfit values
    :rtype: float
    :return: trial step length (alpha)
    """
    i = np.argmin(f)
    p = np.polyfit(x[i-1:i+2], f[i-1:i+2], 2)

    if p[0] <= 0:
        logger.critical(msg.cli("Polynomial line fitting returned a negative "
                                "p[0] value which signifies a negative misfit "
                                "and is not allowed."))
        sys.exit(-1)

    return -p[1] / (2 * p[0])


# The below functions were included in the original SeisFlows package but
# currently serve no purpose in SeisFlows. They are retained for legacy.
def lsq2(x, f):
    """
    Parabolic least squares fit

    :type x: np.array
    :param x: x coordinates
    :type f: np.array
    :param f: y coordinates
    """
    p = np.polyfit(x, f, 2)
    if p[0] > 0:
        return -p[1]/(2*p[0])
    else:
        print -1
        raise Exception()


def nabla(V, h=[]):
    """
    Finite Differences

    Returns sum of first-order spatial derivatives of a function defined on
    a 2D rectangular grid; generalizes Laplacian
    """
    W = np.zeros(V.shape)

    if h==[]:
       h = np.ones((V.ndim, 1))

    # Interior
    W[1:-1, 1:-1] += (V[1:-1, 2:] - V[1:-1, :-2]) / (2.*h[0])
    W[1:-1, 1:-1] += (V[2:, 1:-1] - V[:-2, 1:-1]) / (2.*h[1])

    # Top/bottom edges
    W[0, 1:-1] = ((V[1, 1:-1] - V[0, 1:-1]) / h[1] +
                  (V[0, 2:] - V[0, :-2]) / (2. * h[0]))
    W[-1, 1:-1] = ((V[-1, 1:-1] - V[-2, 1:-1]) / h[1] +
                   (V[-1, 2:] - V[-1, :-2]) / (2. * h[0]))

    # Left/right edges
    W[1:-1, 0] = (V[2:, 0] - V[:-2, 0]) / (2. * h[1]) + (V[1:-1, 1] - V[1:-1, 0]) / h[0]
    W[1:-1, -1] = (V[2:, -1] - V[:-2, -1]) / (2. * h[1]) + (V[1:-1, -1] - V[1:-1, -2]) / h[0]

    # Corners
    W[0, 0] = (V[1, 0] - V[0, 0]) / h[1] + (V[0, 1] - V[0, 0]) / h[0]
    W[0, -1] = (V[1, -1] - V[0, -1]) / h[1] + (V[0, -2] - V[0, -1]) / h[0]
    W[-1, 0] = (V[-2, 0] - V[-1, 0]) / h[1] + (V[-1, 1] - V[-1, 0]) / h[0]
    W[-1, -1] = (V[-1, -1] - V[-2, -1]) / h[1] + (V[-1, -1] - V[-1, -2]) / h[0]

    return W


def nabla2(V, h=[]):
    """
    Finite Differences

    Returns sum of second-order spatial derivatives of a function defined on
    a 2D rectangular grid; generalizes Laplacian
    """
    W = np.zeros(V.shape)

    if h==[]:
       h = np.ones((V.ndim, 1))

    # interior
    W[1:-1,1:-1] += (V[1:-1,2:] -2.*V[1:-1,1:-1] + V[1:-1,:-2])/h[0]**2
    W[1:-1,1:-1] += (V[2:,1:-1] -2.*V[1:-1,1:-1] + V[:-2,1:-1])/h[1]**2

    # left/right edges
    W[0,1:-1] = W[1,1:-1]
    W[-1,1:-1] = W[-2,1:-1]

    # top/bottom edges
    W[0,1:-1] = W[1,1:-1]
    W[-1,1:-1] = W[-2,1:-1]

    # corners
    W[0,0] = (W[0,1] + W[1,0])/2
    W[0,-1] = (W[0,-2] + W[1,-1])/2
    W[-1,0] = (W[-1,1] + W[-2,0])/2
    W[-1,-1] = (W[-1,-2] + W[-2,-1])/2

    return W


def grad(V, h=[]):
    """
    Finite Differences

    Evaluates derivatives on a 2D rectangular grid
    """
    ny, nx = V.shape

    X = np.zeros((ny, nx))
    Y = np.zeros((ny, nx))

    if h==[]:
       h = np.ones((V.ndim, 1))

    # interior
    X[:,1:-1] = (V[:,2:] - V[:,:-2])/(2.*h[0])
    Y[1:-1,:] = (V[2:,:] - V[:-2,:])/(2.*h[1])

    # left/right edges
    X[:,0] = (V[:,1] - V[:,0])/h[1]
    X[:,-1] = (V[:,-1] - V[:,-2])/h[1]

    # top/bottom edges
    Y[0,:] = (V[1,:] - V[0,:])/h[0]
    Y[-1,:] = (V[-1,:] - V[-2,:])/h[0]

    return X,Y


def tv(Z, h=[], epsilon=1.e-6):
    """
    Finite Differences
    """
    nrow = Z.shape[0]
    ncol = Z.shape[1]

    Zx = (Z[:,1:] - Z[:,:-1])/h[0]
    Zy = (Z[1:,:] - Z[:-1,:])/h[1]

    top = np.zeros((nrow, ncol))
    bot = np.zeros((nrow, ncol))

    top[:,1:] += Zx
    top[1:,:] += Zy
    top[ :,-1] += Zx[:,-1]
    top[-1, :] += Zy[-1,:]

    bot[:,1:] += Zx**2
    bot[1:,:] += Zy**2
    bot[ :,-1] += Zx[:,-1]**2
    bot[-1, :] += Zy[-1,:]**2

    return top/(bot + epsilon*bot.max())**0.5


