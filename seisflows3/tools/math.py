#!/usr/bin/env python
"""
Mathematical tools for Seisflows
"""
import numpy as np

from scipy.signal import hilbert as analytic

infinity = np.inf


def gauss2(X, Y, mu, sigma, normalize=True):
    """
    Evaluates Gaussian over points of X, Y

    :type X: np.ndarray
    :param X: x-axis to evaluate gaussian over
    :type Y: np.ndarray
    :param Y: y-axis to evaluate gaussian over
    :type mu: ???
    :param mu: expected value
    :type sigma: ???
    :param sigma: standard deviation
    :type normalize: bool
    :param normalize: normalize the results
    """
    # evaluates Gaussian over X,Y
    D = (sigma[0, 0] * sigma[1, 1]) - (sigma[0, 1] * sigma[1, 0])
    B = np.linalg.inv(sigma)

    X = X - mu[0]
    Y = Y - mu[1]
    Z = (B[0, 0] * X ** 2. +
         B[0, 1] * X * Y +
         B[1, 0] * X * Y +
         B[1, 1] * Y ** 2.)
    Z = np.exp(-0.5 * Z)

    if normalize:
        Z *= (2. * np.pi * np.sqrt(D)) ** -1.

    return Z


def backtrack2(f0, g0, x1, f1, b1=0.1, b2=0.5):
    """
    Safeguarded parabolic backtrack

    Note for equation look to
        Nocedal & Wright, 2006 ??

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
    """
    # Parabolic backtrack
    x2 = -g0 * x1 ** 2 / (2 * (f1 - f0 - g0 * x1))

    # Apply safeguards
    if x2 > b2 * x1:
        x2 = b2 * x1
    elif x2 < b1 * x1:
        x2 = b1 * x1

    return x2


def backtrack3(f0, g0, x1, f1, x2, f2):
    """
    Safeguarded cubic backtrack
    """
    raise NotImplementedError


def polyfit2(x, f):
    """
    Parabolic line fitting

    :type x: np.array
    :param x: x coordinates
    :type f: np.array
    :param f: y coordinates
    """
    # parabolic fit
    i = np.argmin(f)
    p = np.polyfit(x[i-1:i+2], f[i-1:i+2], 2)

    if p[0] > 0:
        return -p[1] / (2 * p[0])
    else:
        print(-1)
        raise Exception()


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


def angle(x, y):
    """
    Determine the angle between two vectors using dot products

    :type x: np.array
    :param x: vector 1
    :type y: np.array
    :param y: vector 2
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
    """
    return np.dot(np.squeeze(x), np.squeeze(y))


def hilbert(w):
    """
    Take the Hilbert transform of some function to get the analytic signal

    :type w: np.array
    :param w: signal data, must be real
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


