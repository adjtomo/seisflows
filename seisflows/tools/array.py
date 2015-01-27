import os

import numpy as np
import scipy.signal as _signal
import scipy.interpolate as _interp


def gridplot(z):
    """Plots values on 2-D rectangular grid."""

    import pylab

    pylab.pcolor(z)
    pylab.show()


def gridsmooth(Z, span):
    """Smooths values on 2-D rectangular grid."""

    import warnings

    warnings.filterwarnings('ignore')

    x = np.linspace(-2.*span, 2.*span, 2.*span + 1.)
    y = np.linspace(-2.*span, 2.*span, 2.*span + 1.)
    (X, Y) = np.meshgrid(x, y)
    mu = np.array([0., 0.])
    sigma = np.diag([span, span])**2.
    F = gauss2(X, Y, mu, sigma)
    F = F/np.sum(F)
    W = np.ones(Z.shape)
    Z = _signal.convolve2d(Z, F, 'same')
    W = _signal.convolve2d(W, F, 'same')
    Z = Z/W
    return Z


def meshplot(x, y, z):
    """Plots values on 2-D unstructured mesh."""
    import pylab

    r = (max(x) - min(x))/(max(y) - min(y))
    rx = r/np.sqrt(1 + r**2)
    ry = 1/np.sqrt(1 + r**2)

    f = pylab.figure(figsize=(10*rx, 10*ry))
    p = pylab.tricontourf(x, y, z, 125)
    pylab.axis('image')
    return f, p


def meshsmooth(x, z, v, span, nx, nz):
    """Smooths values on 2-D unstructured mesh."""

    # construct rectangular grid
    xi = np.linspace(x.min(), x.max(), nx)
    zi = np.linspace(z.min(), z.max(), nz)
    xi, zi = np.meshgrid(xi, zi)
    xi = xi.flatten()
    zi = zi.flatten()

    # go from unstructured 'mesh' to rectangular 'grid'
    meshcoords = np.column_stack([x, z])
    gridcoords = np.column_stack([xi, zi])
    vi = _interp.griddata(meshcoords, v, gridcoords, 'linear')

    # smooth
    xi = np.reshape(xi, (nz, nx))
    zi = np.reshape(zi, (nz, nx))
    vi = np.reshape(vi, (nz, nx))
    vs = gridsmooth(vi, span)

    # back to unstructured mesh
    xi = xi.flatten()
    zi = zi.flatten()
    gridcoords = np.column_stack([xi, zi])
    vs = vs.flatten()
    vs = _interp.griddata(gridcoords, vs, meshcoords, 'linear')

    return vs


def gauss2(X, Y, mu, sigma):
    """Evaluates Gaussian over points of X,Y."""
    # evaluates Gaussian over X,Y
    D = sigma[0, 0]*sigma[1, 1] - sigma[0, 1]*sigma[1, 0]
    B = np.linalg.inv(sigma)
    X = X - mu[0]
    Y = Y - mu[1]
    Z = B[0, 0]*X**2. + B[0, 1]*X*Y + B[1, 0]*X*Y + B[1, 1]*Y**2.
    Z = np.exp(-0.5*Z)
    Z *= (2.*np.pi*np.sqrt(D))**(-1.)
    return Z


def sortrows(a, return_index=False, return_inverse=False):
    """Sorts rows of numpy array."""
    si = np.lexsort(a.T)
    if return_inverse:
        sj = np.argsort(si)

    if return_index and return_inverse:
        return a[si], si, sj
    elif return_index:
        return a[si], si
    elif return_inverse:
        return a[si], sj
    else:
        return a[si]


def uniquerows(a, sort_array=False, return_index=False):
    """Finds unique rows of numpy array."""
    if sort_array:
        if return_index:
            sa, si = sortrows(a, return_index=True)
        else:
            sa = sortrows(a)
    else:
        sa, sj = sortrows(a, return_inverse=True)
    ui = np.ones(len(sa), 'bool')
    ui[1:] = (np.diff(sa, axis=0) != 0).any(axis=1)

    if sort_array:
        ua = sa[ui]
        if return_index:
            ui = si[ui]
    else:
        ua = a[ui[sj]]
        if return_index:
            ui = np.array(range(len(ui)))[ui[sj]]

    if return_index:
        return ua, ui
    else:
        return ua


def loadnpy(filename):
    """Loads numpy binary file."""
    return np.load(filename)


def savenpy(filename, v):
    """Saves numpy binary file."""
    np.save(filename, v)
    os.rename(filename + '.npy', filename)
