import os

import numpy as np
import scipy.signal as _signal
import scipy.interpolate as _interp

from seisflows.tools.math import gauss2


def sortrows(a, return_index=False, return_inverse=False):
    """ Sorts rows of numpy array
    """
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
    """ Finds unique rows of numpy array
    """
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


def stack(*args):
    return np.column_stack(args)


### array input/output

def loadnpy(filename):
    """Loads numpy binary file."""
    return np.load(filename)


def savenpy(filename, v):
    """Saves numpy binary file."""
    np.save(filename, v)
    os.rename(filename + '.npy', filename)




# in the function and variable names, we use 'grid' to describe a set of
# structured coordinates, and 'mesh' to describe a set of unstructured 
# coordinates


def gridsmooth(Z, span):
    """ Smooths values on 2D rectangular grid
    """
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


def meshsmooth(v, mesh, span):
    """ Smooths values on 2D unstructured mesh
    """
    V, grid = mesh2grid(v, mesh)
    nz, nx = V.shape
    W = np.ones((nz, nx))

    # maks nans
    inan = np.isnan(V)
    if np.any(inan):
        V[inan] = 0.
        W[inan] = 0.

    # apply smoother
    V = gridsmooth(V, span)
    W = gridsmooth(W, span)
    V = V/W

    if np.any(inan):
        V[inan] = 0.

    vs = grid2mesh(V, grid, mesh)
    return vs


def mesh2grid(v, mesh):
    """ Interpolates from an unstructured coordinates (mesh) to a structured 
        coordinates (grid)
    """
    x = mesh[:,0]
    z = mesh[:,1]
    lx = x.max() - x.min()
    lz = z.max() - z.min()
    nn = v.size

    nx = np.around(np.sqrt(nn*lx/lz))
    nz = np.around(np.sqrt(nn*lz/lx))
    dx = lx/nx
    dz = lz/nz

    # construct structured grid
    x = np.linspace(x.min(), x.max(), nx)
    z = np.linspace(z.min(), z.max(), nz)
    X, Z = np.meshgrid(x, z)
    grid = stack(X.flatten(), Z.flatten())

    # interpolate to structured grid
    V = _interp.griddata(mesh, v, grid, 'linear')

    # workaround edge issues
    if np.any(np.isnan(V)):
        W = _interp.griddata(mesh, v, grid, 'nearest')
        for i in np.where(np.isnan(V)):
            V[i] = W[i]

    V = np.reshape(V, (nz, nx))
    return V, grid


def grid2mesh(V, grid, mesh):
    """ Interpolates from structured coordinates (grid) to unstructured 
        coordinates (mesh)
    """
    return _interp.griddata(grid, V.flatten(), mesh, 'linear')

