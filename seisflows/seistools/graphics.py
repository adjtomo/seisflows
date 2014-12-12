
import numpy as _np

from seisflows.tools.array import gridplot


def splot(s,h,normalize=True,nxmax=1000,nymax=1000):
    """Plots seismic record section"""
    ny = s.shape[0]
    nx = s.shape[1]
    iy = range(ny)
    ix = range(nx)

    if ny > nymax:
        iy = _ceil(ny,nymax)
        iy = range(0,ny,iy)

    if nx > nxmax:
        ix = _ceil(nx,nxmax)
        ix = range(0,nx,ix)

    for ir in range(h.nr):
        if normalize:
            w = max(abs(s[:,ir]))
            if w > 0:
                s[:,ir] = s[:,ir]/w
            else:
                break

    gridplot(_np.flipud(s[_np.ix_(iy,ix)]))


def wplot(w):
    """Plots seismic waveform"""

    import pylab

    pylab.figure(figsize=(7.5,3))
    line = pylab.plot(w[:,0],w[:,1])

    pylab.setp(line,color='k',linestyle='-')
    pylab.setp(line,color='k',linestyle='', marker='.',markersize=0.5)
    pylab.axis('tight')
    pylab.title(file)

    pylab.show()



### utility functions

def _ceil(x,y):
    z = float(x)/float(y)
    return int(_np.ceil(z))
