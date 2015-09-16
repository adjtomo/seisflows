
import numpy as np


### user supplied preconditioners

class diagonal(object):
    """ User supplied diagonal preconditioner

        Rescales model parameters based on user supplied weights
    """
    def __init__(self):
        """ Loads any required dependencies
        """
        from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths
        import solver

        PAR = SeisflowsParameters()
        PATH = SeisflowsPaths()

        self.path = PATH.PRECOND
        self.load = solver.load
        self.merge = solver.merge


    def __call__(self, q):
        """ Applies preconditioner to given vector
        """
        p = self.merge(self.load(self.path))
        return p*q


### experimental

def fix(A):
    nrow = A.shape[0]
    ncol = A.shape[1]
    for i in range(ncol):
       if A[i,i] < 0.:
           A[:,i] *= -1
    return A


class pca(object):
    """ PCA preconditioner

        Equivalent to a change of material parameters, with choice of
        new parameters based on principle component analysis 
    """
    def __init__(self):
        """ Loads any required dependencies
        """
        import solver

        # solver methods
        self.merge = solver.merge
        self.split = solver.split

        # solver properties
        self.nproc = solver.mesh.nproc
        self.ngll = solver.mesh.ngll
        self.parameters = solver.parameters
 

    def __call__(self, q):
        """ Applies preconditioner to given vector
        """
        old = self.split(q)
        nn = len(self.parameters)

        # compute covariance
        cov = np.zeros((nn,nn))
        for ii,ikey in enumerate(self.parameters):
            for jj,jkey in enumerate(self.parameters):
                for iproc in range(self.nproc):
                    cov[ii,jj] += np.dot(old[ikey][iproc], old[jkey][iproc])

        inv = self.invert(cov)

        if True:
            print 'cov:'
            print cov
            print 'inv:'
            print inv

        # apply preconditioner
        new = {}
        for ii,ikey in enumerate(self.parameters):
            # initialize with zeros
            new[ikey] = []
            for iproc in range(self.nproc):
                ngll = self.ngll[iproc]
                new[ikey] += [np.zeros(ngll)]

            for iproc in range(self.nproc):
                for jj,jkey in enumerate(self.parameters):
                            new[ikey][iproc] += inv[ii,jj]*old[jkey][iproc]

        return self.merge(new)


    def invert(self, cov):
        inv = np.linalg.inv(cov)
        inv /= (np.trace(inv)/len(self.parameters))
        return inv
        

class pca2(pca):
    def invert(self, C):

        d,E = np.linalg.eig(C)
        E = fix(E)
        D = np.diag(d)
        W = np.dot(D**0.1, E)

        return np.linalg.inv(np.dot(W.T, W))


