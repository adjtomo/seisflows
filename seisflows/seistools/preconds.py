
import numpy as np


### user supplied preconditioners

class diagonal(object):
    """ User supplied diagonal preconditioner

        Rescales model parameters based on user supplied weights
    """
    def __init__(self, path=None, solver=None):
        """ Loads any required dependencies
        """
        # path to user supplied preconditioner
        assert path
        self.path = path

        # solver methods
        self.load = solver.load
        self.merge = solver.merge


    def __call__(self,q):
        """ Applies preconditioner to given vector
        """
        p = self.merge(self.load(self.path))
        return p*q


### experimental

class pca(object):
    """ PCA preconditioner

        Equivalent to a change of material parameters, with choice of
        new parameters based on principle component analysis 
    """
    def __init__(self, solver=None):
        """ Loads any required dependencies
        """
        # solver methods
        self.merge = solver.merge
        self.split = solver.split

        # solver properties
        self.nproc = solver.mesh.nproc
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

        # diagonalize
        eigval,eigvec = np.linalg.eig(cov)
        precond = np.dot(eigvec.T, eigvec)
        inv = np.linalg.inv(precond)

        # apply preconditioner
        new = {}
        for ii,ikey in enumerate(self.parameters):
            new[ikey] = []
            for jj,jkey in enumerate(self.parameters):
                    for iproc in range(self.nproc):
                        new[ikey] += [inv[ii,jj]*old[jkey][iproc]]

        return self.merge(new)


