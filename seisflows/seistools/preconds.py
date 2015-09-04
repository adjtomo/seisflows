
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

        # diagonalize
        #eigval,eigvec = np.linalg.eig(cov)
        #precond = np.dot(eigvec.T, eigvec)
        #inv = np.linalg.inv(precond)

        inv = self.invert(cov)

        if True:
            print 'cov:', cov
            print 'inv:', inv

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
    def invert(self, cov):
        inv = np.linalg.inv(cov)
        inv /= (np.trace(inv)/len(self.parameters))
        return np.diag(np.diag(inv))


class pca3(pca):
    def invert(self, cov):
        inv = np.linalg.inv(cov)
        inv /= (np.trace(inv)/len(self.parameters))
        return np.diag(np.diag(inv)**0.5)


class pca4(pca):
    def invert(self, cov):
        nn = len(self.parameters)
        cov += np.trace(cov)/nn * np.eye(nn)
        inv = np.linalg.inv(cov)
        inv /= (np.trace(inv)/len(self.parameters))
        return inv


