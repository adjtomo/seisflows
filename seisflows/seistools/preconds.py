

class diagonal(object):
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


class pca(object):
    def __init__(solver=None):
        """ Loads any required dependencies
        """
        # solver methods
        self.merge = solver.merge
        self.split = solver.split

        # solver properties
        self.nproc = solver.mesh.nproc
        self.parameters = solver.parameters
 

    def __call__(q):
        """ Applies preconditioner to given vector
        """
        r = self.split(q)

        # compute covariance
        c = {}
        for key1 in self.parameters:
            c[key1] = {}
            for key2 in self.parameters:
                c[key1][key2] = np.dot(r[key1], r[key2])

        # diagonalize
        a = eig(c)
        b = inv(a)

        # apply preconditioner
        s = {}
        for key1 in self.parameters:
            for key2 in self.parameters:
                for iproc in range(self.nproc):
                    s[key1][iproc] +=  [d[key2]*r[key][iproc]]
        return s

