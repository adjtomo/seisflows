import numpy as np

from seisflows.seistools import adjoint, misfit

from seisflows.tools import unix
from seisflows.tools.code import Struct
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, custom_import

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver


class double_difference(custom_import('preprocess', 'legacy')):
    """ Data preprocessing class
    """

    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        super(DoubleDifference, self).check()


    def prepare_eval_grad(self, path='.'):
        """ Prepares solver for gradient evaluation by writing residuals and
          adjoint traces
        """
        unix.cd(path)

        d, h = self.load(prefix='traces/obs/')
        s, _ = self.load(prefix='traces/syn/')

        d = self.apply(self.process_traces, [d], [h])
        s = self.apply(self.process_traces, [s], [h])

        r = self.apply(self.write_residuals, [s, d], [h], inplace=False)

        s = self.apply(self.generate_adjoint_traces, [s, d, r], [h])
        self.save(s, h, prefix='traces/adj/')


    def write_residuals(self, s, d, h):
        """ Computes residuals from observations and synthetics
        """
        nr = h.nr
        nt = h.nt
        dt = h.dt

        r = np.zeros((nr, nr))
        for ir in range(nr):
            for jr in range(nr):

                if ir < jr:
                    r[ir, jr] = (misfit.wtime(s[:, ir], s[:, jr], nt, dt) -
                                 misfit.wtime(d[:, ir], d[:, jr], nt, dt))

                elif ir > jr:
                    r[ir, jr] = -r[ir, jr]

                else:
                    r[ir, jr] = 0

        # write residuals
        np.savetxt('residuals', np.sqrt(np.sum(r*r, 0)))

        return np.array(r)


    def generate_adjoint_traces(self, s, d, r, h):
        """ Computes adjoint traces from observed and synthetic traces
        """
        nr = h.nr
        dt = h.dt

        for ir in range(nr):
            nrm = sum((s[:, ir]*s[:, ir])*dt)
            fit = np.sum(r[ir, :])

            s[1:-1, ir] = (s[2:, ir] - s[0:-2, ir])/(2.*dt)
            s[0, ir] = 0.
            s[-1, ir] = 0.
            s[:, ir] *= fit/nrm

        return s


