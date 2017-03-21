
import sys
import numpy as np

from seisflows.config import custom_import, ParameterError

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class steepest_descent(custom_import('optimize', 'base')):
    """ Steepest descent method
    """
    restarted = False

    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        # line search algorithm
        if 'LINESEARCH' not in PAR:
            setattr(PAR, 'LINESEARCH', 'Bracket')

        super(steepest_descent, self).check()


    def setup(self):
        super(steepest_descent, self).setup()


    def compute_direction(self):
        g_new = self.load('g_new')

        precond=self.precond()

        if precond:
            p_new = -precond(g_new)
        else:
            p_new = -g_new

        self.save('p_new', p_new)
        savetxt('s_new', self.dot(g_new, p_new))

        return p_new


    def restart(self):
        pass

