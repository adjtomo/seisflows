
from os.path import join
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import divides, exists
from seisflows.tools.config import ParameterError, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import system
import solver
import optimize
import preprocess
import postprocess


class inversion(object):
    """ Seismic inversion base class.

      Computes iterative model updates using any of a variety of possible 
      settings, and provides a base class on top of which custom inversion 
      strategies can be implemented.

      To allow overriding, separate methods are provided for each step in the
      inversion, including 'evaluate_function', 'evaluate_gradient',
      'compute_direction', and 'line_search' as well as generic 'initialize'
      and 'finalize' methods.

      Calls to forward and adjoint solvers are abstracted through the 'solver'
      interface so that various forward modeling packages can be used
      interchangeably.

      Commands for launching serial or parallel jobs on a workstation or cluster
      are abstracted through the 'system' interface so that the code can be
      easily ported.

      For assistance using this package, please view online documentation, 
      browse comments, or email rmodrak -at- princeton -dot- edu
    """

    def check(self):
        """ Checks parameters, paths, and dependencies
        """

        # check parameters
        if 'BEGIN' not in PAR:
            raise ParameterError(PAR, 'BEGIN')

        if 'END' not in PAR:
            raise ParameterError(PAR, 'END')

        if 'VERBOSE' not in PAR:
            setattr(PAR, 'VERBOSE', 1)

        # check paths
        if 'GLOBAL' not in PATH:
            raise ParameterError(PATH, 'GLOBAL')

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        if 'FUNC' not in PATH:
            setattr(PATH, 'FUNC', join(PATH.GLOBAL, 'evalfunc'))

        if 'GRAD' not in PATH:
            setattr(PATH, 'GRAD', join(PATH.GLOBAL, 'evalgrad'))

        if 'HESS' not in PATH:
            setattr(PATH, 'HESS', join(PATH.GLOBAL, 'evalhess'))

        if 'OPTIMIZE' not in PATH:
            setattr(PATH, 'OPTIMIZE', join(PATH.GLOBAL, 'optimize'))

        # input settings
        if 'DATA' not in PATH:
            setattr(PATH, 'DATA', None)

        if not exists(PATH.DATA):
            assert 'MODEL_TRUE' in PATH

        if 'MODEL_INIT' not in PATH:
            raise ParameterError(PATH, 'MODEL_INIT')

        # output settings
        if 'OUTPUT' not in PATH:
            raise ParameterError(PATH, 'OUTPUT')

        if 'SAVEMODEL' not in PAR:
            setattr(PAR, 'SAVEMODEL', 1)

        if 'SAVEGRADIENT' not in PAR:
            setattr(PAR, 'SAVEGRADIENT', 0)

        if 'SAVEKERNELS' not in PAR:
            setattr(PAR, 'SAVEKERNELS', 0)

        if 'SAVETRACES' not in PAR:
            setattr(PAR, 'SAVETRACES', 0)

        if 'SAVERESIDUALS' not in PAR:
            setattr(PAR, 'SAVERESIDUALS', 0)


    def main(self):
        """ Carries out seismic inversion
        """
        self.setup()

        for self.iter in range(PAR.BEGIN, PAR.END+1):
            optimize.iter = self.iter

            print "Starting iteration", self.iter
            self.initialize()

            print "Computing search direction"
            self.compute_direction()

            print "Computing step length"
            self.line_search()

            self.finalize()
            print ''

            if self.isdone:
                return


    def setup(self):
        """ Lays groundwork for inversion
        """
        # clean scratch directories
        if PAR.BEGIN == 1:
            unix.rm(PATH.GLOBAL)
            unix.mkdir(PATH.GLOBAL)

        # set up optimization
        optimize.setup()

        # set up pre- and post-processing
        preprocess.setup()
        postprocess.setup()

        # set up solver
        if PAR.BEGIN == 1:
            system.run('solver', 'setup',
                       hosts='all')
            return

        if PATH.LOCAL:
             system.run('solver', 'setup',
                        hosts='all')
 

    def initialize(self):
        """ Prepares for next model update iteration
        """
        isready = self.solver_status()
        if not isready:
            print 'Generating synthetics'

            self.prepare_model(path=PATH.GRAD, suffix='new')

            system.run('solver', 'eval_func',
                       hosts='all',
                       path=PATH.GRAD)

            self.sum_residuals(path=PATH.GRAD, suffix='new')


    def compute_direction(self):
        """ Computes search direction
        """
        self.evaluate_gradient()
        optimize.compute_direction()


    def line_search(self):
        """ Conducts line search in given search direction
        """
        optimize.initialize_search()

        for optimize.step in range(1, PAR.SRCHMAX+1):
            isdone = self.search_status()

            if isdone == 1:
                optimize.finalize_search()
                break
            elif isdone == 0:
                optimize.compute_step()
                continue
            elif isdone == -1:
                self.isdone = -1
                print ' line search failed'


    def search_status(self):
        """ Evaluates trial model and checks line search status

          Evaluates trial model by calling evalute_function, which results in
          a forward simulation. After that, queries line search status and keeps
          track of which trial model is the best so far.
        """
        if PAR.VERBOSE:
            print " trial step", optimize.step
        self.evaluate_function()

        isdone, isbest = optimize.search_status()
        if not PATH.LOCAL:
            if isbest and isdone:
                unix.rm(PATH.SOLVER + '_best')
                unix.mv(PATH.SOLVER, PATH.SOLVER + '_best')
            elif isbest:
                unix.rm(PATH.SOLVER + '_best')
                unix.cp(PATH.SOLVER, PATH.SOLVER + '_best')

        return isdone


    def evaluate_function(self):
        """ Calls forward solver and writes misfit
        """
        self.prepare_model(path=PATH.FUNC, suffix='try')

        system.run('solver', 'eval_func',
                   hosts='all',
                   path=PATH.FUNC)

        self.sum_residuals(path=PATH.FUNC, suffix='try')


    def evaluate_gradient(self):
        """ Calls adjoint solver and runs process_kernels
        """
        # adjoint simulation
        system.run('solver', 'eval_grad',
                   hosts='all',
                   path=PATH.GRAD,
                   export_traces=divides(self.iter, PAR.SAVETRACES))

        postprocess.process_kernels(
            path=PATH.GRAD)


    def finalize(self):
        """ Saves results from most recent model update iteration
        """
        if divides(self.iter, PAR.SAVEMODEL):
            self.save_model()

        if divides(self.iter, PAR.SAVEGRADIENT):
            self.save_gradient()

        if divides(self.iter, PAR.SAVEKERNELS):
            self.save_kernels()

        if divides(self.iter, PAR.SAVETRACES):
            self.save_traces()

        if divides(self.iter, PAR.SAVERESIDUALS):
            self.save_residuals()

        # clean up directories for next iteration
        if not PATH.LOCAL:
            unix.rm(PATH.GRAD)
            unix.mv(PATH.FUNC, PATH.GRAD)
            unix.mkdir(PATH.FUNC)

            unix.rm(PATH.SOLVER)
            unix.mv(PATH.SOLVER + '_best', PATH.SOLVER)

        else:
            unix.rm(PATH.GRAD)
            unix.rm(PATH.FUNC)
            unix.mkdir(PATH.GRAD)
            unix.mkdir(PATH.FUNC)

        self.isdone = False


    ### utility functions

    def prepare_model(self, path='', suffix=''):
        """ Writes model in format used by solver
        """
        unix.mkdir(path)
        src = PATH.OPTIMIZE +'/'+ 'm_' + suffix
        dst = path +'/'+ 'model'
        parts = solver.split(loadnpy(src))
        solver.save(dst, parts)


    def sum_residuals(self, path='', suffix=''):
        """ Sums residuals to obtain misfit function value
        """
        src = path +'/'+ 'residuals'
        dst = PATH.OPTIMIZE +'/'+ 'f_' + suffix
        residuals = []
        for file in unix.ls(src):
            fromfile = np.loadtxt(src +'/'+ file)
            residuals.append(fromfile**2.)
        np.savetxt(dst, [np.sum(residuals)])


    def solver_status(self):
        """ Decides if solver prerequisites are in place
        """
        if self.iter <= 1:
            isready = False
        elif PATH.LOCAL:
            isready = False
        else:
            isready = True
        return isready


    def save_gradient(self):
        src = join(PATH.GRAD, 'gradient')
        dst = join(PATH.OUTPUT, 'gradient_%04d' % self.iter)
        unix.mv(src, dst)


    def save_model(self):
        src = PATH.OPTIMIZE +'/'+ 'm_new'
        dst = join(PATH.OUTPUT, 'model_%04d' % self.iter)
        solver.save(dst, solver.split(loadnpy(src)))


    def save_kernels(self):
        src = join(PATH.GRAD, 'kernels')
        dst = join(PATH.OUTPUT, 'kernels_%04d' % self.iter)
        unix.mv(src, dst)


    def save_traces(self):
        src = join(PATH.GRAD, 'traces')
        dst = join(PATH.OUTPUT, 'traces_%04d' % self.iter)
        unix.mv(src, dst)


    def save_residuals(self):
        src = join(PATH.GRAD, 'residuals')
        dst = join(PATH.OUTPUT, 'residuals_%04d' % self.iter)
        unix.mv(src, dst)

