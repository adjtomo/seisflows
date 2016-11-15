
from os.path import join
import sys
import numpy as np

from seisflows.tools import msg
from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import divides, exists
from seisflows.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver
import optimize
import preprocess
import postprocess


class inversion(object):
    """ Waveform inversion base class.

      Peforms iterative nonlinear inversion and provides a base class on top
      of which specialized strategies can be implemented.

      To allow customization, the inversion workflow is divided into generic 
      methods such as 'initialize', 'finalize', 'evaluate_function', 
      'evaluate_gradient', which can be easily overloaded.

      Calls to forward and adjoint solvers are abstracted through the 'solver'
      interface so that various forward modeling packages can be used
      interchangeably.

      Commands for running in serial or parallel on a workstation or cluster
      are abstracted through the 'system' interface.
    """

    def check(self):
        """ Checks parameters and paths
        """

        # check parameters
        if 'BEGIN' not in PAR:
            raise ParameterError(PAR, 'BEGIN')

        if 'END' not in PAR:
            raise ParameterError(PAR, 'END')

        if 'VERBOSE' not in PAR:
            setattr(PAR, 'VERBOSE', 1)

        # scratch paths
        if 'SCRATCH' not in PATH:
            raise ParameterError(PATH, 'SCRATCH')

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        if 'FUNC' not in PATH:
            setattr(PATH, 'FUNC', join(PATH.SCRATCH, 'evalfunc'))

        if 'GRAD' not in PATH:
            setattr(PATH, 'GRAD', join(PATH.SCRATCH, 'evalgrad'))

        if 'HESS' not in PATH:
            setattr(PATH, 'HESS', join(PATH.SCRATCH, 'evalhess'))

        if 'OPTIMIZE' not in PATH:
            setattr(PATH, 'OPTIMIZE', join(PATH.SCRATCH, 'optimize'))

        # input paths
        if 'DATA' not in PATH:
            setattr(PATH, 'DATA', None)

        if 'MODEL_INIT' not in PATH:
            raise ParameterError(PATH, 'MODEL_INIT')

        # output paths
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

        # parameter assertions
        assert 1 <= PAR.BEGIN <= PAR.END

        # path assertions
        if not exists(PATH.DATA):
            assert 'MODEL_TRUE' in PATH
            assert exists(PATH.MODEL_TRUE)

        if not exists(PATH.MODEL_INIT):
            raise Exception()


    def main(self):
        """ Carries out seismic inversion
        """
        optimize.iter = PAR.BEGIN
        self.setup()
        print ''

        while optimize.iter <= PAR.END:
            print "Starting iteration", optimize.iter
            self.initialize()

            print "Computing gradient"
            self.evaluate_gradient()

            print "Computing search direction"
            self.compute_direction()

            print "Computing step length"
            self.line_search()

            self.finalize()
            self.clean()

            optimize.iter += 1
            print ''


    def setup(self):
        """ Lays groundwork for inversion
        """
        # clean scratch directories
        if PAR.BEGIN == 1:
            unix.rm(PATH.SCRATCH)
            unix.mkdir(PATH.SCRATCH)

            preprocess.setup()
            postprocess.setup()
            optimize.setup()

        if PATH.DATA:
            print 'Copying data' 
        else:
            print 'Generating data' 

        system.run('solver', 'setup', 
                   hosts='all')


    def initialize(self):
        """ Prepares for next model update iteration
        """
        self.write_model(path=PATH.GRAD, suffix='new')

        print 'Generating synthetics'
        system.run('solver', 'eval_func',
                   hosts='all',
                   path=PATH.GRAD)

        self.sum_residuals(path=PATH.GRAD, suffix='new')


    def compute_direction(self):
        """ Computes search direction
        """
        optimize.compute_direction()


    def line_search(self):
        """ Conducts line search in given search direction
        """
        optimize.initialize_search()

        while True:
            self.iterate_search()

            if optimize.isdone:
                optimize.finalize_search()
                break
            elif optimize.step_count < PAR.STEPMAX:
                optimize.compute_step()
                continue
            else:
                retry = optimize.retry_status()
                if retry:
                    print ' Line search failed\n\n Retrying...'
                    optimize.restart()
                    self.line_search()
                    break
                else:
                    print ' Line search failed\n\n Aborting...'
                    sys.exit(-1)


    def iterate_search(self):
        """ First, calls self.evaluate_function, which carries out a forward 
          simulation given the current trial model. Then calls
          optimize.update_status, which maintains search history and checks
          stopping conditions.
        """
        if PAR.VERBOSE > 0:
            print " trial step", optimize.step_count+1

        self.evaluate_function()
        optimize.update_status()


    def evaluate_function(self):
        """ Performs forward simulation to evaluate objective function
        """
        self.write_model(path=PATH.FUNC, suffix='try')

        system.run('solver', 'eval_func',
                   hosts='all',
                   path=PATH.FUNC)

        self.sum_residuals(path=PATH.FUNC, suffix='try')


    def evaluate_gradient(self):
        """ Performs adjoint simulation to evaluate gradient
        """
        system.run('solver', 'eval_grad',
                   hosts='all',
                   path=PATH.GRAD,
                   export_traces=divides(optimize.iter, PAR.SAVETRACES))

        postprocess.write_gradient(
            path=PATH.GRAD)

        src = join(PATH.GRAD, 'gradient')
        dst = join(PATH.OPTIMIZE, 'g_new')
        savenpy(dst, solver.merge(solver.load(src, suffix='_kernel')))


    def finalize(self):
        """ Saves results from current model update iteration
        """
        system.checkpoint()

        if divides(optimize.iter, PAR.SAVEMODEL):
            self.save_model()

        if divides(optimize.iter, PAR.SAVEGRADIENT):
            self.save_gradient()

        if divides(optimize.iter, PAR.SAVEKERNELS):
            self.save_kernels()

        if divides(optimize.iter, PAR.SAVETRACES):
            self.save_traces()

        if divides(optimize.iter, PAR.SAVERESIDUALS):
            self.save_residuals()


    def clean(self):
        """ Cleans directories in which function and gradient evaluations were
          carried out
        """
        unix.rm(PATH.GRAD)
        unix.rm(PATH.FUNC)
        unix.mkdir(PATH.GRAD)
        unix.mkdir(PATH.FUNC)


    def write_model(self, path='', suffix=''):
        """ Writes model in format used by solver
        """
        unix.mkdir(path)
        src = PATH.OPTIMIZE +'/'+ 'm_' + suffix
        dst = path +'/'+ 'model'
        parts = solver.split(loadnpy(src))
        solver.save(dst, parts)


    def sum_residuals(self, path='', suffix=''):
        """ Returns sum of squares of residuals
        """
        src = path +'/'+ 'residuals'
        dst = PATH.OPTIMIZE +'/'+ 'f_' + suffix
        residuals = []
        for file in unix.ls(src):
            fromfile = np.loadtxt(src +'/'+ file)
            residuals.append(fromfile**2.)
        np.savetxt(dst, [np.sum(residuals)])


    def save_gradient(self):
        src = join(PATH.GRAD, 'gradient')
        dst = join(PATH.OUTPUT, 'gradient_%04d' % optimize.iter)
        unix.mv(src, dst)


    def save_model(self):
        src = PATH.OPTIMIZE +'/'+ 'm_new'
        dst = join(PATH.OUTPUT, 'model_%04d' % optimize.iter)
        solver.save(dst, solver.split(loadnpy(src)))


    def save_kernels(self):
        src = join(PATH.GRAD, 'kernels')
        dst = join(PATH.OUTPUT, 'kernels_%04d' % optimize.iter)
        unix.mv(src, dst)


    def save_traces(self):
        src = join(PATH.GRAD, 'traces')
        dst = join(PATH.OUTPUT, 'traces_%04d' % optimize.iter)
        unix.mv(src, dst)


    def save_residuals(self):
        src = join(PATH.GRAD, 'residuals')
        dst = join(PATH.OUTPUT, 'residuals_%04d' % optimize.iter)
        unix.mv(src, dst)



