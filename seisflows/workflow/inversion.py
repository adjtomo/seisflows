
import numpy as np

from seisflows.tools import unix
from seisflows.tools.arraytools import loadnpy, savenpy
from seisflows.tools.codetools import divides, glob, irange, join
from seisflows.tools.configure import getclass, ParameterObject

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')

system = getclass('system',PAR.SYSTEM)()
solver = getclass('solver',PAR.SOLVER)()



class inversion(object):
  """ Seismic inversion base class.

    Computes iterative model updates using a variety of possible settings, and
    provides a base class on top of which custom inversion strategies can be
    implemented.

    To allow overriding, separate methods are provided for each step in the
    inversion, including 'evaluate_function', 'evaluate_gradient', 
    'compute_direction', and 'line_search' as well as generic 'initialize'
    and 'finalize' methods. 

    Calls to forward and adjoint solvers are abstracted through the 'solver'
    interface so that various forward modeling packages can be used 
    interchangably.

    Commands for launching serial or parallel jobs on a workstation or cluster
    are abstracted through the 'system' interface so that the code can be
    easily ported.

    For assistance using this package, please browse comments or email
    rmodrak -at- princeton -dot- edu
  """


  def __init__(self):
    """ Class constructor
    """
    self.iter = 0

    # check user supplied parameters
    if 'BEGIN' not in PAR:
	raise Exception

    if 'END' not in PAR:
	raise Exception

    if 'SAVEMODELS' not in PAR:
	setattr(PAR,'SAVEMODELS',1)

    if 'SAVEKERNELS' not in PAR:
        setattr(PAR,'SAVEKERNELS',0)

    if 'SAVETRACES' not in PAR:
        setattr(PAR,'SAVETRACES',0)

    self.configure_paths()


  def main(self):
    """ Carries out seismic inversion
    """
    print 'Setting up...'
    print 'Preparing data...\n'
    self.setup()

    for self.iter in irange(PAR.BEGIN,PAR.END):
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
    self.optimize = getclass('optimize',PAR.OPTIMIZE)()
    self.postprocess = getclass('postprocess',PAR.POSTPROCESS)()

    # prepare solver
    isready = self.solver_status()
    if not isready:
      system.run( solver.setup,
          hosts='all' )

    if PAR.BEGIN==1:
      # prepare starting model
      parts = solver.load(PATH.OUTPUT+'/'+'model_init')
      src = PATH.OUTPUT+'/'+'model_init'
      dst = PATH.OPTIMIZE+'/'+'m_new'
      savenpy(dst,solver.merge(parts))


  def initialize(self):
    """ Prepares for next model update iteration
    """
    isready = self.solver_status()
    if not isready:
      print 'Generating synthetics'

      self.prepare_model(path=PATH.GRAD, suffix='new')

      # forward simulation
      system.run( solver.evaluate_func,
          hosts='all',
          path=PATH.GRAD )

      self.write_misfit(path=PATH.GRAD, suffix='new')


  def compute_direction(self):
    """ Computes search direction
    """
    if PAR.SCHEME in ['sd','cg','qn']:
      self.evaluate_gradient()
      self.optimize.compute_direction()

    elif PAR.SCHEME in ['gn','tn']:
      self.optimize.initialize_newton()
      for self.ilcg in irange(1,PAR.LCGMAX):
	self.apply_hessian()
	isdone = self.optimize.solve_newton()
	if isdone:
	  break


  def line_search(self):
    """ Conducts line search in given search direction
    """
    self.optimize.initialize_search()

    for self.step in irange(1,PAR.SRCHMAX):
      isdone = self.search_status()

      if isdone==1:
        self.optimize.finalize_search()
	break
      elif isdone==0:
        self.optimize.compute_step()
        continue
      elif isdone==-1:
        self.isdone = -1
        print ' line search failed'


  def search_status(self):
    """ Checks line search status
    """
    if PAR.VERBOSE: print " trial step", self.step
    self.evaluate_function()
    isdone, isbest = self.optimize.search_status()

    if not PATH.LOCAL:
      if isbest and isdone:
        unix.mv(PATH.SCRATCH,PATH.SCRATCH+'_best')
      elif isbest:
        unix.rm(PATH.SCRATCH+'_best')
        unix.cp(PATH.SCRATCH,PATH.SCRATCH+'_best')

    return isdone


  def evaluate_function(self):
    """ Calls forward solver and writes misfit
    """
    self.prepare_model(path=PATH.FUNC, suffix='try')

    # forward simulation
    system.run( solver.evaluate_func, 
          hosts='all',
          path=PATH.FUNC )

    self.write_misfit(path=PATH.FUNC, suffix='try')


  def evaluate_gradient(self):
    """ Calls adjoint solver and runs process_kernels
    """
    # adjoint simulation
    system.run( solver.evaluate_grad, 
          hosts='all',
          path=PATH.GRAD,
          export_traces=divides(self.iter,PAR.SAVETRACES) )

    self.postprocess.process_kernels()


  def apply_hessian(self):
     """ Computes the action of the Hessian on a given vector through 
       appropriate solver call
     """
     self.prepare_model(path=PATH.HESS, suffix='lcg')

     # solver calls
     system.run( solver.action_hess, 
         hosts='all',
         path = PATH.HESS,
         hessian=PAR.SCHEME )

     self.postprocess.process_kernels(suffix='lcg')


  def finalize(self):
    """ Saves results from model update last iteration
    """
    if divides(self.iter,PAR.SAVEMODELS):
	self.save_model()

    if divides(self.iter,PAR.SAVEKERNELS):
        self.save_kernels()

    if divides(self.iter,PAR.SAVETRACES):
        self.save_traces()

    # clean up directories for next iteration
    if not PATH.LOCAL:
      unix.rm(PATH.GRAD)
      unix.mv(PATH.FUNC,PATH.GRAD)
      unix.mkdir(PATH.FUNC)

      unix.rm(PATH.SCRATCH)
      unix.mv(PATH.SCRATCH+'_best',PATH.SCRATCH)

    else:
      unix.rm(PATH.GRAD)
      unix.rm(PATH.FUNC)
      unix.mkdir(PATH.GRAD)
      unix.mkdir(PATH.FUNC)

    self.isdone = False


  def configure_paths(self):
    """ Configures directory structure
    """
    PATH.SOLVER = join(PATH.GLOBAL,'solver')
    PATH.OPTIMIZE = join(PATH.GLOBAL,'optimize')
    PATH.POSTPROCESS = join(PATH.GLOBAL,'postprocess')

    PATH.FUNC = join(PATH.SOLVER,'func')
    PATH.GRAD = join(PATH.SOLVER,'grad')
    PATH.HESS = join(PATH.SOLVER,'hess')
    PATH.MESH = join(PATH.SOLVER,'mesh')

    PATH.OUTPUT = join(PATH.SUBMIT_DIR,'output')
    unix.mkdir(PATH.OUTPUT)

    PATH.SCRATCH = join(PATH.GLOBAL,'scratch')
    if PATH.LOCAL: PATH.SCRATCH = join(PATH.LOCAL,'scatch')



  ### utility functions

  def prepare_model(self,path='',suffix=''):
      """ Writes model in format used by solver
      """
      unix.mkdir(path)
      src = PATH.OPTIMIZE+'/'+'m_'+suffix
      dst = path+'/'+'model'
      parts = solver.split(loadnpy(src))
      solver.save(dst,parts)


  def write_misfit(self,path='',suffix=''):
      """ Sums residuals to obtain misfit function value
      """
      src = path+'/'+'residuals'
      dst = PATH.OPTIMIZE+'/'+'f_'+suffix
      residuals = []
      for file in unix.ls(src):
        fromfile = np.loadtxt(src+'/'+file)
        residuals.append(fromfile**2.)
      np.savetxt(dst,[np.sum(residuals)])


  def solver_status(self):
      """ Decides if solver prerequisites are in place
      """
      if PAR.BEGIN == 1 and self.iter == 0:
        isready = False
      elif self.iter == 1:
	isready = False
      elif PATH.LOCAL:
	isready = False
      else:
	isready = True
      return isready


  def save_model(self):
      src = PATH.OPTIMIZE+'/'+'m_new'
      dst = join(PATH.OUTPUT,'model_%04d'%self.iter)
      solver.save(dst,solver.split(loadnpy(src)))


  def save_traces(self):
      src = glob(join(PATH.GRAD,'kernels/sum*'))
      dst = join(PATH.OUTPUT,'kernels_%04d'%self.iter)
      unix.mkdir(dst)
      unix.mv(src,dst)


  def save_traces(self):
      src = join(PATH.GRAD,'traces')
      dst = join(PATH.OUTPUT,'traces_%04d'%self.iter)
      unix.mv(src,dst)


  def save_residuals(self):
      src = join(PATH.GRAD,'residuals')
      dst = join(PATH.OUTPUT,'residuals_%04d'%self.iter)
      unix.mv(src,dst)


