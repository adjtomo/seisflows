#!/bin/env python

from seisflows.tools import unix
from seisflows.tools.config import loadclass, loadvars, ConfigObj, ParameterObj, Null

OBJ = ConfigObj('SeisflowsObjects')
PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

register = OBJ.register

null = Null()


class run:
  """ Tests system interface.
  """

  def check(self):
    if 'NTASK' not in PAR:
        raise Exception

    if 'NPROC' not in PAR:
        setattr(PAR,'NPROC',1)

    if 'VERBOSE' not in PAR:
        setattr(PAR,'VERBOSE',0)


  def main(self):
    system.run('workflow', 'hello', hosts='all')
    print ''


  def hello(self):
    """ Sends hello message from compute node
    """
    import time
    time.sleep(1)
    print 'Hello from', system.getnode()
    print ''


# run test
if __name__ == '__main__':

    PAR.update(loadvars('parameters','.'))
    PATH.update(loadvars('paths','.'))

    workflow = run()
    workflow.check()
    register('workflow',workflow)

    system = loadclass('system',PAR.SYSTEM)()
    system.check()
    register('system',system)

    register('preprocess',null)
    register('solver',null)
    register('postprocess',null)
    register('optimize',null)

    system.submit(workflow)

