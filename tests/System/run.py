#!/usr/bin/python

from seisflows.tools import unix
from seisflows.tools.configtools import getclass, ParameterObject

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')

system = getclass('system',PAR.SYSTEM)()


class run:
  """ Tests system interface.
  """

  def __init__(self):

    # check parameters
    if 'NTASK' not in PAR:
        raise Exception

    if 'NPROC' not in PAR:
        setattr(PAR,'NPROC',1)

    if PAR.VERBOSE != 0:
        setattr(PAR,'VERBOSE',0)


  def main(self):
    system.run(self.hello, hosts='all')
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

    system.submit(run)

