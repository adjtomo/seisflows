#!/usr/bin/python

from seisflows.tools import unix
from seisflows.tools.configure import getclass, ParameterObject

import parameters
import paths

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')

system = getclass('system',PAR.SYSTEM)()


class run:
  """ Tests system interface.
  """

  def main(self):
    system.run(self.hello, hosts='all')
    print ''


  def hello(self):
    # send hello message from compute node
    import time
    time.sleep(1)
    print 'Hello from', system.getnode()
    print ''


# run test
if __name__ == '__main__':

    system.submit(run)

