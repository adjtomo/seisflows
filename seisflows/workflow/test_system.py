
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system


class test_system:
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

