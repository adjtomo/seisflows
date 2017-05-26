
import sys
import time

from random import random
from seisflows.config import ParameterError

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']


class test_fault_tolerance:
    """ Tests system interface
    """

    def check(self):
        if 'NTASK' not in PAR:
            raise Exception

        if 'NPROC' not in PAR:
            setattr(PAR, 'NPROC', 1)

        if 'VERBOSE' not in PAR:
            setattr(PAR, 'VERBOSE', 0)

        if 'FAILRATE' not in PAR:
            setattr(PAR, 'FAILRATE', 0.667)

        # assertions
        assert 0. <= PAR.FAILRATE < 1.


    def main(self):
        system.run('workflow', 'hello', 
            hosts='all')

        print ''


    def main(self):
        system.run('workflow', 'hello',
            hosts='all',
            msg='Hello from %d')

        print ''


    def hello(self, msg='Hello from %d'):
        """ Prints hello message
        """
        time.sleep(1)

        if random() < FAILRATE:
            print 'task failed...'
            sys.exit(-1)

        try:
            print msg % system.taskid()+1 
        except:
            print msg

        print ''

