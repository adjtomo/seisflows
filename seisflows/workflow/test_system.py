
import sys
import time

from seisflows.config import ParameterError

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']


class test_system:
    """ Tests system interface
    """

    def check(self):
        if 'NTASK' not in PAR:
            raise Exception

        if 'NPROC' not in PAR:
            setattr(PAR,'NPROC',1)

        if 'VERBOSE' not in PAR:
            setattr(PAR,'VERBOSE',0)


    def main(self):
        system.run('workflow', 'hello', 
            hosts='all')

        system.run('workflow', 'hello',  
            hosts='head', 
            msg=None)

        print ''


    def hello(self, msg='Hello'):
        """ Sends hello message from compute node
        """
        if not msg:
            return

        time.sleep(1)
        print '%s from %d' % (msg, system.taskid()+1)
        print ''

