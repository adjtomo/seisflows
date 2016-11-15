
import sys

from seisflows.config import save

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class base(object):
    """ Abstract base class
    """

    def check(self):
        raise NotImplementedError('Must be implemented by subclass.')

    def submit(self):
        raise NotImplementedError('Must be implemented by subclass.')

    def run(self):
        raise NotImplementedError('Must be implemented by subclass.')

    def getnode(self):
        raise NotImplementedError('Must be implemented by subclass.')

    def checkpoint(self):
        save()

