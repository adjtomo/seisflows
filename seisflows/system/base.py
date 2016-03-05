
from os.path import abspath, join

from seisflows.tools.config import SeisflowsObjects, SeisflowsParameters, SeisflowsPaths

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


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
        for obj in [SeisflowsParameters(), SeisflowsPaths(), SeisflowsObjects()]:
           obj.save(PATH.OUTPUT)

