
from os.path import join
 
from seisflows.tools.config import SeisflowsObjects, SeisflowsParameters, SeisflowsPaths

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()
OBJ = SeisflowsObjects()



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

    def save_objects(self):
        OBJ.save(PATH.OUTPUT)

    def save_parameters(self):
        PAR.save(PATH.OUTPUT)

    def save_paths(self):
        PATH.save(PATH.OUTPUT)

