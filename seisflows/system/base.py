
from seisflows.config import save

class base(object):
    """ Abstract base class
    """

    def check(self):
        """ Checks parameters and paths
        """
        raise NotImplementedError('Must be implemented by subclass.')



    def submit(self):
        """ Submits workflow
        """
        raise NotImplementedError('Must be implemented by subclass.')



    def run(self, classname, method, **kwargs):
        """ Executes the following task:
              classname.method(*args, **kwargs)
        """
        raise NotImplementedError('Must be implemented by subclass.')



    def taskid(self):
        """ Provides a unique identifier for each running task
        """
        raise NotImplementedError('Must be implemented by subclass.')



    def checkpoint(self):
        """ Writes information to disk so workflow can be resumed in case of
          interruption        
        """
        save()

