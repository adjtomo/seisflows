
from seisflows.config import save

class base(object):
    """ Workflow abstract base class
    """


    def check(self):
        """ Checks parameters and paths
        """
        raise NotImplementedError('Must be implemented by subclass.')



    def main(self):
        """ Main routine

          Execution of a workflow is equivalent to stepping through 
          workflow.main
        """
        raise NotImplementedError('Must be implemented by subclass.')


    def checkpoint(self):
        """ Writes information to disk so workflow can be resumed following a
          break
        """
        save()


