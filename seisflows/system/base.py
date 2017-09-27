

from os.path import join

from seisflows.config import save, saveobj
from seisflows.tools import unix


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



    def run(self, classname, method, *args, **kwargs):
        """ Runs task multiple times
        """
        raise NotImplementedError('Must be implemented by subclass.')



    def run_single(self, classname, method, *args, **kwargs):
        """ Runs task a single time
        """
        raise NotImplementedError('Must be implemented by subclass.')



    def taskid(self):
        """ Provides a unique identifier for each running task
        """
        raise NotImplementedError('Must be implemented by subclass.')



    def checkpoint(self, path, classname, method, args, kwargs):
        """ Writes information to disk so tasks can be executed remotely
        """
        argspath = join(path, 'kwargs')
        argsfile = join(argspath, classname+'_'+method+'.p')
        unix.mkdir(argspath)
        saveobj(argsfile, kwargs)
        save()

