#!/usr/bin/env python
"""
This is the base class seisflows.system.Base
This class provides the core utilities interaction with HPC systems which must
be overloaded by subclasses
"""
import os

from seisflows.config import save, saveobj
from seisflows.tools import unix


class Base(object):
    """
    Abstract base class
    """
    def check(self):
        """
        Checks parameters and paths
        """
        raise NotImplementedError('Must be implemented by subclass.')

    def submit(self):
        """
        Submits workflow
        """
        raise NotImplementedError('Must be implemented by subclass.')

    def run(self, classname, method, *args, **kwargs):
        """
        Runs task multiple times
        """
        raise NotImplementedError('Must be implemented by subclass.')

    def run_single(self, classname, method, *args, **kwargs):
        """
        Runs task a single time
        """
        raise NotImplementedError('Must be implemented by subclass.')

    def taskid(self):
        """
        Provides a unique identifier for each running task
        """
        raise NotImplementedError('Must be implemented by subclass.')

    def checkpoint(self, path, classname, method, args, kwargs):
        """
        Writes information to disk so tasks can be executed remotely

        :type path: str
        :param path: path the kwargs
        :type classname: str
        :param classname: name of the class to save
        :type method: str
        :param method: the specific function to be checkpointed
        :type kwargs: dict
        :param kwargs: dictionary to pass to object saving
        """
        argspath = os.path.join(path, "kwargs")
        argsfile = os.path.join(argspath, f"{classname}_{method}.p")
        unix.mkdir(argspath)
        saveobj(argsfile, kwargs)
        save()

