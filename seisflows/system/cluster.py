#!/usr/bin/env python3
"""
The Cluster class provides the core utilities interaction with HPC systems
which must be overloaded by subclasses for specific workload managers, or
specific clusters.

.. warning::
    The Cluster class is an abstract base class for the Systems module which and
    MUST be overwritten by system-specific child classes, it cannot be used to
    run jobs by itself.
"""
import subprocess
from seisflows.system.workstation import Workstation


class Cluster(Workstation):
    """
    [system.cluster] Generic or common HPC/cluster interfacing commands

    :type walltime: int
    :param walltime: maximum job time in minutes for the master SeisFlows
        job submitted to cluster
    :type tasktime: int
    :param tasktime: maximum job time in minutes for each job spawned by
        the SeisFlows master job during a workflow. These include, e.g.,
        running the forward solver
    :type environs: str
    :param environs: Optional environment variables to be provided in the
        following format VAR1=var1,VAR2=var2... Will be set using
        os.environs
    """
    __doc__ = Workstation.__doc__ + __doc__

    def __init__(self, walltime=10, tasktime=1, environs="", **kwargs):
        """Instantiate the Cluster System class"""
        super().__init__(**kwargs)

        self.walltime = walltime
        self.tasktime = tasktime
        self.environs = environs

    def submit(self, workflow, submit_call=None):
        """
        Main insertion point of SeisFlows onto the compute system.

        .. rubric::
            $ seisflows submit

        .. note::
            The expected behavior of the submit() function is to:
            1) run system setup, creating directory structure,
            2) execute workflow by submitting workflow.main()

        :type submit_call: str
        :param submit_call: the command line workload manager call to be run by
            subprocess. These need to be passed in by specific workload manager
            subclasses.
        """
        self.setup()
        workflow.checkpoint()
        # check==True: subprocess will wait for workflow.main() to finish
        subprocess.run(submit_call, shell=True, check=True)

    def run(self, classname, method, single=False, **kwargs):
        """
        Runs a task multiple times in parallel

        .. note::
            The expected behavior of the run() function is to: submit N jobs to
            the system in parallel. For example, in a simulation step, run()
            submits N jobs to the compute system where N is the number of
            events requiring an adjoint simulation.

        :type classname: str
        :param classname: the class to run
        :type method: str
        :param method: the method from the given `classname` to run
        :type single: bool
        :param single: run a single-process, non-parallel task, such as
            smoothing the gradient, which only needs to be run by once.
            This will change how the job array and the number of tasks is
            defined, such that the job is submitted as a single-core job to
            the system.
        """
        raise NotImplementedError('Must be implemented by subclass.')

    def taskid(self):
        """
        Provides a unique identifier for each running task. This is
        compute system specific.

        :rtype: int
        :return: this function is expected to return a unique numerical
            identifier.
        """
        raise NotImplementedError('Must be implemented by subclass.')

