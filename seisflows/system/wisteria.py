#!/usr/bin/env python3
"""
Wisteria is the University of Tokyo Fujitsu brand high performance computer.
Wisteria runs on the Fujitsu/PJM job scheduler.

.. note::

    - Wisteria has two node groups, Odyssey (compute nodes) and Aquarius 
      (data/learning nodes w/ GPU)
    - Odyssey has 7680 nodes with 48 cores/node
    - Aquarius has 45 nodes with 36 cores/node
    - Aquarius also contains 8x Nvidia A100

.. note:: Wisteria Caveat 1     
                                                 
    On Wisteria you cannot submit batch jobs from compute nodes and you cannot
    SSH from compute nodes (Manual 5.13), so the master job must be 
    run from the login node or the pre-post node (Manual 5.2.3)

.. note:: Wisteria Caveat 2

    On Wisteria, the login node Conda environment is not inherited by compute
    nodes, and `module`/`conda` are not usable until the relevant modules are
    loaded. If `conda_env` is set, `submit_workflow`/`run_functions` are
    routed through the shared `runscripts/conda_activate-wisteria` wrapper
    (handling both entry points), which loads `conda_modules` (defaulted
    from `gpu`, see Parameters below) and activates Conda before running the
    corresponding SeisFlows entry point.

.. note:: Wisteria Caveat 3

    On Wisteria, command line arguments for the `submit` and `run` script,
    normally input like '--key value' interfere with the batch submission cmd
    `pjsub`. So instead we use the `pjsub` '-x' flag which allows us to set
    environment variables. We use these in place of command line arguments
"""
import os
from seisflows import ROOT_DIR
from seisflows.system.fujitsu import Fujitsu


class Wisteria(Fujitsu):
    """
    System Wisteria
    ---------------
    University of Tokyo HPC Wisteria, running Fujitsu job scheduler

    Parameters
    ----------
    :type group: str
    :param group: User's group for allocating and charging resources. In the 
        pjsub script this is the '-g' option.
    :type rscgrp: str
    :param rscgrp: the resource group (i.e., partition) to submit jobs to. In 
        the pjsub script this is the '-L rscgrp' option.
        Available `rscgrp`s for Wisteria are:

        - debug-o: Odyssey debug, 30 min max, [1, 144] nodes available
        - short-o: Odyssey short, 8 hr. max, [1, 72] nodes available
        - regular-o: Odyssey regular, 24-48 hr. max, [1, 2304] nodes available
        - priority-o: Odyssey priority, 48 hr. max, [1, 288] nodes available

        - debug-a: Aquarius debug, 30 min max, [1, 1] nodes available
        - short-a: Aquarius short, 2 hr. max, [1, 2] nodes available
        - regular-a: Aquarius regular, 24-48 hr. max, [1, 8] nodes available

        - share-debug: Aquarius GPU debug, 30 min max, 1, 2, 4 GPU available
        - share-short: Aquarius GPU short queue, 2 hr. max, 1, 2, 4 GPU avail.
        - share: Aquarius GPU-exclusive, available 1, 2 and 4 GPU, select 
          using `gpu`
    :type gpu: int
    :param gpu: if not None, tells SeisFlows to use the GPU version of SPECFEM,
        the integer value of `gpu` will set the number of requested GPUs for a
        simulation on system (i.e., #PJM -L gpu=`gpu`). Required if using
        GPU-exclusive rscgrps's
    :type conda_env: str
    :param conda_env: name or full path of the Conda environment SeisFlows
        is installed in. If given, job submission/execution is routed
        through `runscripts/conda_activate-wisteria`, which loads
        `conda_modules` and activates this environment on the compute node
        before running SeisFlows (see Caveat 2). Required on Wisteria unless
        Conda is made available to compute node jobs some other way.
    :type conda_modules: str or list
    :param conda_modules: comma-separated string, or list, of `module load`
        targets required before Conda/MPI are usable on the compute node
        (in addition to Conda's own module, which is always loaded). If not
        given, defaults to a GPU- or CPU-appropriate module set based on
        `gpu`: 'cuda/12.2,gcc,ompi' if `gpu` else 'intel,impi'.

    Paths
    -----

    ***
    """
    __doc__ = Fujitsu.__doc__ + __doc__

    def __init__(self, group=None, rscgrp=None, gpu=None, submit_to=None,
                 conda_env=None, conda_modules=None, **kwargs):
        """Wisteria init"""
        super().__init__(**kwargs)

        self.group = group
        self.rscgrp = rscgrp
        self.gpu = gpu
        self.submit_to = submit_to or self.rscgrp
        self.conda_env = conda_env

        if conda_modules is None:
            # Preserves the module sets used by the original (now retired)
            # 'custom_run-wisteria'/'custom_submit-wisteria' scripts, which
            # branched on a 'GPU_MODE' variable derived from `gpu`
            conda_modules = "cuda/12.2,gcc,ompi" if self.gpu else "intel,impi"
        elif not isinstance(conda_modules, str):
            conda_modules = ",".join(conda_modules)
        self.conda_modules = conda_modules

        # Wisteria resource groups and their cores per node
        self._rscgrps = {
                # Node-occupied resource allocation (Odyssey)
                "debug-o": 48, "short-o": 48, "regular-o": 48, "priority-o": 48,
                # Node-occupied resource allocation (Aquarius)
                "debug-a": 36, "short-a": 36, "regular-a": 36,
                # GPU-exclusive resource allocation. Share will give you access
                # to N GPUs based on `gpu`
                "share-debug": 1, "share-short": 2, "share": 1,
                }

    def check(self):
        """
        Checks parameters and paths
        """
        super().check()

        assert(self.conda_env is not None), (
            f"Wisteria requires `conda_env` to be set (name or full path of "
            f"the Conda environment SeisFlows is installed in), since "
            f"compute node jobs do not otherwise have access to Conda "
            f"(see Caveat 2)")

    @property
    def submit_workflow(self):
        """See `_entry_point`. Overwrites `Cluster.submit_workflow`"""
        return self._entry_point("submit")

    @property
    def run_functions(self):
        """See `_entry_point`. Overwrites `Cluster.run_functions`"""
        return self._entry_point("run")

    def _entry_point(self, name):
        """
        Returns the command used to invoke a SeisFlows entry point script
        ('submit' or 'run'), routed through the shared
        `conda_activate-wisteria` wrapper script, which loads
        `conda_modules`, activates `conda_env`, and then executes the real
        entry point script (see Caveat 2 and Caveat 3 -- unlike e.g. Miyabi/
        PBS, arguments are NOT forwarded here since PJM does not tolerate
        them; the wrapper derives them itself for each entry point).

        :type name: str
        :param name: entry point script name, 'submit' or 'run'
        :rtype: str
        :return: command used to invoke the given entry point
        """
        wrapper = os.path.join(ROOT_DIR, "system", "runscripts",
                               "conda_activate-wisteria")
        return f"{wrapper} {name} {self.conda_env} {self.conda_modules}"

