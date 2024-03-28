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
    nodes, so it requires custom `submit` and `run` script which first load
    the correct modules, and then run the corresponding script

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
    :type gpu: int
    :param gpu: if not None, tells SeisFlows to use the GPU version of SPECFEM, 
        the integer value of `gpu` will set the number of requested GPUs for a 
        simulation on system (i.e., #PJM -L gpu=`gpu`)

    Paths
    -----

    ***
    """
    __doc__ = Fujitsu.__doc__ + __doc__

    # Overwrites submit and run call locations to provide custom scripts that
    # are used to deal with the issue that the Conda environment is not 
    # automatically inherited on a compute node
    submit_workflow = os.path.join(ROOT_DIR, "system", "runscripts",
                                   "custom_submit-wisteria")   
    run_functions = os.path.join(ROOT_DIR, "system", "runscripts", 
                                 "custom_run-wisteria")   

    def __init__(self, group=None, rscgrp=None, gpu=None, **kwargs):
        """Wisteria init"""
        super().__init__(**kwargs)

        self.group = group
        self.rscgrp = rscgrp
        self.gpu = gpu

        # Wisteria resource groups and their cores per node
        self._rscgrps = {
                # Node-occupied resource allocation (Odyssey)
                "debug-o": 48, "short-o": 48, "regular-o": 48, "priority-o": 48,
                # Node-occupied resource allocation (Aquarius)
                "debug-a": 36, "short-a": 36, "regular-a": 36, 
                # GPU-exclusive resource allocation
                "share-debug": 1, "share-short": 2, "share": 5
                }

    @property
    def run_call_header(self):
        """Override run call header to allow for GPU version requirements"""
        if self.gpu:
            _call = super().run_call_header.replace(
                    f"-L node={self.nodes}",
                    f"-L gpu={self.gpu}"
                    )
            return _call
        else:
            return super().run_call_header


