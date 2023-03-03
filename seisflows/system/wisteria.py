#!/usr/bin/env python3
"""
Wisteria is the University of Tokyo Fujitsu brand high performance computer.
Wisteria runs on the Fujitsu/PJM job scheduler.

.. notes::

    - Wisteria has two node gruops, Odyssey (compute nodes) and Aquarius 
      (data/learning nodes w/ GPU)
    - Odyssey has 7680 nodes with 48 cores/node
    - Aquarius has 45 nodes with 36 cores/node
"""
import os
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

    Paths
    -----

    ***
    """
    __doc__ = Fujitsu.__doc__ + __doc__


    def __init__(self, group=None, rscgrp=None, **kwargs):
        """Wisteria init"""
        super().__init__(**kwargs)

        self.group = group
        self.rscgrp = rscgrp

        # Wisteria resource groups and their cores per node
        self._rscgrps = {
                "debug-o": 48, "short-o": 48, "regular-o": 48, "priority-o": 48,
                "debug-a": 48, "short-a": 48, "regular-a": 48
                }

