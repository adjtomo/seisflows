#!/usr/bin/env python3
"""
Miyabi is a supercomputer operated by the Joint Center for Advanced High
Performance Computing (JCAHPC), a collaboration between the University of
Tsukuba and the University of Tokyo. Miyabi runs on the PBS Professional job
management system and therefore overloads the `Pbs` System module.
Miyabi-specific parameters and functions are defined here.

Miyabi is split into two independently-scheduled subsystems which do NOT
share a login node or a CPU architecture, so Users must build/run their
solver against the matching subsystem:

- Miyabi-G: 1,120 GPU-accelerated compute nodes, each with a single NVIDIA
  GH200 Grace-Hopper Superchip (1x Arm Neoverse V2 CPU, 72 cores, + 1x H100
  GPU). Login node (Login-G) and compute node CPU architecture is aarch64.
- Miyabi-C: 190 CPU-only compute nodes, each with two Intel Xeon MAX 9480
  CPUs (56 cores each, 112 cores/node). Login node (Login-C) and compute
  node CPU architecture is x86_64.

Information on Miyabi can be found in the Miyabi User's Guide (this class
was written against v1.8, Fujitsu Limited, 2025), and on the JCAHPC website:
https://miyabi.jcahpc.jp/

.. note:: Miyabi Caveat 1

    Miyabi's job management system does not support multibyte characters.
    Job names, project codes, and queue names passed to 'qsub' must be
    plain ASCII (User's Guide Sec. 5.3.1), enforced here in `check()`.

.. note:: Miyabi Caveat 2

    Only "node-occupied use" queues (whole compute nodes) are supported by
    this module. Miyabi-G additionally provides "MIG use" queues which
    allocate fractional (1/4) GH200 GPUs; these are not currently supported
    here since SPECFEM-style multi-node MPI workflows are expected to use
    whole nodes.

.. note:: Miyabi Caveat 3

    Submitting a job to a project/allocation that has run out of "tokens"
    (Miyabi's core-hour equivalent) will cause the job to fail at
    submission or execution time. Token usage can be checked on the login
    node with the `show_token` command; SeisFlows does not manage this.

.. note:: Miyabi Caveat 4

    Like Wisteria, Miyabi's compute nodes do NOT inherit the login node's
    Conda environment, so a bare 'qsub'-submitted job will not have access
    to the Python environment used to install SeisFlows. If `conda_env` is
    set, `submit_workflow`/`run_functions` are routed through the shared
    `runscripts/conda_activate-miyabi` wrapper (handling both entry points),
    which activates Conda (loaded as a module by default on Miyabi) before
    handing off to the real SeisFlows entry point. If `conda_env` is not
    set, jobs are submitted directly and must already have Conda activated
    some other way (e.g., inherited via `-V`, if your site allows it).

    This mechanism, and this module's defaults more generally, are
    currently only set up/tested for Miyabi-G (GPU). Miyabi-C support
    should work by choosing a '-c' `queue`, but has not been exercised.
"""
import os
from seisflows import ROOT_DIR
from seisflows.system.pbs import Pbs


class Miyabi(Pbs):
    """
    System Miyabi
    -------------
    JCAHPC HPC Miyabi (Univ. of Tsukuba / Univ. of Tokyo), PBS Professional
    based system

    Parameters
    ----------
    :type queue: str
    :param queue: Name of the (node-occupied use) queue used for job
        submission. `select` (i.e., number of nodes, calculated internally
        from `nproc`/`node_size`) determines which of the following
        queues is ultimately used when submitting to a routing queue
        ('regular-g'/'regular-c'); Users may also submit directly to a
        specific queue. Available queues are:

        Miyabi-G (GPU/aarch64, 72 cores + 1x GH200 GPU per node):
            - debug-g: 30 min max, [1, 16] nodes
            - short-g: 8 hr max, [1, 8] nodes
            - regular-g: routes to one of the following based on `select`
                - small-g: 48 hr max, [1, 16] nodes
                - medium-g: 48 hr max, [17, 64] nodes
                - large-g: 48 hr max, [65, 128] nodes
                - x-large-g: 24 hr max, [129, 256] nodes

        Miyabi-C (CPU/x86_64, 112 cores per node):
            - debug-c: 30 min max, [1, 4] nodes
            - short-c: 8 hr max, [1, 2] nodes
            - regular-c: routes to one of the following based on `select`
                - small-c: 48 hr max, [1, 16] nodes
                - medium-c: 48 hr max, [17, 32] nodes
                - large-c: 48 hr max, [33, 64] nodes

        (see User's Guide Table 5-5 and Table 5-7 for the full/current
        listing, including token/budget-restricted node limits.)
    :type group: str
    :param group: Project code used to charge tokens for job execution, set
        with the PBS '-W group_list=<group>' option. Required by Miyabi for
        any job that consumes tokens.
    :type submit_to: str
    :param submit_to: (Optional) queue used to submit the main/master job,
        which is a serial Python task that controls the workflow. Likely
        this should be 'debug-g'/'debug-c' for small jobs. If not given,
        defaults to `queue`.
    :type conda_env: str
    :param conda_env: (Optional) name or full path of the Conda environment
        SeisFlows is installed in. If given, job submission/execution is
        routed through `runscripts/conda_activate-miyabi`, which activates
        this environment on the compute node before running SeisFlows (see
        Caveat 4). Required on Miyabi unless Conda is made available to
        compute node jobs some other way.

    Paths
    -----
    ***
    """
    __doc__ = Pbs.__doc__ + __doc__

    def __init__(self, mpiexec="mpirun", queue="debug-g", group=None,
                 submit_to=None, conda_env=None, **kwargs):
        """Miyabi init"""
        super().__init__(**kwargs)

        self.mpiexec = mpiexec
        self.queue = queue
        self.group = group
        self.submit_to = submit_to or self.queue
        self.conda_env = conda_env

        # Node-occupied use queues only, see Caveat 2 above. Node sizes
        # (cores/node) are taken from the hardware specs in User's Guide
        # Sec. 1.2 (Table 1-1, Table 1-2); queue names come from Sec. 5.2.1
        # (Table 5-5, Table 5-7). Queue-specific node count/walltime limits
        # are documented above but not hard-enforced here, since they may be
        # revised by JCAHPC independently of this module.
        self._queues = {
            # Miyabi-G: NVIDIA GH200 (Arm Neoverse V2, 72 cores/node)
            "debug-g": 72, "short-g": 72, "regular-g": 72,
            "small-g": 72, "medium-g": 72, "large-g": 72, "x-large-g": 72,
            # Miyabi-C: 2x Intel Xeon MAX 9480 (56+56 cores/node)
            "debug-c": 112, "short-c": 112, "regular-c": 112,
            "small-c": 112, "medium-c": 112, "large-c": 112,
        }

    def check(self):
        """
        Checks parameters and paths
        """
        super().check()

        assert(self.group is not None), (
            f"Miyabi requires a project code to charge tokens for job "
            f"execution, set with parameter `group` "
            f"(PBS '-W group_list=<group>')")

        # Miyabi Caveat 1: job management system does not support multibyte
        # (i.e., non-ASCII) characters anywhere in the qsub directives
        for val, name in [(self.title, "title"), (self.group, "group"),
                           (self.queue, "queue"),
                           (self.submit_to, "submit_to")]:
            assert(str(val).isascii()), (
                f"Miyabi's job management system does not support "
                f"multibyte characters, but `system.{name}`=='{val}' "
                f"contains non-ASCII characters"
            )

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
        ('submit' or 'run'). Miyabi's compute nodes do not inherit the
        login node's Conda environment (Caveat 4), so if `conda_env` is
        set, the call is routed through the shared
        `conda_activate-miyabi` wrapper script, which activates the given
        Conda environment before executing the real entry point script
        (with all further arguments forwarded to it unchanged). If
        `conda_env` is not set, the entry point script is called directly.

        :type name: str
        :param name: entry point script name, 'submit' or 'run'
        :rtype: str
        :return: command used to invoke the given entry point
        """
        if not self.conda_env:
            return os.path.join(ROOT_DIR, "system", "runscripts", name)

        wrapper = os.path.join(ROOT_DIR, "system", "runscripts",
                               "conda_activate-miyabi")
        return f"{wrapper} {name} {self.conda_env}"
