SPECFEM3D Example
=================

Example 4 runs an inversion between a homogeneous halfspace and the same
halfspace with its velocities changed by a few percent, using
SPECFEM3D_Cartesian. By default it is small enough to run on a single compute
node (28 cores).

.. code:: bash

    seisflows examples run 4 --specfem3d_repo /path/to/specfem3d

Nothing is downloaded or compiled. ``--specfem3d_repo`` must point to a
SPECFEM3D repository that you have compiled already, and which contains the
``EXAMPLES/applications/homogeneous_halfspace`` example. SeisFlows needs the
executables ``xmeshfem3D``, ``xgenerate_databases``, ``xspecfem3D``,
``xcombine_sem``, ``xsmooth_sem_pde`` and ``xcombine_vol_data_vtk`` in ``bin/``,
and says which are missing if any are.

What the example does
---------------------

1. Makes ``specfem3d_workdir/`` from the SPECFEM3D example, with the mesh
   and the number of cores set for your ``--nproc``.
2. Runs the SPECFEM3D mesher and database generator to make the starting
   model (``OUTPUT_FILES_INIT``), the example's homogeneous halfspace
   (vp 2800 m/s, vs 1500 m/s, rho 2300 kg/m3).
3. Does the same with the velocities (not the density) scaled by
   ``--perturb`` percent to make the target model (``OUTPUT_FILES_TRUE``).
4. Writes the SeisFlows parameter file, and with ``run`` rather than
   ``setup``, submits the workflow.

The domain is 134 x 134 x 60 km with absorbing boundaries on five sides. The
events are 100 s long. All three components (``XYZ``) are used, with a
traveltime misfit.

.. warning::

    Step 2 and 3 run SPECFEM3D with MPI on ``--nproc`` cores. On a cluster, do
    this from an interactive session on a compute node
    (e.g., ``salloc -N 1`` or ``qsub -I``), not on the login node.

Arguments
---------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Argument
     - Meaning
   * - ``--specfem3d_repo``
     - Compiled SPECFEM3D repository. Defaults to ``specfem3d/`` in the
       current directory.
   * - ``--nproc``
     - Cores per simulation, default 28. Sets how the mesh is divided (see
       below).
   * - ``--ntask``
     - Number of events, default 1. The first is at the center of the domain,
       the rest are evenly spaced on a ring around it.
   * - ``--nsta``
     - Number of stations, default 25. They are a grid of equidistant points on
       the surface that fills the domain except for a margin along the edges
       (25 is 5 x 5). If ``nsta`` is not a square, the last row is incomplete.
   * - ``--niter``
     - Number of iterations, default 1.
   * - ``--perturb``
     - Percent change of the target velocities, default 5. Use a negative
       number for slower, e.g., ``--perturb -5``.
   * - ``--mpiexec``
     - MPI executable used to make the models, default ``mpirun``. On Chinook
       use ``mpiexec``.
   * - ``--system``
     - SeisFlows system module, default ``workstation``. See below.
   * - ``--par``
     - SeisFlows parameters to set, as ``KEY=VALUE``. See below.
   * - ``--par_yaml``
     - YAML file of SeisFlows parameters to set. See below.

Fitting a node
--------------

The mesh is divided between the cores as evenly as possible along its two
horizontal directions, and the number of elements follows, so that they stay
about 3.7 km wide:

.. list-table::
   :header-rows: 1

   * - ``--nproc``
     - Slices
     - Elements
   * - 4
     - 2 x 2
     - 36 x 36 x 16
   * - 24
     - 6 x 4
     - 36 x 36 x 16
   * - 28
     - 7 x 4
     - 35 x 36 x 16
   * - 36
     - 6 x 6
     - 36 x 36 x 16
   * - 40
     - 8 x 5
     - 40 x 35 x 16

Choose ``--nproc`` to be the size of a node of the partition or queue that you
will use (e.g., 28 on Chinook), or smaller.

Running on a cluster
--------------------

By default the example runs on the machine that it is run on (the
``workstation`` system). To run on a cluster, choose its system module with
``--system`` and set the parameters that it requires with ``--par``. The
parameters depend on the system: ``seisflows par`` lists them once the
parameter file exists, and the example stops with a suggestion if one of yours
does not exist.

.. code:: bash

    seisflows examples run 4 --specfem3d_repo ~/specfem3d --mpiexec mpiexec \
        --system chinook --par partition=t1small submit_to=debug \
        walltime=30 tasktime=40 blocking=False

Settings that you use again, such as everything that your cluster needs, are
better kept in a YAML file, which is also where you can put any other
parameter of the workflow. ``--par`` is applied after it, so it can be used to
change one value.

.. code:: yaml

    # chinook.yaml
    partition: t1small
    submit_to: debug
    walltime: 30
    tasktime: 40
    blocking: False

.. code:: bash

    seisflows examples run 4 --specfem3d_repo ~/specfem3d --mpiexec mpiexec \
        --system chinook --par_yaml chinook.yaml

``blocking: False`` lets the main job finish while the simulations are
queued and run, which is how a small queue can be used for it (see
`running on Chinook <running_on_chinook.html>`__ and `running on Miyabi
<running_on_miyabi.html>`__ for the other parameters of those systems).

The example's own parameters can be changed in the same way, e.g.,
``--par step_count_max=3``. To only set up and look at the parameter file
before running, use ``seisflows examples setup 4`` instead of ``run``, edit
``parameters.yaml`` (``seisflows par`` or ``seisflows swap system chinook``
also work) and run ``seisflows submit``.
