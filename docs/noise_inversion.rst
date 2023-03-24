Ambient Noise Adjoint Tomography
================================

SeisFlows contains workflows that allows Users to use numerical solvers
to invert for Empirical Greens Functions (EGFs) derived from ambient noise
cross-correlations. We call this ambient noise adjoint tomography (ANAT).

The methdology of this workflow follows the work of Wang et al., that is, it
assumes noise sources are uniformly distributed and treats a given station `R`
as a virtual source. At each virtual source, a point force is applied in
orthogonal coordinates and the waveforms recorded at corresponding components
for all other stations are treated as Synthetic Greens Functions (SGFs).

EGFs and SGFs are then compared to generate adjoint sources which can be used
for kernel generation and iterative inversions.


Setup
-----

The following section describes how a User must set up their working directory
and parameter file to run an ANAT inversion. To set up your parameter file,
use the `noise_inversion` workflow module.

.. code:: bash

    seisflows par workflow noise_inversion
    seisflows configure

Data Path
`````````

Parameter `path_data` points to where your EGFs are stored,  `path_data` must
be set to a directory containing your EGFs, separated by station name and EGF
type.

.. note::

    The ANAT workflow is currently only set up to work with real data (EGF) and
    for a few different EGF combinations: ZZ, RR and TT. See !!! REFERENCE !!!
    for interpretations of these EGF.

An example directory structure should follow your `STATIONS` file. If your
STATIONS file looks like this:

.. code::

  S000000 AA 2.43610e+05 2.78904e+05 0.0 0.0
  S000001 AA 3.38981e+05 1.77849e+05 0.0 0.0
  S000002 AA 1.64438e+05 2.94733e+05 0.0 0.0

then you must have a corresponding `path_data` which has subdirectories
formatted NN.SSS/ (N=network, S=station), and further subdirectories separating
EGF kernels. An example directory looks like:

.. code:: bash

      EGF/
      ├── AA.S000000/
      |   ├── ZZ/
      |   |   ├── AA.S000001.BXY.semd
      |   |   ├── AA.S000002.BXY.semd
      │   |   └── ...
      |   ├── RR/
      │   └── TT/
      ├── AA.S000001/
      |   ├── ZZ/
      |   ├── RR/
      │   └── TT/
      └── AA.S000002/
          └── ...


Sources
```````

Similar to a normal workflow, SeisFlows requires all sources be predefined in
your working directories `DATA/` directory. This requires that all stations
in your `STATIONS` file have a corresponding source (if there is a corresponding
EGF).

For the example `STATIONS` file above, SeisFlows expects corresponding source
files formatted `<SOURCE_PREFIX>_NNSSS`. For SPECFEM3D/3D_GLOBE, sources must
be FORCESOLUTIONS. For SPECFEM2D, sources must be SOURCE.

.. code:: bash

     specfem_workdir/
     └── DATA/
         ├── Par_file
         ├── FORCESOLUTION_AAS000000
         ├── FORCESOLUTION_AAS000001
         ├── FORCESOLUTION_AAS000002
         └── ...

Each force solution should have the same coordinates as it's corresponding
station so that the force is input at the station location.

SeisFlows will take care of setting the correct force vector directions
internally, so Users only have to set the locations.

You can take advantage of a built-in Python function for converting a `STATIONS`
file into corresponding source files. The function !!! AUTOLINK HERE!!!  simply
scans through a `STATIONS` file and replaces coordinates in a corresponding
source file with station locations.

.. code:: python

    from seisflows.tools.specfem import covert_stations_to_sources

    convert_stations_to_sources(stations_file="./DATA/STATIONS",
                                source_file="./DATA/FORCESOLUTION",
                                source_type="FORCESOLUTION",
                                output_dir="./")

Kernel Selection
````````````````

The parameter `kernels` is unique to the noise workflow and determines which
SGF SeisFlows should create, this is directly tied to what EGF data you have
available.

Current available options are ZZ (vertical-vertical),
TT (transverse-tranvserse), and RR (radial-radial). The kernel parameter must
be input as a comma-separated list.

.. code:: bash

    seisflows par kernels ZZ,TT

In the parameter file this looks like

.. code:: yaml

    kernels: ZZ,TT



