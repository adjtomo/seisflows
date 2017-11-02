
1. Download SeisFlows
---------------------

To run SeisFlows you'll need a Unix system with Python 2.7, Numpy, Scipy, Obspy and standard Unix utilities.  After these prerequisites are in place, from the command line type::
 
        mkdir ~/packages
        cd ~/packages
        git clone https://github.com/rmodrak/seisflows.git

If you prefer a location other than ``~/packages`` , modify the commands above and below accordingly.


2. Set environment variables
----------------------------

Add the following lines to ``~/.bashrc`` (modify accordingly, if you are using a shell other than bash)::

        export PATH=$PATH:~/packages/seisflows/scripts
        export PYTHONPATH=~/packages/seisflows
 

Don't forget to update any open shells::

        source ~/.bashrc
 

 

3. Run "system" test
---------------------

 
Run the following test to make sure everything is working::

        ~/packages/seisflows/tests/run_test_system


If a ''hello'' message is displayed, the test was successful.

 

 

4. Run nonlinear optimization test
----------------------------------


Run the following test to make sure everything is working::

        ~/packages/seisflows/tests/run_test_optimize


If the optimization problems are solved in 60 iterations or fewer, the test was successful.

 

 

5. Configure and compile SPECFEM2D
----------------------------------

First, download SPECFEM2D from GitHub::

        cd ~/packages
        git clone --recursive --branch devel https://github.com/geodynamics/specfem2d.git specfem2d-d745c542
        cd specfem2d-d745c542
        git checkout d745c542

For now, it is important to work with the exact version specified above (``d745c542``). This is necessary because, unlike SPECFEM3D and 3D_GLOBE, SPECFEM2D development is sometimes a bit haphazard, with frequent interface changes.


Next, configure and compile SPECFEM2D using ifort (preferred) or some other fortran compiler::

        cd ~/packages/specfem2d-d745c542
        ./configure FC=ifort
        make all

(Since `make` by itself does not compile all the required utilities, be sure to remember to type `make all`.)  For troubleshooting any compilation issues, please view SPECFEM2D's manual and GitHub issues page.
 


6. Set up checkerboard test
---------------------------

Download the starting model and other input files required for the waveform inversion checkerboard test.  For simplicity, let's assume the checkerboard working directory will be placed in ``~/tests`` (if you prefer a different location, then modify the following commands accordingly)::
 
        mkdir ~/tests/
        cd ~/tests/
        wget --recursive --no-parent --no-host-directories --cut-dirs=2 --reject "index.html*" http://tigress-web.princeton.edu/~rmodrak/2dAcoustic/


A directory ``~/tests/checkers`` is now being created.  Among other files, ``parameters.py`` and ``paths.py`` are being downloaded.

After the download completes, make sure that all paths specified in ``paths.py``  are correct.  For example, if you compiled SPECFEM2D somewhere other than ``~/packages/specfem2d-d745c542``, you will need to modify the ``SPECFEM2D_BIN`` entry accordingly. 

Next, take a minute to view the ``parameters.py`` file and note the close similarity between the first set of parameters and the `directory structure <https://github.com/PrincetonUniversity/seisflows/tree/master/seisflows>`_ of the SeisFlows repository.

Note: File hosting services are provided by my alma mater.  The download server may become temporarily unavailable due to system maintenance or permanently unavailable due to expiration of my account.

 
7. Run checkerboard test in serial
----------------------------------

To run the checkerboard test type::

        sfclean ; sfrun

within ``~/tests/checkers``.

For now, the inversion will run only a single event on only a single processor.  Once we verify that everything is working correctly in this case, we can move on to multiple events and multiple processors by modifying ``parameters.py``.



8. Run checkerboard test in parallel
-----------------------------------------
On a laptop or desktop with multiple cores, the work of an inversion can be carried out in parallel.  To run the checkerboard example in parallel over events (that is, with multiple event simulations running at the same time on different cores), make the following changes to ``parameters.py``:

- to invert all available events instead of just one event, change ``NTASK`` from ``1`` to ``25``
- change ``SYSTEM`` from ``serial`` to ``multithreaded``
- add a parameter ``NPROCMAX`` and set it to the number of cores available on your machine.

Besides running in parallel over events, the work of an individual event simulation can be parallelized over model regions. See the SPECFEM3D user manual for more information. Both parallelization over events and over model regions can be used at the same time under SeisFlows.  The current example, however, illustrates only event parallelism.

Besides ``serial`` and ``multithreaded`` settings for running SeisFlows on laptops and desktops, there are also PBS, SLURM, and LSF options for running on clusters. See `here <http://seisflows.readthedocs.org/en/latest/usage/usage.html#system-configuration>`_ for more information.


9. Visualize inversion results
------------------------------

Visualization requires software such as Pylab, Matlab, or Paraview.

With any such software, one approach for plotting SPECFEM2D models or kernels is to interpolate from the unstructured numerical mesh on which the model parameters are defined to a uniform rectangular grid.  The Pylab script `plot2d <http://tigress-web.princeton.edu/~rmodrak/visualize/plot2d>`_ illustrates this approach.


Another method is to compute a Delaunay triangulation and plot the model or kernel over the unstructured mesh itself.  A Pylab script `plot2d_delaunay <http://tigress-web.princeton.edu/~rmodrak/visualize/plot2d_delaunay>`_ is available for illustration.

To plot results from the checkerboard example using ``plot2d``, run the following command from the working directory::

          plot2d output/model_init/proc000000_x.bin \
                 output/model_init/proc000000_z.bin \
                 output/model_0001/proc000000_vs.bin

(The command line syntax is the same for the other script.)  For either script to work, Pylab must be installed and the Pylab backend properly configured. If you prefer visualization software other than Pylab, feel free to use the above scripts for reference in writing your plotting own tools. 


10. Creating your own examples
------------------------------
It may be clear by now that with SeisFlows, wave simulations must be performed using an external software package such as SPECFEM2D or SPECFEM3D.  The ability to interface with external solvers ensures flexibility, and the choice of SPECFEM as a default option gives access to cutting-edge meshing and hardware accelaration capabilities.  However, the use of external package also creates additional work for the user because, to carry set up one's own inversion, one must become familiar not only with the SeisFlows package, but also with a separate solver package.  

To move beyond the above checkerboard test case, familiarity with how to set up simulations with SPECFEM--in paricular with how to create models in SPECFEM's idionsyncratic binary format--is essential.  `Issue #83 <https://github.com/rmodrak/seisflows/issues/83>`_ may be helpful in this regard.  Trying the two other `examples available for download <https://github.com/rmodrak/seisflows/blob/master/docs/index.rst#examples-available-for-download>`_ may also be useful.
