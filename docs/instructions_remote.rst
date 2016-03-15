
1. Download SeisFlows
---------------------

To run SeisFlows you'll need a Unix system with Python 2.7, Numpy, Scipy, and standard Unix utilities.  After these prerequisites are in place, from the command line type::
 
        mkdir /home/packages
        cd /home/packages
        git clone https://github.com/PrincetonUniversity/seisflows.git

If you prefer a location other than ``/home/packages`` , modify the commands above and below accordingly.


2. Set environment variables
----------------------------

Add the following lines to ``/home/.bash_profile`` (modify accordingly, if you are using a shell other than bash)::

        export PATH=$PATH:/home/packages/seisflows/scripts
        export PYTHONPATH=$PYTHONPATH:/home/packages/seisflows
 

Don't forget to update any open shells::

        source /home/.bash_profile
 

 

3. Run "system" test
---------------------

 
Run the following test to make sure everything is working::

        cd /home/packages/seisflows/tests/integration/system
        ./clean.py; ./run.py


If a ''hello'' message is displayed, the test was successful.

 

 

4. Run nonlinear optimization test
----------------------------------


Run the following test to make sure everything is working::

        cd /home/packages/seisflows/tests/integration/optimize
        ./clean.py; ./run.py


If the optimization problem is solved in 50 iterations or fewer, the test was successful.

 

 

5. Configure and compile SPECFEM2D
----------------------------------

First, download SPECFEM2D from github::

        cd /home/packages
        git clone --recursive --branch devel https://github.com/geodynamics/specfem2d.git specfem2d-d745c542
        cd specfem2d-d745c542
        git checkout d745c542

For now, it is important to work with the exact version specified above (``d745c542``). This is necessary because, unlike SPECFEM3D and 3D_GLOBE, SPECFEM2D development is sometimes a bit haphazard, with frequent interface changes.


Next, configure and compile SPECFEM2D using ifort (preferred) or gfortran::

        cd /home/packages/specfem2d-d745c542
        ./configure FC=ifort
        make all

For troubleshooting any compilation issues, please view the SPECFEM2D manual and github issues page.
 


6. Set up checkerboard test
---------------------------

Download the starting model and other input files required for the waveform inversion checkerboard test.  For simplicity, let's assume the checkerboard working directory will be placed in ``/home/tests`` (if you prefer a different location, then modify the following commands accordingly)::
 
        mkdir /home/tests/
        cd /home/tests/
        wget --recursive --no-parent --no-host-directories --cut-dirs=2 --reject "index.html*" http://tigress-web.princeton.edu/~rmodrak/2dAcoustic/


A directory ``/home/tests/checkers`` is now being created.  Among other files, ``parameters.py`` and ``paths.py`` are being downloaded.

After the download completes, make sure that all paths specified in ``paths.py``  are correct.  For example, if you compiled SPECFEM2D somewhere other than ``/home/packages/specfem2d-d745c542``, you will need to modify the ``SPECFEM2D_BIN`` entry accordingly.

 
7. Run checkerboard test in serial
----------------------------------

To run the checkerboard test type::

        sfclean ; sfrun

within ``/home/tests/checkers``.

For now, the inversion will run only a single event on only a single processor.  Once we verify that everything is working correctly in this case, we can move on to multiple events and multiple processors by modifying ``parameters.py``.



8. Run checkerboard test in parallel
-----------------------------------------
On a laptop or desktop with multiple cores, the work of an inversion can be carried out in parallel.  To run the checkerboard example in parallel over events (that is, with multiple event simulations running at the same time on different cores), make the following changes to ``parameters.py``:

- to invert all available events instead of just one event, change ``NTASK`` from ``1`` to ``25``
- change ``SYSTEM`` from ``serial`` to ``multithreaded``
- add a parameter ``NPROCMAX`` and set it to the number of cores available on your machine.

Besides running in parallel over events, the work of an individual event simulation can be parallelized over model regions. See the SPECFEM3D user manual for more information. Both parallelization over events and over model regions can be used at the same time under SeisFlows.  The current example, however, only illustrates parallelization over events.

Besides ``serial`` and ``multithreaded`` settings for running SeisFlows on laptops and desktops, there are also PBS, SLURM, and LSF options for running on clusters. See `here <http://seisflows.readthedocs.org/en/latest/manual/manual.html#system-configuration>`_ for more information.


9. Visualize inversion results
------------------------------

Visualization requires software such as Pylab, Matlab, or Paraview.

With any such software, one approach for plotting SPECFEM2D models or kernels is to interpolate from the unstructured numerical mesh on which the model parameters are defined to a uniform rectangular grid.  The Pylab script `plot2d <http://tigress-web.princeton.edu/~rmodrak/visualize/plot2d>`_ illustrates this approach.


Another method is to compute a Delaunay triangulation and plot the model or kernel over the unstructured mesh itself.  A Pylab script `plot2d_delaunay <http://tigress-web.princeton.edu/~rmodrak/visualize/plot2d_delaunay>`_ is available for illustration.

To plot results from the checkerboard example using ``plot2d``, run the following command from the working directory::

          plot2d output/model_init/proc000000_x.bin \
                 output/model_init/proc000000_z.bin \
                 output/model_0001/proc000000_vs.bin

(The command line syntax is the same for the other script, ``plot2d_delaunay``.)  For either script to work, Pylab must be installed and the Pylab backend properly configured. If you prefer visualization software other than Pylab, feel free to use the above scripts for reference in writing your plotting tools. 
