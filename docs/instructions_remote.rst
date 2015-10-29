
1. Download SeisFlows
---------------------

To run seisflows you'll need a Unix system with Python and standard Unix utilities.  From the command line type::
 
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
        ./clean; ./run


If a ''hello'' message is displayed, the test was successful.

 

 

4. Run "optimization" test
--------------------------


Run the following test to make sure everything is working::

        cd /home/packages/seisflows/tests/integration/optimize
        ./clean; ./run


If the optimization problem is solved in 50 iterations or fewer, the test was successful.

 

 

5. Configure and compile SPECFEM2D
----------------------------------

First, download SPECFEM2D from github (to avoid possible version differences, let's use version ``d745c542``)::

        cd /home/packages
        git clone --recursive --branch devel https://github.com/geodynamics/specfem2d.git specfem2d-d745c542
        cd specfem2d-d745c542
        git checkout d745c542


Next, configure and compile SPECFEM2D using ifort (preferred) or gfortran::

        cd /home/packages/specfem2d-d745c542
        ./configure FC=ifort
        make all

For troubleshooting any compilation issues, please view the SPECFEM2D manual and github issues page.
 


6. Set up checkerboard test
-------------------------------

Download the starting model and other input files required for the waveform inversion checkerboard test.  For simplicity, let's assume the checkerboard working directory be placed in ``/home/tests`` (if you prefer a different location, then modify the following commands accordingly)::
 
        mkdir /home/tests/
        cd /home/tests/
        wget --recursive --no-parent --no-host-directories --cut-dirs=2 --reject "index.html*" http://tigress-web.princeton.edu/~rmodrak/2dAcoustic/


A directory ``home/tests/checkers`` is now being created.  Among other files, ``parameters.py`` and ``paths.py`` are being downloaded to this directory.

After the download completes, make sure that all paths specified in ``paths.py``  are correct.  For example, if you compiled SPECFEM2D somewhere other than ``/home/packages/specfem2d-d745c542``, you will need to modify the ``SPECFEM2D_BIN`` entry accordingly.

 
7. Run checkerboard test in serial
--------------------------------------

To run the checkboard test type::

        sfclean ; sfrun

within ``/home/tests/checkers``.

For now, the inversion will run only a single event on only a single processor.  Once we verify that everything is working correctly in thise case, we can move on to multiple events and multiple processors by modifying ``parameters.py``.



8. Run checkerboard test in parallel
-----------------------------------------
On a laptop or desktop with multiple cores, the work of an inversion can be carried out in parallel.  To run the checkerboard example in parallel over events (that is, with multiple event simulations running at the same on different cores), make the following changes to ``parameters.py``:

- to invert all available events instead of just one event, change ``NTASK`` from ``1`` to ``25``
- change the ``SYSTEM`` entry from ``serial`` to ``parallel``
- change the ``NPROCMAX`` entry to the number of cores available on your machine.

Besides running in parallel over events, the work of an individual event simulation can be parallelized over model regions. See the SPECFEM3D user manual for more information. Both parallelization over events and over model regions can be used at the same time under SeisFlows.  The current example, however, only illustrates the former.

Besides ``serial`` and ``parallel`` settings for running SeisFlows on laptops and desktops, there are also ``pbs``, ``slurm``, ``lsf`` options for running on computer clusters. See `here <http://seisflows.readthedocs.org/en/latest/manual/manual.html#system-configuration>`_ for more information.

