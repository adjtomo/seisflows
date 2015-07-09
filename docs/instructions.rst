1. Download seisflows
---------------------

To run seisflows you'll need a Unix system with Python and standard utilities.  From the command line type::
 
        mkdir /home/packages
        cd /home/packages
        git clone https://github.com/PrincetonUniversity/seisflows.git

If you prefer a location other than ``/home/packages`` , modify the commands above and below accordingly.


2. Set environment variables
----------------------------

Add the following lines to ``/home/.bash_profile``, if you are using bash (or modify accordingly, if you are using a different shell)::

        export PATH=$PATH:/home/packages/seisflows/scripts
        export PYTHONPATH=$PYTHONPATH:/home/packages/seisflows
 

Don't forget to update any open shells::

        source /home/.bash_profile
 

 

3. Run "system" test
---------------------

 
Run the following test to make sure everything is working::

        cd /home/packages/seisflows/tests/integration/system
        ./clean; ./run


If an hello message is displayed, then the test was successful and it's alright to move on.

 

 

4. Run "optimization" test
--------------------------


Run the following test to make sure everything is working::

        cd /home/packages/seisflows/tests/integration/system
        ./clean; ./run


If the optimization problem is solved in 50 iterations or fewer, the test was successful and we can move on.

 

 

5. Configure and compile SPECFEM2D
----------------------------------

First, download SPECFEM2D from github (to avoid possible version differences, let us use version ``d745c542``)::

        cd /home/packages
        git clone --recursive --branch devel https://github.com/geodynamics/specfem2d.git specfem2d-d745c542
        cd specfem2d-d745c542
        git checkout d745c542


Next, configure and compile SPECFEM2D using ifort (preferred) or gfortran::

        cd /home/packages/specfem2d-d745c542
        ./configure FC=ifort
        make all
 
If there are no compilation errors, then proceed to the next step.


6. Set up FWI checkerboard test
-------------------------------

Download the starting model and other input files required for the full waveform inversion (FWI) checkerboard test.  For simplicity, let us assume these files will be placed in ``/home/tests`` (if you prefer a different location, then modify the following commands accordingly)::
 
        mkdir /home/tests/
        cd /home/tests/
        wget --recursive --no-parent --no-host-directories --cut-dirs=2 --reject "index.html*" http://tigress-web.princeton.edu/~rmodrak/Examples2d/


After the download completes, make sure that all paths specified in ``home/tests/checkers/paths.py``  are correct.  For example, if you compiled SPECFEM2D somewhere other than ``/home/packages/specfem2d-d745c542``, you will need to modify the entry SPECFEM2D binary file entry accordingly.

 
7. Run FWI checkerboard test in serial
--------------------------------------

To run the checkboard test, simply type::

        sfclean ; sfrun

within the directory ``/home/tests/checkers``.

For now, the inversion will run on the local host using only a single event and only a single processor.  Once we verify that everything is working correctly in this simple case, we can move to the case of multiple events and multiple processors by changing the 'system' setting in ``/home/tests/checkers/parameters.py``.



8. Run FWI checkerboard test in parrallel
-----------------------------------------
If you have access to a multicore laptop or desktop, SeisFlows can be used to carry out an inversion running multiple wavefield simulations in parallel.  To run the FWI checkerboard example in this manner, change the ``SYSTEM`` entry in ``parameters.py`` from ``serial`` to ``parallel`` and change the ``NPROCMAX`` entry to the number of processors available.
