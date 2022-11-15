Running on a Cluster (UAF Chinook)  
==================================

This docs page introduces Users to running research-grade problems on high 
performance computers (HPC). It is currently targeted at a specific university
cluster but may be expanded to other systems as necessary. Instruction should 
hopefully be generalizable to other clusters, although Users may need to 
`write their own custom interface <extending.html>`__. 

`Chinook <https://uaf-rcs.gitbook.io/uaf-rcs-hpc-docs/hpc#chinook>`__ is 
University of Alaska Fairbank’s (UAF) high performance computer. Chinook 
is an Intel machine running the SLURM workload manager and Rocky Linux 8. 
Chinook is operated by Research Computing Systems (RCS).

.. note:: 

    These instructions were written to be followed along during a group meeting 
    at UAF and therefore go into some minute details that may not be relevant 
    for all.

.. note::
    
    Resources last accessed Nov. 14, 2022

0) Access Chinook
~~~~~~~~~~~~~~~~~

For those following along in-person, we will access Chinook via SSH and then 
access the updated chinook by SSH'ing into Chinook04 which contains the 
latest OS update for Chinook. You should be met by this cool fish upon 
successful login to Chinook04:

.. parsed-literal:: 

           /`-._
         _/,.._/          dP""b8 88  88 88 88b 88  dP"Yb   dP"Yb  88  dP
      ,-'   ,  `-:,.-')  dP   `" 88  88 88 88Yb88 dP   Yb dP   Yb 88odP  
     : o ):';     _  {   Yb      888888 88 88 Y88 Yb   dP Yb   dP 88"Yb  
      `-.  `' _,.-\`-.)   YboodP 88  88 88 88  Y8  YbodP   YbodP  88  Yb 
         `\\``\,.-'    



1) Install Conda
~~~~~~~~~~~~~~~~

First we need to install Conda, the Python package managaer. RCS has 
installation instructions related to installing Minoconda here:

https://uaf-rcs.gitbook.io/uaf-rcs-hpc-docs/third-party-software/miniconda

2) Install or load adjTomo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You have two options here for grabbing the adjTomo softwater suite. 1) Easiest 
would be to load a pre-installed Conda environment:

.. code:: bash
    
    conda activate /import/c1/ERTHQUAK/bhchow/REPOS/miniconda3/envs/adjtomo    

2) A more flexible solution would be to create your own Conda environment and 
install software yourself. For that you will have to follow the instructions on 
the `main docs page <index.html#installation>`__. 

If you go with Option 2, make sure you activate your conda environment before 
proceeding.

3) Set up a SPECFEM2D/3D/3D_GLOBE working directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SeisFlows requires a pre-established SPECFEM working directory, with:

a) binary executables (configured and compiled), 
b) an appropriate DATA/ directory containing source, stations and Par_file, and 
c) an initial (and target) model defined using one of SPECFEM's internally defined model formats.

.. note::
    
    If you already have a valid directory where you run forward simulations, 
    you can likely skip the following subsections 3a -- 3d.


3a) Configure and compile
`````````````````````````

You will need to clone SPECFEM2D/3D/3D_GLOBE (you choose the flavor), configure
and compile the code. Below are instructions specifically for Chinook. 

Other clusters will have different compiler options and requirements that are 
machine/OS specific so it is difficult to write a generalized set of 
instructions.


Here we choose SPECFEM3D and compile using the Intel compiler suite:

.. code:: bash
    
    mkdir $CENTER1/REPOS  # Center1 is our working filesystem
    cd $CENTER1/REPOS
    git clone --branch devel https://github.com/SPECFEM/specfem3d.git
    cd specfem3d
    module load intel  # latest Intel compiler suite
    ./configure F90=ifort FC=ifort MPIFC=mpiifort CC=icc MPICC=mpiicc --with-mpi 
    make all

3b) Generate appropriate DATA/ directory
``````````````````````````````````````````

Here you can choose to set your own mesh and model parameters to suit your 
research problem. For the sake of simplicity we will use the homogeneous 
halfspace model located in the EXAMPLES/ directory to generate our starting
model.

We will also work in a separate SPECFEM working directory (outside the cloned
repository) to keep things clean and manageable.


.. code:: bash

    mkdir $CENTER1/work/specfem3d_workdir  # clean working directory
    cd $CENTER1/work/specfem3d_workdir
    ln -s $CENTER1/REPOS/specfem3d/bin .  # making sure we have the executables
    cp -r $CENTER1/REPOS/specfem3d/EXAMPLES/homogeneous_halfspace/DATA .
    cp -r $CENTER1/REPOS/specfem3d/EXAMPLES/homogeneous_halfspace/meshfem3D_files ./DATA
    mkdir OUTPUT_FILES


3c) Dealing with multiple sources
`````````````````````````````````

One key difference that needs to be addressed is that SeisFlows requires sources
be tagged. For example, if you want to run 10 events in your inversion
you will need to individually tag each event with the appropriate format.

In SPECFEM3D our source prefix will be 'CMTSOLUTION'. If we have multiple 
CMTSOLUTIONS, then one easy way to differentiate them would be to name them e.g.: 
CMTSOLUTION_1, CMTSOLUTION_2, ..., CMTSOLUTION_N. These tags could also 
refer to event ids or origin times, it's up to the user.

`Here is one example of the naming schema used in a published study. 
<https://github.com/bch0w/spectral/tree/master/nzatom/cmtsolutions>`__

For this example, since we don't have multiple sources to choose from, we will
simply copy our example CMTSOLUTION and rename:

.. code:: bash

    cd $CENTER1/work/specfem3d_workdir/DATA
    mv CMTSOLUTION CMTSOLUTION_01  # source 1 is the example default 
    cp CMTSOLUTION_01 CMTSOLUTION_02  # source 2 is the same as source 1
    ln -s CMTSOLUTION_01/ CMTSOLUTION  # not necessary but aesthetically pleasing

3d) Create Initial (and target) models
``````````````````````````````````````

Now we'll run SPECFEM to generate our mesh and model. This is the same procedure 
you would follow if running a forward simulation in SPECFEM, except we will not
run the solver. 

We need a slurm-specific SBATCH script to run our executables. You can find `example SBATCH scripts for Chinook here <https://github.com/bch0w/simutils/blob/master/cluster/runscripts/chinook/specfem3d/>`__. I will use two files from this directory, `run_xmeshfem3d.sh` and `run_xgenerate_databases.sh`.

.. note::
    
    SPECFEM2D and SPECFEM3D_GLOBE do not require the `xgenerate_databases` step

.. code:: bash

    sbatch run_xmeshfem3d.sh  # generates mesh files
    sbatch run_xgenerate_databases.sh  # generates model files

By the end we want to have a number of binary (.bin) files that contain our
model. These should be located in the local path:  

.. code:: bash

    ls OUTPUT_FILES/DATABASES_MPI  # should contain vp, vs, and rho files

Finally, we need to set the `model` parameter in the SPECFEM Par_file to 'gll'.
This will tell future runs of SPECFEM to read the model we just created, 
rather than trying to define it from internal parameters:

.. code:: bash
    
    seisflows sempar -P DATA/Par_file model gll

Have a look at the `command line tool docs page <command_line_tool.html>`__ 
for more information on the command line tools available for SeisFlows.


4) Setting up a SeisFlows working directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We are now ready to run SeisFlows. We just have to set up a working directory
and point the parameter file at the correct locations such that SeisFlows can
find our SPECFEM working directory. 

I will run SeisFlows in a separate directory to keep things clean. 

.. code:: bash

    mkdir $CENTER1/work/seisflows_workdir
    cd $CENTER1/work/seisflows_workdir
    seisflows setup  # creates a template parameters.yaml file

Have a look at the `parameter file docs page <parameter_file.html>`__ for 
more information on how the file is structured.

4a) SeisFlows parameter file
```````````````````````````````

You can look at the generated parameter file to see what the template version 
looks like (using a text editor or cat). We will simply overwrite some of the
base starting parameters to suit our current use case. Use the ``seisflows par``
command to do this quickly on the command line.

SeisFlows already contains a pre-built Chinook interface (based on a general 
SLURM interface). You can use ``seisflows print modules`` to see all valid 
system (and other modules) choices. 

.. code:: bash

    seisflows print modules

If you do not see your own system (for non-Chinook users) supported, you will 
need to follow the instructions on 
`writing your own system-subclass <https://seisflows.readthedocs.io/en/devel/extending.html>`__

Here we overwrite some default parameters to set up the base modules for our 
workflow:

.. code:: bash

    seisflows par system chinook  # chinook system interface
    seisflows par solver specfem3d  # specfem3d cartesian version
    seisflows par preprocess null  # turn OFF preprocessing for now
    seisflows par optimize null  # turn OFF optimization 


By default we are running a ``forward`` workflow, which simply runs forward
simulations en-masse. In following sections we will swap over to an inversion
workflow.

4b) Configuring the parameter file
````````````````````````````````````

Each choice of base module (i.e., workflow system, solver, preprocess, optimize)
comes with it's own distinct set of parameters. SeisFlows therefore 
dynamically generates a parameter file based on User choices for the base 
modules and the appropriate source code doc strings. 

We can configure our parameter file with:

.. code:: bash

    seisflows configure

Have a look at your parameter file now to see all the module-specific parameters 
that have been instantiated.


4c) Checking the parameter file
`````````````````````````````````

As with SPECFEM, the parameter file in SeisFlows controls the entire package, 
and all the parameters that have been set using the ``seisflows configure`` 
command are applicable to your current workflow. 

.. warning::

    It is up to a prospetive user to carefully read and understand what each 
    parameter does. I have tried to make the docstrings as comprehensive as 
    possible, but things do slip through the cracks. If you find that a certain 
    parameter is not well explained, ambiguous, etc. please open up a GitHub 
    issue or PR with clarifying changes.

Each module in SeisFlows has a ``check`` function which it uses to determine
parameter validity. 

Users can use this ``check`` function to quickly determine missing,
inappropriate, or invalid parameters in their parameter file.

.. code:: bash

    seisflows check

You can use this method to fix parameters one by one until no errors are 
raised, after which you should be confident that you are able to run your 
workflow.

Following the parameter errors raised, you will have to change the following:

.. code:: bash

    # Changing paths to tell SeisFlows where to find SPECFEM
    seisflows par path_specfem_bin ${CENTER1}/work/specfem3d_workdir/bin
    seisflows par path_specfem_data ${CENTER1}/work/specfem3d_workdir/DATA
    seisflows par path_model_init ${CENTER1}/work/specfem3d_workdir/OUTPUT_FILES/DATABASES_MPI

Based on docstrings, I know I will also want to set the following parameters 
in order to suit my current research problem:

.. code:: bash

    # Changing parameters to suit our workflow
    seisflows par ntask 2  # two events, corresponding to two CMTSOLUTIONS
    seisflows par tasktime 5  # walltime for individual simulations
    seisflows par walltime 20  # walltime for the entire workflow
    seisflows par nproc 4  # to match the SPECFEM parameter of the same name
    seisflows par export_traces True  # save seismograms to disk 


5) Submit the main job
~~~~~~~~~~~~~~~~~~~~~~~~~

SeisFlows operates using a serial, single-core main job submitted to a 
compute node. This main job will act like `you`, the researcher:

Through the pre-defined Chinook/SLURM system interface, the main job already 
knows how to:

- submit jobs (using sbatch), 
- monitor the queue (using sacct)
- book keep SPECFEM and manage the filesystem
- stop jobs if any errors occur

To submit the main job, we simply run:

.. code:: bash

    seisflows submit

Now that we have submitted the workflow, the main job will run en-masse
forward simulations. In other words, it runs two forward simulations 
corresponding to the two CMTSOLUTIONS we have in our DATA/ directory.

.. note::

    On Chinook, in order to keep the main partition clean, all master jobs are 
    submitted to the 'debug' node by default. This is hardcoded into the Chinook 
    implementation. Future work may place the main job on the login node as well.


6) Inspecting SeisFlows
~~~~~~~~~~~~~~~~~~~~~~~~~~

Have a look at the `working directory docs page <working_directory.html>`__ 
for an explanation of the directories and files being generated.

Monitor the job queue to see the master job and all spawned compute jobs 
that get submitted to the system using the `squeue` or `sacct` commands.

- The main log is writing to ``sflog.txt``
- Each spawned job is logging to a unique file in ``logs/``
- Each source has it's own working directory in ``scratch/solver/``

a) Recovering from job failures
````````````````````````````````

SeisFlows has a state file (`sfstate.txt`) that tracks the progress of your 
inversion. Each main workflow function (e.g., forward simulations) constitute a
'checkpoint' in the workflow. If a function completes sucessfully, it is 
labeled 'completed'. Jobs which fail are labelled 'failed'.

If your job fails (e.g., due to walltime), you can simply run 
``seisflows submit`` again, and SeisFlows will know to skip over the already 
completed tasks, saving computational cost.

.. note::
    Currently, SeisFlows does not know how to track individually completed jobs. 
    E.g., for a two event workflow, one event completes a successful forward 
    simulation, but the other one fails for unknown reason. Currently SeisFlows 
    will need to re-run ALL forward simulations. In the future I hope to 
    include some more detailed checkpointing to avoid this.

b) SeisFlows debug mode
`````````````````````````

SeisFlows has a debug mode, which is simply an IPython environment with all
SeisFlows modules and parameters loaded. This allows the User to step through
code while debugging or developing. 

This is especially useful when you are looking at source code (trying to 
figure out a bug), and you want to know "what is this variable?", or 
"what does this function return?". You can figure that out with:

.. code:: bash
    
    seisflows debug


7) Modifying for a synthetic inversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Great! This is essentially the standard method of operating SeisFlows: 
manually setting up your SPECFEM directory, tooling the parameter file, and
submitting your job.

But what if you now want to run a synthetic inversion to compare synthetic 
seismograms from two very similar models? How do you get from here to there?

.. note::
        
    It is a good idea to either clear out your current working directory, or
    start a new one, before proceeding with a separate workflow. To delete all
    non-essential files, you can run `seisflows clean -f`. Otherwise to start
    a new working directory, you can simply copy over the parameter file to
    a new directory. 

.. note::
     
    If you decide to copy over the parameter file (from previous note), make 
    sure you update your paths! 


a) Swap modules in the parameter file
``````````````````````````````````````

SeisFlows ``swap`` allows Users to swap out valid modules without disturbing 
the remainder of the parameter file. That is, if we wanted to swap out 
our 'forward' workflow for an 'inversion' workflow, we can do:

.. code:: bash

    seisflows swap workflow inversion

If you look at your parameter file now, you will see a suite of new parameters
that control an inversion workflow.

This is the same for swapping from SPECFEM3D -> SPECFEM3D_GLOBE or choosing 
preprocessing parameters.

The inversion workflow requires a corresponding `preprocess` and `optimize` 
module. We can set these to the preferred classes `pyaflowa` and `LBFGS`. Again
have a look at the output of `seisflows print modules` for all choices.

.. code:: bash
    
    seisflows swap preprocess pyaflowa
    seisflows swap optimize LBFGS


b) Generate your target model
````````````````````````````````

The inversion workflow requires data. Since we have decided to do a synthetic
inversion, SeisFlows requires a target model. If we were doing a real-data
inversion, SeisFlows would require waveform data.

We'll set up our target model as a slightly altered homogeneous halfspace model
to keep things simple:

.. code:: bash

    cd $CENTER1/work/specfem3d_workdir
    mv OUTPUT_FILES OUTPUT_FILES_INIT  # setting aside our initial model
    cd DATA/meshfem3D_files
    mv Mesh_Par_file Mesh_Par_file_init  # setting aside initial mesh
    cp Mesh_Par_file_init Mesh_Par_file_true
    ln -s Mesh_Par_file_true Mesh_Par_file  # ensuring mesh name is correct
    
Here you need to manually: 

1) open up the `Mesh_Par_file` file, 
2) scroll down to the `'Domain materials'` section (around Line 86) and 
3) edit the material parameters to your choosing.

I will increase velocities by 10%, that is Vp: 2800 -> 3020 m/s and Vs: 
1500 -> 1650 m/s.

And now we need to run the SPECFEM binaries again to generate our target model

.. code:: bash

    cd $CENTER1/work/specfem3d_workdir
    mkdir OUTPUT_FILES_TRUE
    ln -s OUTPUT_FILES_TRUE OUTPUT_FILES
    seisflows sempar -P DATA/Par_file model default  # make sure SPECFEM reads the model from the mesh
    sbatch run_xmeshfem3D.sh
    sbatch run_xgenerate_databases.sh
    seisflows sempar -P DATA/Par_file model gll  # reset for seisflows run


c) Set inversion-specific parameters
`````````````````````````````````````

Again we can use `seisflows check` to see what new parameters we need to set, 
which are introduced by the 3 new modules we have (workflow, preprocess, 
optimize).

.. code:: bash
    
    seisflows check

Following the 'check'list we will need to change the folowing parameters

.. code:: bash

    seisflows par data_case synthetic  # synthetic inversion (no data)
    seisflows par path_model_true ${CENTER1}/work/specfem3d_workdir/OUTPUT_FILES_TRUE/DATABASES_MPI


We'll also set the following parameters:

.. code:: bash

    seisflows par path_model_init ${CENTER1}/work/specfem3d_workdir/OUTPUT_FILES_INIT/DATABASES_MPI  # to deal with the fact that we renamed this directory
    seisflows par materials elastic  # update vp and vs
    seisflows par end 2  # stop after iteration 2 is finished

d) SeisFlows submit
````````````````````

Again we run `submit` to submit our workflow. Make sure you have cleared out 
your previous run or switched to a new working directory!

.. code:: bash

    seisflows submit
 
