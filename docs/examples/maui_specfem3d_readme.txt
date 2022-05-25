SEISFLOWS3 DEVEL VERSION INVERSION EXAMPLE README FOR MAUI
STARTED: 27/08/2021
LAST UPDATED: 17/11/2021
CONTACT: Bryant Chow (bhchow@alaska.edu)

################################################################################

INTRODUCTION

################################################################################

This README is an introduction to the SeisFlows3 (seisflows) package. 
It discusses the command line interface (CLI), the structure of a seisflows
working directory, and runs a basic inversion.

The steps cannot be run out of order, but some are optional, which means they 
provide information but contain no critical commands to run.

################################################################################

STEP 1: Environment setup

################################################################################

(1) SeisFlows needs to be run in a Conda environment. This must be done on the
Maui ancillary node. To get there you need to SSH from maui using:

    $ ssh -X w-mauivlab01.maui.nesi.org.nz

(2) We will use my pre-built Conda environment to run SeisFlows3. We first need
to load the Anaconda module on slurm.

    $ module load Anaconda3/5.2.0-GCC-7.1.0
    $ source activate /nesi/project/gns03247/PyPackages/conda_envs/seisflows

NOTE: You can also create your own conda environment using the text file 
provided which lists all the packages in my environment. 

Do not run the following 2 commands, they're just for informative purposes.
    
    $ cd /scale_wlg_nobackup/filesets/nobackup/gns03247/bchow/seisflows/examples/TUTORIAL/misc/conda_envs
    $ conda create --name seisflows seisflows_conda_env.txt

(3) You should now be able to access the CLI by typing 'seisflows' in the terminal
which will list out the available commands in the Seisflows package

    $ seisflows

usage: seisflows [-h] [-w [WORKDIR]] [-p [PARAMETER_FILE]]
                 [--path_file [PATH_FILE]]
                 {setup,modules,configure,init,submit,resume,restart,clean,par,check,convert,reset,inspect,debug,edit}
                 ...

================================================================================

                     SeisFlows: Waveform Inversion Package

================================================================================

optional arguments:
  -h, --help            show this help message and exit
  -w [WORKDIR], --workdir [WORKDIR]
                        The SeisFlows working directory, default: cwd
  -p [PARAMETER_FILE], --parameter_file [PARAMETER_FILE]
                        Parameters file, default: 'parameters.yaml'
  --path_file [PATH_FILE]
                        Legacy path file, default: 'paths.py'

command:
  Available SeisFlows arguments and their intended usages

    setup               Setup working directory from scratch
    modules             Print available module names
    configure           Fill parameter file with defaults
    init                Initiate working environment
    submit              Submit initial workflow to system
    resume              Re-submit previous workflow to system
    restart             Remove current environment and submit new workflow
    clean               Remove active working environment
    par                 View and edit parameter file
    check               Check state of an active environment
    convert             Convert model file format
    reset               Reset underlying machinery
    inspect             View inheritenace and ownership
    debug               Start interactive debug environment
    edit                Open source code file in text editor

'seisflows [command] -h' for more detailed descriptions of each command.


################################################################################

STEP 2: Directory setup

################################################################################

(1) We can use the CLI to create a new SF3 working directory. Let's make our own 
(empty) directory to follow along. 

   $ cd /scale_wlg_nobackup/filesets/nobackup/gns03247/bchow/seisflows/tutorial/
   $ mkdir work
   $ cd work
   $ seisflows setup

(2) A file named 'parameters.yaml' (par file) has been created. This is an 
ASCII file that will be your main point of control over SF3. 

(3) If you open 'parameters.yaml' you will see a few modules which define the 
main components of SF3. 

    $ cat parameters.yaml

Each of these modules controls a separate functionality of SF3, for example the
'preprocess' module controls how observed and synthetic waveforms are compared.

To view the available module choices defined by SF3, run the modules command
    
    $ seisflows modules

NOTE: You can see that each module is split into seisflows and seisflows-super.
The 'seisflows' modules are the base or default package modules, the
'seisflows-super' modules are more specific packages which supercede the base
modules. For example when running SF3 on Maui, we choose the 'maui' module
which inherits attributes from the 'slurm_lg' system.

(4) Because we're working on Maui and interested in performing an inversion 
we need to change some of the default module choices in the par file.

We can use the 'par' command which just directly edits the parameter file 

    $ seisflows par -h  # help command to see what it does
    $ seisflows par system maui
    $ seisflows par preprocess pyatoa

(4) Now run 'configure' to auto fill in the par file based on our module choices.

    $ seisflows configure

(5) Open up the 'parameters.yaml' file again to see the completed par file.
Commented areas contain docstrings explaining what each of the parameters does. 

NOTE: Parameters that need to be filled in by the user are denoted with
 '!!! REQUIRED PARAMETER !!!'

    $ cat parameters.yaml

(6) To save some time, I have pre-filled a parameters.yaml file, you can see 
which values have been set to what by using vimdiff. NOTE :q to quit vimdiff

    $ vimdiff parameters.yaml ../utils/parameters.yaml

(7) The utils directory which contains this parameters.yaml file has material
we'll need to run our inversion, so lets copy it into our working directory and 
then use the pre-filled parameter file.

    $ rsync -av ../utils . 
    $ cp -f utils/parameters.yaml .


################################################################################

STEP 3: DATA AND SPECFEM3D WORKING DIRECTORY (OPTIONAL)

################################################################################

* This step explains how to setup the DATA directory for SeisFlows3,
which tells SF3 and SPECFEM3D what sources and stations need to be included
in the inversion. 

* When running earthquake forward simulations, SPECFEM3D requires a DATA/ 
directory which contains (among other things):

    (a) CMTSOLUTION
    (b) STATIONS
    (c) Par_file

* Similarly, SeisFlows3 requires a DATA/ directory with the same files.
These files are distributed by SeisFlows3 to each of the SPECFEM3D working 
directories under the hood.

* One key difference is that SeisFlows3 requires N CMTSOLUTION files,
where N is the number of sources you want to simulate in parallel. N should
also match the parameters.yaml value: NTASK

You can check NTASK with the following command

    $ seisflows par NTASK

* Each CMTSOLUTION file should be different, and have a unique suffix. E.g., use 
earthquake IDs, or numerical ordering (e.g., CMTSOLUTION_001, CMTSOLUTION_002, 
..., CMTSOLUTION_N)

(1) The data, and true and initial models for this problem are already generated

    $ ls -l utils/specfem3d

In this example we will perform a synthetic inversion starting from a 
homogeneous halfspace (hh) initial model, and trying to recover a layered
halfspace true model

(2) The DATA/ directory already contains all the files (a) -- (c)
    
    $ ls -l utils/specfem3d/DATA

NOTE: Whatever is in this DATA/ directory will get copied into all the SeisFlows 
working so don't keep large files here (e.g., tomo files). You can check the 
location of the data directory with

    $ seisflows par specfem_data

(3) Within the DATA/ directory we also have the mesh files which defines a
North Island, NZ domain with ~30k elements.

    $ ls -l utils/specfem3d/DATA/meshfem3D_files

(4) Also provided is the DATABASES_MPI/ directory which defines the starting
velocity model in terms of the Specfem3D binary (.bin) files. 

    $ ls -l utils/specfem3d/OUTPUT_FILES_HH/DATABASES_MPI

The layered halfspace has a similar directory

    $ ls -l utils/specfem3d/OUTPUT_FILES_LAYERED/DATABASES_MPI


################################################################################

STEP 4: INITATE SEISFLOWS3 AND DEBUG MODE (OPTIONAL)

################################################################################

* We first want to initiate a SeisFlows3 working directory to have a look at 
an SF3 working state

(1) The init command starts a SeisFlows workflow but does not submit any jobs.

    $ sf init
    $ ls output

(2) 'init' created the output/ directory which contains .p (pickle) files, and 
.json (ascii) files. These files represent the state of the SeisFlows workflow. 

You can see the parameter .json file matches our input yaml file

    $ cat output/seisflows_parameters.json

* One of the most important tools in SeisFlows is its debug mode. This 
interactive IPython environment lets you explore an active workflow/

(3) Activate the debug mode, and when prompted type 'n' to get to IPython

    $ sf debug
    > n

(4) In this IPython environment, the entire SeisFlows architecture is loaded 
into memory. We won't do anything with it yet, but to demonstrate this we can 
look at the PAR and PATH variables

    > print(PAR)   # parameters are stored in a dict object
    > print(PATH)  # paths are also stored as dicts

* We can also access individual SeisFlows modules by calling their generic name

    In [6]: workflow
    Out[6]: <seisflows.workflow.inversion.Inversion at 0x2aab3bb20550>

(5) We will use debug mode more in later tutorials. To exit just type exit() at 
each IPython prompt

    In [7]: exit()
    > exit()

(6) So that we have a blank slate to run our inversion, we will run the 'clean'
command to get rid of all the created SeisFlows components

    $ sf clean
    > y

################################################################################

STEP 6: RUN THE INVERSION

################################################################################

* SeisFlows is an automated tool, which means that once we have set it up 
(parameters.yaml), we can just hit run and watch the output.

(1) We can run SeisFlows with the 'submit' command. When prompted, type 'y'
to submit the master job

    $ sf submit
    $ y

(2) Now that we have submitted the master job, we can check that it is running 
on the cluster

    $ squeue --account=gns03247 --cluster=maui,maui_ancil

(3) We can also track the progress of the inversion by watching the output
logs. NOTE: ctrl+c to quit the tail command

    $ tail -f *log

(4) The inversion will likely take some time (<.5 hr) but should proceed 
automatically until it completes one iteration. 

While it runs, we can have a look at a completed example

    $ cd ../example

################################################################################

STEP 7: UNDERSTANDING THE SEISFLOWS WORKING DIRECTORY

################################################################################

* Once your inversion starts running, it will create a directory structure
to keep track of the workflow and all the data generated. It should look
something like this:

(seisflows) [ancil] chowbr@w-mauivlab01 [example] $ ls
error-14907303.log  output/              output.logs/  output.stats/     scratch/
logs.old/           output-14907303.log  output.optim  parameters.yaml@  utils

Descriptions of the file structure:

(a) output-*.log
    ASCII text file where SeisFlows will print out status updates as it step 
    through the workflow. The number is the cluster job number

(b) error-*.log
    ASCII text file where SeisFlows will print traceback information if
    a workflow exits unexpectedly. If your job crashes, look here to figure out
    why

(c) output.optim
    ASCII text file where SeisFlows will print out information about the
    optimization algorithm, including iteration, step count and misfit info.

(d) output.logs/
    When individual processes are run (e.g., a simulation), they produce their
    own output files which are written here. Files are named after job numbers

(e) logs.old/
    If you re-run a workflow (e.g., using $ sf resume), new log files (a and b)
    will be created. The old ones are moved here so that previous run 
    information is retained

(f) output.stats/
    List-like stats information related to the optimization algorithm are
    written out here, to be used for e.g., plotting or analysis

(g) output/
    Location where all permament results are stored. These include the pickle
    files and parameter files which store the state of the workflow.
    Similarly gradients, updated velocity models, waveforms etc. are stored 
    here if requested by the User.

(h) scratch/
    The location of all the large, temporary working files required to run the 
    workflow. Each of the SeisFlows modules has its own directory in scratch/

    Many of these files are cleaned after each iteration.

    (1) evalfunc 
    Whenever SF3 evalutes the misfit/objective function, it stores temporary
    files here, these include the model being evaluated, and waveform residuals

    (2) evalgrad
    Whenever SF3 evalutes the gradient, it stores temporary files here,
    including event and misfit kernels, and the model being evaluated

    (3) optimize
    SF3 stores optimization algorithm information here as numpy .npy files.
    These include the current model, perturbed model, search direction, 
    gradient, and misfit.

    (4) preprocess
    SF3 stores waveform information here. If using the 'Pyatoa' preprocess 
    module, Pyatoa will store ASDFDataSets, processing logs (ASCII), and
    waveform figures + maps as .pdf files.

    (5) solver
    Each individual event gets its own directory in the solver scratch directory
    Each of these directories is simply a SPECFEM3D working directory, 
    containing the model databases, specfem binaries, and DATA/ directory.

    The 'mainsolver' is simply the first alphabetical directory, and is 
    typically used for SPECFEM tasks that don't need to be run by all of the 
    processes (e.g., smoothing)

    (6) system
    Any system-dependent files go here, for SLURM systems this is empty

################################################################################

STEP 8: ANALYZING THE INVERSION

################################################################################

* Once the inversion finishes running we can look at the completed directory
to understand how the inversion progressed

    $ cd ../work

(1) First we look at the output log 

    $ cat output*log

We can see that the inversion ran a forward and adjoint simulation, computed a
search direction with L-BFGS, and performed a line search with 2 trial steps

Within the line search we can see that the model was perturbed in order to
reduce the misfit.

We have also stopped the inversion before the 'finalize' function could be run, 
which means we can look at the scratch/ directory in a state before some temp
files are removed

(2) If we check the evalgrad directory we can see what temporary files SF3
saves while calculating gradient information

    $ ls scratch/evalgrad

If we look closer we can see SF3 is saving the .bin kernel files. These are 
large files SF3 deletes them after each iteration

    $ ls scratch/evalgrad/gradient

(3) Pyatoa will save figures and preprocessing logs as text files. We can look 
at the logs here but you will need to use rsync to view the .pdfs on 
your own machine, or use a program like 'evince' when on Maui 
(since we are ssh'd from maui, maui ancil cannot access the display)

    $ cat scratch/preprocess/logs/i01s00_2013p617227.log

The Pyatoa log file is event dependent and contains processing information
for each station 

(4) Now we will finish the inversion by resuming from the 'finalize' function.

!!! WARNING: Once you complete this step, the evalfunc and evalgrad directories
will be deleted.

First we set the RESUME_FROM parameter to 'finalize' to tell SF3 that we want to 
continue with the same inversion from a specific function

    $ seisflows par resume_from finalize

We also need to remove the STOP_AFTER parameter, which originally told SF3 to 
stop the inversion after the line search

    $ seisflows par stop_after ''

And finally we resume the workflow

    $ seisflows resume
    > y

This step should only take a few moments

(5) Once the master job stops running we can see that the final model from
iteration 1 has been saved to the output/ directory

    $ ls -l output

This model is saved in the SPECFEM3D binary (.bin) format and can be visualized
using .vtk files as normal. It is also stored as a vector in scratch/optimize
so if you'd like to continue running the inversion, this newly updated
model will be the starting model for the next tteration.


################################################################################

CONCLUSION

################################################################################

We have now run a single iteration from an SF3 inversion. Have a look around the 
directory and get a feel for where files are placed and how things are 
organized, this file structure is mostly universal for SF3. 

Some things you can do now are play around with the true and starting model,
the number of earthquakes in the inversion, or the types of preprocessing 
used.
