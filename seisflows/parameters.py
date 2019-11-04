"""
SEISFLOWS MODULES
"""
WORKFLOW='thrifty_inversion_nz'      # inversion, migration, modeling
SOLVER='specfem3d_nz'               # specfem2d, specfem3d
SYSTEM='maui_lg'                    # serial, pbs, slurm
OPTIMIZE='LBFGS'                    # steepest_descent, LBFGS, NLCG
LINESEARCH='Backtrack'              # Bracket, Backtrack
PREPROCESS='base'                   # base
POSTPROCESS='base'                  # base

"""
SIMULATION PARAMETERS
"""
CASE='Synthetic'        # 'Synthetic' for syn-syn, 'Data' for syn-data
MATERIALS='Elastic'        # 'Elastic': Vp and Vs, 'Acoustic': Vp only
DENSITY='Constant'         # If 'Variable', density will be updated
PRECOND=None

"""
WORKFLOW
"""
BEGIN=1                 # first iteration
END=1                   # last iteration
NSRC=2                  # number of sources

"""
SYSTEM
"""
NTASK=NSRC             # number of tasks
NPROC=80               # number of processers
NODESIZE=40            # number of cores per node (set by system)
WALLTIME=180           # master job walltime
TASKTIME=25            # maximum job time for each slave job 

# 'maui_lg' SYSTEM
ACCOUNT='nesi00263'                # NeSI account name 
MAIN_CLUSTER='maui'                # cluster to run simulations on
MAIN_PARTITION='nesi_research'     # partition of simulation cluster
ANCIL_CLUSTER='maui_ancil'         # cluster to run data processing on
ANCIL_PARTITION='nesi_prepost'     # partition of processing cluster
ANCIL_TASKTIME=20                  # for shorter tasktimes (default=TASKTIME)
NODES=2                            # number of nodes to occupy on cluster
CPUS_PER_TASK=1                    # available for multithreading (default=1)
SLURMARGS='--hint=nomultithread'   # request the entire node, requires 4 NODE
WITH_OPENMP=False                  # use openMP for job submission

"""
PREPROCESSING
"""
FORMAT='ascii'          # data file format

# 'BASE' PREPROCESS (ONLY FOR SYNTHETIC/SYNTHETIC EXAMPLES)
CHANNELS='z'            # data channels to be used
NORMALIZE=0             # normalize tracesi
FILTER=''               # highpass, lowpass, bandpass
MUTE=0                  # mute direct arrival
MUTECONST=0.            # mute constant (for muting early arrivals)
MUTESLOPE=0.            # mute slope (for muting early arrivals)

"""
PYATOA - REQUIRES EXTERNAL PACKAGE, NOT USED BY SEISFLOWS
"""
# PREPROCESSING
SYNTHETIC_UNIT = "disp"  # what units are Specfem3D synthetics; [disp, vel, acc]
MIN_PERIOD = 10.         # minimum period to filter all data in seconds
MAX_PERIOD = 30.         # maximum period
FILTER_CORNERS = 4       # corners for bandpass filter, default=4
ROTATE_TO_RTZ = False    # rotate waveforms to radial, transverse, vertical
FIX_WINDOWS = False      # fix misfit windows between iterations

# PYFLEX, PYADJOINT
PYFLEX_MAP = "hikurangi"        # config map for Pyflex; see pyatoa.plugins
ADJ_SRC_TYPE = "cc_hikurangi"   # misfit function to use in Pyadjoint
WINDOW_AMP_RATIO = 0.2          # extra amplitude threshold criteria; [0:1]

# OUTPUT FILES
SNAPSHOT = True             # periodically copy .h5 files incase corruption
WRITE_MISFIT_JSON = True    # save misfit information into a .json file 
CREATE_SRCRCV_VTK = True    # generate .vtk files for srcrcv locations

# PLOTTING
PLOT_WAVEFORMS = True       # plot waveform figures during workflow
PLOT_SRCRCV_MAPS = True     # plot source receiver maps during workflow
PLOT_MISFIT_MAPS = False    # plot misfit maps during workflow
COMBINE_IMGS = True         # combine waveforms and srcrcv maps into a .pdf
PURGE_WAVEFORMS = True      # if COMBINE_IMGS; delete waveforms after combine
PURGE_TILES = True          # if COMBINE_IMGS; delete intermediate wave/map imgs

"""
POSTPROCESSING
"""
TASKTIME_SMOOTH=3      # scales SYSTEM.TASKTIME for longer smoothing, default=1
SMOOTH_H=20000.          # smoothing Gaussian std. in horizontal, meters
SMOOTH_V=2000.          # smoothing Gaussian std. in vertical, meters
SCALE=1.                # scaling factor

""" 
OPTIMIZATION
"""
STEPCOUNTMAX=7         # maximum allowable trial step lengths (default=10, 
                        # minimum=3)
STEPINIT=0.25           # <<< This doesn't exist? step length safeguard
STEPFACTOR=0.75         # <<< This doesn't exist?
STEPLENINIT=0.05        # initial step, fraction of current model (default=0.05)
STEPLENMAX=0.2          # max step, fraction of current model (default=0.5)

# 'LBFGS' OPTIMIZATION (defaults if not set)
# LBFGSMAX=''           # periodic restart invterval (default=infinity)
# LBFGSMEM=''           # LBFGS memory, used for restarts (default=3)
# LBFGSTHRES=''         # descent direction threshold (default=0.0)

"""
SOLVER
"""
NT=10000               # number of time steps
DT=0.03                # time step
F0=.1                  # dominant frequency (SPECFEM2D only)

"""
SAVE OUTPUTS
"""
SAVEGRADIENT=False
SAVEKERNELS=False
SAVEMODEL=True
SAVE_AS="vector"      # save as "vector" (.npy) or "binary" (.bin) or "both"
SAVERESIDUALS=False
SAVETRACES=False

"""
CREATE VTKS FOR VISUALIZATION
False if not saved, otherwise list of quantities to combine
"""
CREATE_EVENT_KERNEL_VTK=["vp_kernel", "vs_kernel"]
CREATE_SUM_NOSMOOTH_KERNEL_VTK=["vp_kernel", "vs_kernel"]
CREATE_GRADIENT_KERNEL_VTK=["vp_kernel", "vs_kernel"]

"""
ADVANCED
"""
ENVIRONS=''
SOLVERIO='fortran_binary'

