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
MISFIT='Waveform'
MATERIALS='Elastic'
DENSITY='Variable'
PRECOND=None

"""
WORKFLOW
"""
BEGIN=2                 # first iteration
END=2                   # last iteration
NREC=20                 # number of receivers
NSRC=2                  # number of sources

"""
SYSTEM
"""
NTASK=NSRC             # number of tasks
NPROC=144              # number of processers
NODESIZE=80            # number of cores per node (set by system)
WALLTIME=240           # master job walltime
TASKTIME=30            # maximum job time for each slave job 

# 'maui_lg' SYSTEM
ACCOUNT='nesi00263'                # NeSI account name 
MAIN_CLUSTER='maui'                # cluster to run simulations on
MAIN_PARTITION='nesi_research'     # partition of simulation cluster
ANCIL_CLUSTER='maui_ancil'         # cluster to run data processing on
ANCIL_PARTITION='nesi_prepost'     # partition of processing cluster
ANCIL_TASKTIME=int(NREC*0.2)       # for shorter tasktimes (default=TASKTIME)
# NODES=4
CPUS_PER_TASK=1                    # available for multithreading (default=1)
# SLURMARGS='--hint=nomultithread'   # request the entire node, requires 4 NODE


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
POSTPROCESSING
"""
SMOOTH_H=10000.          # smoothing radius in the horizontal, meters
SMOOTH_V=4000.          # smoothing radius in the vertical, meters
SCALE=1.                # scaling factor

""" 
OPTIMIZATION
"""
STEPCOUNTMAX=3         # maximum allowable trial step lengths (default=10, 
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
NT=30000               # number of time steps
DT=0.00775             # time step
F0=.1                  # dominant frequency (SPECFEM2D only)

"""
SAVE OUTPUTS (1 TO SAVE, 0 TO DISCARD)
"""
SAVEGRADIENT=0
SAVEKERNELS=0
SAVEMODEL=1
SAVERESIDUALS=0
SAVETRACES=0

"""
ADVANCED
"""
ENVIRONS=''
SOLVERIO='fortran_binary'

