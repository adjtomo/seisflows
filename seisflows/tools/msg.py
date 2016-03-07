

###

WarningOverwrite = """

WARNING: Data from previous workflow found in working directory.

To delete data and start a new workflow type:
  sfclean; sfrun

To resume existing workflow type:
  sfresume
"""


ReceiverError_SPECFEM = """

ERROR READING RECEIVERS

    Error reading receivers.

"""


SourceError_SPECFEM = """

ERROR READING SOURCES

    In DIRECTORY, there must be one or more files matching WILDCARD.

    DIRECTORY:  "%s"
    WILDCARD:  "%s"

"""


ParameterWarning_SPECFEM = """

PARAMETER WARNING

    There is a conflict between parameters.

    SPECFEM Parameter:  "%s"
    Old Value:  %s
    Overwriting with:  %s

"""


DataFormatWarning = """

DATA FORMAT WARNING

    reader format: %s
    writer format: %s

    Incompatible file formats may result in job failure or other problems.

"""

ReaderError = """

READER ERROR

   Seismic data reader not found.

   PAR.READER must correspond to an entry in seisflows.seistools.readers

"""


WriterError = """

WRITER ERROR

   Seismic data writer not found.

   PAR.WRITER must correspond to an entry in seisflows.seistools.writers

"""



###

TaskError_LSF = """

TASK ERROR

    Task failed:  %s.%s

    For more information, see output.lsf/%s

    Stopping workflow...

"""


TaskError_PBS = """

TASK ERROR

    Task failed:  %s.%s

    For more information, see output.pbs/%s

    Stopping workflow...

"""


TaskError_SLURM = """

TASK ERROR

    Task failed:  %s.%s

    For more information, see output.slurm/%s

    Stopping workflow...

"""


###

obspyImportError = """

DEPENDENCY ERROR

    The current data processing workflow requires OBSPY.  Please install it and
    try again, or consider using legacy data processing class 
    'preprocess.legacy', which lacks this dependency.

"""

mpiError1 = """

DEPENDENCY ERROR

    The following system configuration requires MPI4PY:

        system.%s

    Please install MPI4PY and try again, or consider choosing a different system
    configuration.

"""


mpiError2 = """

SYSTEM CONFIGURATION ERROR

    The following system configuration can be used only with single-core
    solvers:

        system.%s

    If your solver does in fact run on a single core, then set NPROC <- 1.

    If your solver requires multiple cores, then consider using lsf_lg, pbs_lg,
    or slurm_lg system interfaces instead.

"""


###

ImportError1 = """

SEISFLOWS IMPORT ERROR

    Please check that "custom_import" utility is being used as follows:

        custom_import(name1, name2)

    The resulting full dotted name "seisflows.name1.name2" must correspond to a
    module in the SeisFlows package.

"""


ImportError2 = """

SEISFLOWS IMPORT ERROR

    custom_import(name1, name2)

    Please check that "name1" is one of the following

        workflow
        solver
        optimize
        preprocess
        postprocess
        system

"""


ImportError3 = """

SEISFLOWS IMPORT ERROR

    The following module was not found in the SeisFlows package:

        seisflows.%s.%s

    Please check user-supplied %s parameter.

"""


ImportError4 = """

SEISFLOWS IMPORT ERROR

    By convention, SeisFlows module 

        seisflows.%s.%s

    must contain a class named

        %s

"""

