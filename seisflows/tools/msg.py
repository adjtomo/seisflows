
from seisflows.tools.code import Struct
globals()['PBS'] = Struct()
globals()['SLURM'] = Struct()
globals()['SPECFEM'] = Struct()


###

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


###

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

