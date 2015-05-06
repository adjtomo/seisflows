
from seisflows.tools.code import Struct
globals()['PBS'] = Struct()
globals()['SLURM'] = Struct()
globals()['SPECFEM'] = Struct()


###

SPECFEM.ReceiverError = """
ERROR READING RECEIVERS

    Error reading receivers.

"""


SPECFEM.SourceError = """
ERROR READING SOURCES

    In DIRECTORY, there must be one or more files matching WILDCARD.

    DIRECTORY:  "%s"
    WILDCARD:  "%s"

"""


SPECFEM.ParameterWarning = """
PARAMETER WARNING

    There is a conflict between parameters.

    SPECFEM Parameter:  "%s"
    Old Value:  %s
    Overwriting with:  %s

"""


###

PBS.TaskError = """
SYSTEM ERROR

    Task failed:  %s.%s

    For more information, see output.pbs/%s

    Stopping workflow...

"""


SLURM.TaskError = """
SYSTEM ERROR

    Task failed:  %s.%s

    For more information, see output.slurm/%s

    Stopping workflow...

"""


