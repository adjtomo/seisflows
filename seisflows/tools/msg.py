
from seisflows.tools.code import Struct
globals()['PBS'] = Struct()
globals()['SLURM'] = Struct()
globals()['SPECFEM'] = Struct()


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

