
from seisflows.tools.code import Struct


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


obspyImportError = """
OBSPY IMPORT ERROR

    The current workflow requires obspy.  Install it and try again, or consider 
    using legacy data processing class 'preprocess.legacy', which lacks this
    requirement.

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

