

###

TaskError = """
TASK ERROR

    Task failed:  %s.%s

    For more information, see output.slurm/%s

    Stopping workflow...
"""


###

SourceInputError = """
SOURCE INPUT ERROR

    In DIRECTORY, there must be one or more files matching WILDCARD.

    DIRECTORY:  "%s"
    WILDCARD:  "%s"
"""


###

SolverParameterOverwriteWarning = """
SOLVER PARAMETER OVERWRITE WARNING

    There is a conflict between SeisFlows and solver parameters.

    Solver Parameter:  "%s"
    Old Value:  %s
    New Value:  %s
"""
