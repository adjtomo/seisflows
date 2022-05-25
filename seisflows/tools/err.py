#!/usr/bin/env python3
"""
Custom errors for Seisflows
"""


class ParameterError(ValueError):
    """
    A new ValueError class which explains the Parameter's that threw the error
    """
    def __init__(self, *args):
        if len(args) == 0:
            msg = "Bad parameter."
            super(ParameterError, self).__init__(msg)
        elif len(args) == 1:
            msg = f"Bad parameter: {args[0]}"
            super(ParameterError, self).__init__(msg)
        elif args[1] not in args[0]:
            msg = f"{args[1]} is not defined."
            super(ParameterError, self).__init__(msg)
        elif key in obj:
            msg = f"{args[0]} has bad value: {args[1].__getattr__(args[0])}"
            super(ParameterError, self).__init__(msg)

class CheckError(ValueError):
    """
    An error called by the Check functions within each module, that returns the
    name of the class that raised the error, as well as the parameter in
    question.
    """
    def __init__(self, cls, par):
        """
        CheckError simply returns a print message
        """
        msg = f"{cls.__class__.__name__} requires parameter {par}"
        super(CheckError, self).__init__(msg)

