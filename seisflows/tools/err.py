#!/usr/bin/env python
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

