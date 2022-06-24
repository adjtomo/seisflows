#!/usr/bin/env python3
"""
Core class definitions for SeisFlows. Defines some unique class objects that
are used heavily during a Seisflows workflow.
"""
import sys
import logging


class Base(object):
    """
    Defines the core Base object for all SeisFlows modules. All modules MUST
    inherit from the Base object to work properly. This Base class essentially
    dictates the required structure of a SeisFlows class.
    """
    def __init__(self):
        """
        SeisFlows instantiates its required parameters through the
        SeisFlowsPathsParameters class, which scaffolds a rigid framework of
        how parameters and paths should be defined by the program. This is
        then used to build the parameter file dynamically.
        """
        self.required = SeisFlowsPathsParameters()

    def module(self, name):
        """
        Access globally stored SeisFlows modules located in sys.modules

        :rtype: Class or Dict or None
        :return: Returns a SeisFlows module or Dictionary containing paths or
            parameters. Else None if the chosen module has not been instantiated
        """
        try:
            mod = sys.modules[f"seisflows_{name}"]
        except KeyError:
            self.logger.warning(f"seisflows_{name} has not been instantiated")
            mod = None
        return mod

    @property
    def logger(self):
        """
        An instance specific logger which imprints inheritance information into
        the log statements, making it easier to debug functions with
        multiple inheritance
        """
        return logging.getLogger(
            self.__class__.__name__).getChild(self.__class__.__qualname__)

    @property
    def par(self):
        """
        Quick access SeisFlows parameters from sys.modules. Throws a warning
        if parameters have not been instantiated

        :rtype: Dict or None
        :return: Returns a Dictionary with instantiated parameters, or None if
            parameters have not been instantiated
        """
        return self.module("parameters")

    @property
    def path(self):
        """
        Quick access SeisFlows paths from sys.modules. Throws a warning
        if paths have not been instantiated

        :rtype: Dict or None
        :return: Returns a Dictionary with instantiated paths, or None if
            parameters have not been instantiated
        """
        return self.module("paths")

    def check(self, validate=True):
        """
        General check() function for each module to check the validity of the
        user-input parameters and paths
        """
        if validate:
            self.required.validate()

        # Example of a check statement
        # assert(self.par.PARAMETER == example_value), f"Parameter != example"

    def setup(self):
        """
        A placeholder function for any initialization or setup tasks that
        need to be run once at the beginning of any workflow.
        """
        return

    def finalize(self):
        """
        A placeholder function for any finalization or tear-down tasks that
        need to be run at the end of any iteration or workflow.
        """
        return


class Dict(dict):
    """
    A dictionary replacement which allows for easier parameter access through
    getting and setting attributes. Also has some functionality to make string
    printing prettier
    """
    def __str__(self):
        """Pretty print dictionaries and first level nested dictionaries"""
        str_ = ""
        try:
            longest_key = max([len(_) for _ in self.keys()])
            for key, val in self.items():
                str_ += f"{key:<{longest_key}}: {val}\n"
        except ValueError:
            pass
        return str_

    def __repr__(self):
        """Pretty print when calling an instance of this object"""
        return(self.__str__())

    def __getattr__(self, key):
        """Attribute-like access of the internal dictionary attributes"""
        try:
            return self[key]
        except KeyError:
            raise AttributeError(f"{key} not found in Dict")

    def __setattr__(self, key, val):
        """Setting attributes can only be performed one time"""
        self.__dict__[key] = val


class Null:
    """
    A null object that always and reliably does nothing
    """
    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, *args, **kwargs):
        return self

    def __nonzero__(self):
        return False

    def __getattr__(self, key):
        return self

    def __setattr__(self, key, val):
        return self

    def __delattr__(self, key):
        return self


class SeisFlowsPathsParameters:
    """
    A class used to simplify defining required or optional paths and parameters
    by enforcing a specific structure to their entry into the environment.

    .. note::
        if a path or parameter is optional it requires a default value.
    """
    default_par = "REQUIRED PARAMETER"
    default_path = "REQUIRED PATH"

    def __init__(self, base=None):
        """
        We simply store paths and parameters as nested dictioanries. Due to the
        use of inheritance, the class can be passed to itself on initialization
        which means paths and parameters can be adopted from base class

        :type base: seisflows.config.DefinePathsParameters
        :param base: paths and parameters from abstract Base class that need to
            be inherited by the current child class.
        """
        self.parameters, self.paths = {}, {}
        if base:
            self.parameters.update(base.parameters)
            self.paths.update(base.paths)

    def par(self, parameter, required, docstr, par_type, default=None):
        """
        Add a parameter to the internal list of parameters

        :type parameter: str
        :param parameter: name of the parameter
        :type required: bool
        :param required: whether or not the parameter is required. If it is not
            required, then a default value should be given
        :type docstr: str
        :param docstr: Short explanatory doc string that defines what the
            parameter is used for.
        :type par_type: class or str
        :param par_type: the parameter type, used for doc strings and also
            parameter validation
        :param default: default value for the parameter, can be any type
        """
        if required:
            default = self.default_par
        if type(par_type) == type:
            par_type = par_type.__name__
        self.parameters[parameter] = {"docstr": docstr, "required": required,
                                      "default": default, "type": par_type}

    def path(self, path, required, docstr, default=None):
        """
        Add a path to the internal list of paths

        :type path: str
        :param path: name of the parameter
        :type required: bool
        :param required: whether or not the path is required. If it is not
            required, then a default value should be given
        :type docstr: str
        :param docstr: Short explanatory doc string that defines what the
            path is used for.
        :type default: str
        :param default: default value for the path

        """
        if required:
            default = self.default_path
        self.paths[path] = {"docstr": docstr, "required": required,
                            "default": default}

    def validate(self, paths=True, parameters=True):
        """
        Set internal paths and parameter values into sys.modules. Should be
        called by each modules check() function.

        Ensures that required paths and parameters are set by the User in the
        parameter file and that default values are stored for any optional
        paths and parameters which are not explicitely set.

        :type paths: bool
        :param paths: validate the internal path values
        :type parameters: bool
        :param parameters: validate the internal parameter values
        :raises ParameterError: if a required path or parameter is not set by
            the user.
        """
        if paths:
            sys_path = sys.modules["seisflows_parameters"]
            for key, attrs in self.paths.items():
                if attrs["required"] and (key not in sys_path):
                    raise KeyError(f"{sys_path} has no key {key}")
                elif key not in sys_path:
                    setattr(sys_path, key, attrs["default"])

        if parameters:
            sys_par = sys.modules["seisflows_paths"]
            for key, attrs in self.parameters.items():
                if attrs["required"] and (key not in sys_par):
                    raise ValueError(sys_par, key)
                elif key not in sys_par:
                    setattr(sys_par, key, attrs["default"])

