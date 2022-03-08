#!/usr/bin/env python
"""
!!! This superclass is work in progress !!!

This is the subclass seisflows.solver.specfem3d_globe
This class provides utilities for the Seisflows solver interactions with
Specfem3D Globe. It inherits all attributes from seisflows3.solver.Base,
and overwrites these functions to provide specified interaction with Specfem3D.
"""
import os
import sys
import logging
from glob import glob

import seisflows3.plugins.solver.specfem3d_globe as solvertools
from seisflows3.tools.seismic import getpar, setpar, Model, Minmax
from seisflows3.plugins.io import loadbypar, copybin, loadbin, savebin
from seisflows3.tools import unix
from seisflows3.tools.seismic import call_solver
from seisflows3.tools.tools import Struct, exists
from seisflows3.config import (ParameterError, custom_import,
                              SeisFlowsPathsParameters)

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]

system = sys.modules["seisflows_system"]


class Specfem3DGlobe(custom_import("solver", "specfem3d")):
    """
    Python interface to Specfem3D Cartesian. This subclass inherits functions
    from seisflows3.solver.specfem3d.Specfem3D

    !!! See base class for method descriptions !!!
    """
    # Class-specific logger accessed using self.logger
    logger = logging.getLogger(__name__).getChild(__qualname__)

    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters(super().required)

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths on top of the base class checks.
        """
        if validate:
            self.required.validate()
        super().check(validate=False)

        assert(PAR.MATERIALS.upper() in ["ISOTROPIC", "ANISOTROPIC"])

        # Overwrite the base class parameters based on the material choice
        self.parameters = []
        if PAR.MATERIALS.upper() == "ISOTROPIC":
            self.parameters += ["vp", "vs"]
        elif PAR.MATERIALS.upper() == "ANISOTROPIC":
            self.parameters += ["vpv", "vph", "vsv", "vsh", "eta"]

    def load(self, path, prefix='reg1_', suffix='', verbose=False):
        """
        Reads SPECFEM model or kernel

        Models are stored in Fortran binary format and separated into
        multiple files according to material parameter and processor rank.
        """
        model = Model(self.parameters)
        minmax = Minmax(self.parameters)

        for iproc in range(self.mesh_properties.nproc):
            # read database files
            keys, vals = loadbypar(path, self.parameters, iproc, prefix,
                                   suffix)
            for key, val in zip(keys, vals):
                model[key] += [val]

            minmax.update(keys, vals)

        if verbose:
            minmax.write(path, logpath=PATH.SUBMIT)

        return model

    def save(self, path, model, prefix='reg1_', suffix=''):
        """
        Writes SPECFEM3D_GLOBE transerverly isotropic model

        :type path: str
        :param path:
        :type model
        :param model:
        :type prefix: str
        :param prefix: prefix that begins the name of the model parameters
        :type suffix: str
        :param suffix: that follow the name of model parameters
        """
        unix.mkdir(path)

        for iproc in range(self.mesh_properties.nproc):
            for check_key in ["vpv", "vph", "vsv", "vsh", "eta"]:
                if check_key in self.parameters:
                    savebin(model[key][iproc], path, iproc, prefix+key+suffix)
                elif 'kernel' in suffix:
                    pass
                else:
                    src = PATH.OUTPUT + '/' + 'model_init'
                    dst = path
                    copybin(src, dst, iproc, prefix+key+suffix)

            if 'rho' in self.parameters:
                savebin(model['rho'][iproc], path, iproc, prefix+'rho'+suffix)
            elif 'kernel' in suffix:
                pass
            else:
                src = PATH.OUTPUT + '/' + 'model_init'
                dst = path
                copybin(src, dst, iproc, prefix+'rho'+suffix)

    def check_mesh_properties(self, path=None, parameters=None):
        """
        Determine if Mesh properties are okay for workflow

        :type path: str
        :param path: path to the mesh file
        """
        if not hasattr(self, '_mesh_properties'):
            if path is None:
                path = PATH.MODEL_INIT

            if parameters is None:
                parameters = self.parameters

            nproc = 0
            ngll = []
            while True:
                dummy = loadbin(path, nproc, 'reg1_' + parameters[0])
                ngll += [len(dummy)]
                nproc += 1
                if not exists(
                        os.path.join(path,
                                     f"proc{nrpoc}_reg1_{parameters[0]}.bin")):
                    break

            self._mesh_properties = Struct([['nproc', nproc],
                                            ['ngll', ngll]]
                                           )

    def rename_data(self):
        """
        Works around conflicting data filename conventions

        Specfem3D's uses different name conventions for regular traces
        and 'adjoint' traces
        """
        files = glob(os.path.join(self.cwd, "traces", "adj", "*sem.ascii")
        unix.rename("sem.ascii", "sem.ascii.adj", files)

    def initialize_adjoint_traces(self):
        """
        Setup utility: Creates the "adjoint traces" expected by SPECFEM

        !!! This probably doesnt work

        Note:
            Adjoint traces are initialized by writing zeros for all channels.
            Channels actually in use during an inversion or migration will be
            overwritten with nonzero values later on.
        """
        super().initialize_adjoint_traces()

        # workaround for  SPECFEM's use of different name conventions for
        # regular traces and 'adjoint' traces
        if PAR.FORMAT.upper() in ['ASCII', 'ascii']:
            files = glob(os.path.join(self.cwd, "traces", "adj", "*sem.ascii")
            unix.rename("sem.ascii", "adj", files)

    @property
    def data_wildcard(self):
        """
        Returns a wildcard identifier for synthetic data

        :rtype: str
        :return: wildcard identifier for channels
        """
        if PAR.FORMAT.upper() == "ASCII":
            return f"*.?X?.sem.ascii"

    @property
    def data_filenames(self):
        """
        Returns the filenames of all data, either by the requested components
        or by all available files in the directory.

        :rtype: list
        :return: list of data filenames
        """
        unix.cd(os.path.join(self.cwd, "traces", "obs"))

        if PAR.FORMAT.upper() == "ASCII":
            return sorted(glob("*.???.sem.ascii"))



