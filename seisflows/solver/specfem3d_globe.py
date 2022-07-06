#!/usr/bin/env python3
"""
This is the subclass seisflows.solver.specfem3d_globe
This class provides utilities for the Seisflows solver interactions with
Specfem3D Globe. It inherits all attributes from seisflows.solver.specfem3d,
and overwrites these functions to provide specified interaction with Specfem3D.

SPECFEM3D_Globe specfic notes:
    - does not allow SU seismogram outputs, only ASCII, SAC, ASDF, 3D_Array
"""
import os
from glob import glob

from seisflows.solver.specfem3d import Specfem3D
from seisflows.tools import unix


class Specfem3DGlobe(Specfem3D):
    """
    Python interface to Specfem3D Globe. A very simple overload of Specfem3D

    See class `seisflows.solver.specfem3d.Specfem3D` for a more detailed
    explanation of methods and attributes of this class
    """
    def __init__(self):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.
        """
        super().__init__()

    def data_wildcard(self, comp="?"):
        """
        Returns a wildcard identifier for synthetic data

        :rtype: str
        :return: wildcard identifier for channels
        """
        if self.par.FORMAT.upper() == "SU":
            raise NotImplementedError("SU file access is still a WIP")
        elif self.par.FORMAT.upper() == "ASCII":
            return f"*.?X{comp}.sem.ascii"

    def load(self, path, prefix="reg1_", suffix="",  parameters=None):
        """
        Reads SPECFEM model or kernel

        .. note::
            SPECFEM3D_Globe meshes are broken into 3 regions. 
            Region 1 == Crust + Mantle
            Region 2 == Outer core
            Region 3 == Inner core

        .. warning::
            Currently SeisFlows only considers the crust + mantle in Globe
            simulations
        

        :type path: str
        :param path: directory from which model is read
        :type prefix: str
        :param prefix: optional filename prefix
        :type suffix: str
        :param suffix: optional filename suffix, eg '_kernel'
        :type parameters: list
        :param parameters: material parameters to be read
            (if empty, defaults to self.parameters)
        :rtype: dict
        :return: model or kernels indexed by material parameter and
            processor rank, ie dict[parameter][iproc]
        """
        parameters = parameters or self.parameters

        model = Model(parameters)
        minmax = Minmax(parameters)

        for iproc in range(self.mesh_properties.nproc):
            # read database files based on parameters
            keys, vals = loadbypar(path, self.parameters, iproc, prefix,
                                   suffix)
            for key, val in zip(keys, vals):
                model[key] += [val]

            minmax.update(keys, vals)

        return model

    def save(self, path, model, prefix="reg1_", suffix=""):
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
                    src = self.path.OUTPUT + '/' + 'model_init'
                    dst = path
                    copybin(src, dst, iproc, prefix+key+suffix)

            if 'rho' in self.parameters:
                savebin(model['rho'][iproc], path, iproc, prefix+'rho'+suffix)
            elif 'kernel' in suffix:
                pass
            else:
                src = self.path.OUTPUT + '/' + 'model_init'
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
                path = self.path.MODEL_INIT

            if parameters is None:
                parameters = self.parameters

            nproc = 0
            ngll = []
            while True:
                dummy = loadbin(path, nproc, 'reg1_' + parameters[0])
                ngll += [len(dummy)]
                nproc += 1
                if not os.path.exists(
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
        files = glob(os.path.join(self.cwd, "traces", "adj", "*sem.ascii"))
        unix.rename("sem.ascii", "sem.ascii.adj", files)

    def initialize_adjoint_traces(self):
        """
        Setup utility: Creates the "adjoint traces" expected by SPECFEM

        !!! This probably doesnt work

        .. note::
            Adjoint traces are initialized by writing zeros for all channels.
            Channels actually in use during an inversion or migration will be
            overwritten with nonzero values later on.
        """
        super().initialize_adjoint_traces()

        # workaround for  SPECFEM's use of different name conventions for
        # regular traces and 'adjoint' traces
        if self.par.FORMAT.upper() in ['ASCII', 'ascii']:
            files = glob(os.path.join(self.cwd, "traces", "adj", "*sem.ascii"))
            unix.rename("sem.ascii", "adj", files)




