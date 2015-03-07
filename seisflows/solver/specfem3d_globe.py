
import subprocess
from glob import glob
from os.path import join

import numpy as np

import seisflows.seistools.specfem3d_globe as solvertools
from seisflows.seistools.io import savebin
from seisflows.seistools.shared import getpar, setpar
from seisflows.seistools.io import loadbypar, loadbyproc, savebin, \
    applymap, ModelStruct, MinmaxStruct

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists
from seisflows.tools.config import loadclass, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class specfem3d_globe(loadclass('solver', 'base')):
    """ Python interface for SPECFEM3D_GLOBE

      See base class for method descriptions
    """

    # By default, use transversely isotropic model updates.
    parameters = []
    parameters += ['rho']
    parameters += ['vpv']
    parameters += ['vph']
    parameters += ['vsv']
    parameters += ['vsh']
    parameters += ['eta']


    def check(self):
        """ Checks parameters and paths
        """
        super(specfem3d_globe, self).check()


    def generate_data(self, **model_kwargs):
        """ Generates data
        """
        self.generate_mesh(**model_kwargs)

        unix.cd(self.getpath)
        setpar('SIMULATION_TYPE', '1')
        setpar('SAVE_FORWARD', '.true.')
        self.mpirun('bin/xspecfem3D')

        unix.mv(self.data_wildcard, 'traces/obs')
        self.export_traces(PATH.OUTPUT, 'traces/obs')


    def generate_mesh(self, model_path=None, model_name=None, model_type='gll'):
        """ Performs meshing and database generation
        """
        assert(model_name)
        assert(model_type)

        self.initialize_solver_directories()
        unix.cd(self.getpath)

        if model_type == 'gll':
            assert (exists(model_path))
            unix.cp(glob(model_path +'/'+ '*'), self.model_databases)
            self.mpirun('bin/xmeshfem3D')
            self.export_model(PATH.OUTPUT +'/'+ model_name)

        else:
            raise NotImplementedError


    ### model input/output

    def load(self, path, parameters=None, mapping=None, prefix='reg1_', suffix='', verbose=False):
        """ reads SPECFEM model

          Models are stored in Fortran binary format and separated into multiple
          files according to material parameter and processor rank.

          Optionally, 'mapping' can be used to convert on the fly from one set
          of material parameters to another.
        """
        if not parameters:
            parameters = self.parameters

        model = ModelStruct(parameters, mapping)
        minmax = MinmaxStruct(parameters, mapping)

        for iproc in range(PAR.NPROC):
            keys, vals = loadbypar(path, parameters, iproc, prefix, suffix)

            # keep track of min, max
            minmax.update(keys, vals)

            # apply optional mapping
            if mapping:
                keys, vals = applymap(vals, mapping)

            # append latest values
            for key, val in zip(keys, vals):
                model[key] += [val]

        if verbose:
            minmax.write(path, logpath=PATH.SUBMIT)

        return model


    def save(self, path, model, prefix='reg1_', suffix=''):
        """ writes SPECFEM3D_GLOBE transerverly isotropic model
        """
        unix.mkdir(path)

        for iproc in range(PAR.NPROC):
            for key in ['vpv', 'vph', 'vsv', 'vsh', 'eta']:
                if key in self.parameters:
                    savebin(model[key][iproc], path, iproc, prefix + key + suffix)
                else:
                    self.copybin(path, iproc, key)

            if 'rho' in self.parameters:
                savebin(model['rho'][iproc], path, iproc, prefix + 'rho' + suffix)
            elif self.density_scaling:
                raise NotImplementedError
            else:
                self.copybin(path, iproc, prefix + 'rho')


    ### low-level solver interface

    def forward(self):
        """ Calls SPECFEM3D_GLOBE forward solver
        """
        solvertools.setpar('SIMULATION_TYPE', '1')
        solvertools.setpar('SAVE_FORWARD', '.true.')
        self.mpirun('bin/xspecfem3D')
        unix.mv(self.data_wildcard, 'traces/syn')


    def adjoint(self):
        """ Calls SPECFEM3D_GLOBE adjoint solver
        """
        solvertools.setpar('SIMULATION_TYPE', '3')
        solvertools.setpar('SAVE_FORWARD', '.false.')
        unix.rm('SEM')
        unix.ln('traces/adj', 'SEM')
        self.mpirun('bin/xspecfem3D')


    ### utility functions

    @property
    def data_wildcard(self):
        return glob('OUTPUT_FILES/*.sem.ascii')

    @property
    def model_databases(self):
        return join(self.getpath, 'OUTPUT_FILES/DATABASES_MPI')

    @property
    def source_prefix(self):
        return 'CMTSOLUTION'


