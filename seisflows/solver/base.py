#!/usr/bin/env python
"""
This is the base class seisflows.solver.Base
This class provides the core utilities for the Seisflows solver interactions
"""
import os
import sys
import numpy as np
from glob import glob
from functools import partial
from seisflows.plugins import solver_io
from seisflows.tools import msg, unix
from seisflows.tools.err import ParameterError
from seisflows.tools.seismic import Container, call_solver
from seisflows.tools.tools import Struct, diff, exists


# Seisflows configuration
PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
preprocess = sys.modules['seisflows_preprocess']


class Base(object):
    """
    This base class provides an interface through which solver simulations can
    be set up and run and a parent class for the following subclasses:

    SPECFEM2D
    SPECFEM3D
    SPECFEM3D_GLOBE

    Note:
    This class supports only acoustic and isotropic elastic inversions.
    For additional options, see github.com/rmodrak/seisflows-multiparameter

    Function descriptors:

    eval_func, eval_grad, apply_hess

        These methods deal with evaluation of the misfit function or its
        derivatives.  Together, they provide the primary interface through which
        SeisFlows interacts with SPECFEM2D/3D

    forward, adjoint

        These methods allow direct access to low-level SPECFEM2D/3D components,
        providing an alternative interface through which to interact with the
        solver

    steup, generate_data, generate_model

        One-time operations performed at the beginning of inversion or migration

    initialize_solver_directories, initialize_adjoint_traces

        SPECFEM2D/3D requires a particular directory structure in which to run
        and particular file formats for models, data, and parameter files. These
        methods help put in place all these prerequisites

    load, save

        For reading and writing SPECFEM2D/3D models and kernels. On the disk,
        models and kernels are stored as binary files, and in memory, as
        dictionaries with different keys corresponding to different material
        parameters

    split, merge

        Within the solver routines, it is natural to store models as
        dictionaries. Within the optimization routines, it is natural to store
        models as vectors. Two methods, 'split' and 'merge', are used to convert
        back and forth between these two representations

    combine, smooth

        Utilities for combining and smoothing kernels

    generate_data, generate_mesh, eval_fwd, forward, adjoint,
    data_filenames, model_databases, kernel_databases, source_prefix

        !!! Required functions which must be implemented by subclass !!!

    """
    # Determine material properties to be used in workflow
    assert 'MATERIALS' in PAR
    assert 'DENSITY' in PAR

    parameters = []
    if PAR.MATERIALS == 'Elastic':
        parameters += ['vp']
        parameters += ['vs']
    elif PAR.MATERIALS == 'Acoustic':
        parameters += ['vp']

    if PAR.DENSITY == 'Variable':
        parameters += ['rho']

    def check(self):
        """
        Checks parameters and paths
        """
        # number of processors per simulation
        if "NPROC" not in PAR:
            raise ParameterError(PAR, "NPROC")

        # Format used by SPECFEM for reading, writing models
        # (currently, SPECFEM offers "fortran_binary", "adios")
        if "SOLVERIO" not in PAR:
            setattr(PAR, "SOLVERIO", "fortran_binary")

        if "VERBOSE" not in PAR:
            setattr(PAR, "VERBOSE", True)

        # Required: Solver scratch path
        if "SCRATCH" not in PATH:
            raise ParameterError(PATH, "SCRATCH")

        # To override placing solver in scratch directory
        if "LOCAL" not in PATH:
            setattr(PATH, "LOCAL", None)

        # To override the location of solver directory
        if "SOLVER" not in PATH:
            if PATH.LOCAL:
                setattr(PATH, "SOLVER", os.path.join(PATH.LOCAL, "solver"))
            else:
                setattr(PATH, "SOLVER", os.path.join(PATH.SCRATCH, "solver"))

        # Required: Solver input paths
        if "SPECFEM_BIN" not in PATH:
            raise ParameterError(PATH, "SPECFEM_BIN")

        if "SPECFEM_DATA" not in PATH:
            raise ParameterError(PATH, "SPECFEM_DATA")

        # Ensure parameters set properly
        assert self.parameters != []
        assert hasattr(solver_io, PAR.SOLVERIO)
        assert hasattr(self.io, "read_slice")
        assert hasattr(self.io, "write_slice")

    def setup(self):
        """ 
        Prepares solver for inversion or migration.
        Sets up directory structure expected by SPECFEM and copies or generates
        seismic data to be inverted or migrated

        Note:
            As input for an inversion or migration, users can choose between
            providing data, or providing a target model from which data are
            generated on the fly.
            In the former case, a value for PATH.DATA must be supplied;
            in the latter case, a value for PATH.MODEL_TRUE must be provided.
        """
        # Clean up for new inversion
        unix.rm(self.cwd)

        # If Data provided by user
        if PATH.DATA:
            # Copy user supplied data
            self.initialize_solver_directories()
            unix.cp(src=glob(os.path.join(PATH.DATA, self.source_name, '*')),
                    dst="traces/obs/")
        else:
            # Generate data on the fly using the true model
            if PAR.VERBOSE:
                print("generating data in solver/base")
            self.generate_data(model_path=PATH.MODEL_TRUE,
                               model_name='model_true',
                               model_type='gll')

        # Prepare initial model
        self.generate_mesh(model_path=PATH.MODEL_INIT,
                           model_name='model_init',
                           model_type='gll')

        self.initialize_adjoint_traces()

    def clean(self):
        """
        Clean up the run directory and set up for new run
        """
        unix.cd(self.cwd)
        unix.rm('OUTPUT_FILES')
        unix.mkdir('OUTPUT_FILES')

    def generate_data(self, *args, **kwargs):
        """
        Generates data based on a given model

        !!! Must be implemented by subclass !!!
        """
        raise NotImplementedError

    def generate_mesh(self, *args, **kwargs):
        """
        Performs meshing and database generation

        !!! Must be implemented by subclass !!!
        """
        raise NotImplementedError

    def eval_func(self, path='', export_traces=False, write_residuals=True):
        """
        High level solver interface

        Performs forward simulations and evaluates the misfit function

        :type path: str
        :param path: directory from which model is imported
        :type export_traces: bool
        :param export_traces: if True, save traces to OUTPUT.
            if False, discard traces
        :type write_residuals: bool
        :param write_residuals: calculate and export residuals
        """
        unix.cd(self.cwd)
        self.import_model(path)
        self.forward()

        if write_residuals:
            preprocess.prepare_eval_grad(self.cwd)
            self.export_residuals(path)

    def eval_fwd(self, *args, **kwargs):
        """
        High level solver interface

        Performs forward simulations needed for misfit function evaluation

        !!! Must be implemented by subclass !!!
        """
        raise NotImplementedError

    def eval_grad(self, path='', export_traces=False):
        """
        High level solver interface

        Evaluates gradient by carrying out adjoint simulations.
        Note: A function evaluation must already have been carried out.

        :type path: str
        :param path: directory from which model is imported
        :type export_traces: bool
        :param export_traces: if True, save traces to OUTPUT.
            if False, discard traces
        """
        unix.cd(self.cwd)
        self.adjoint()
        self.export_kernels(path)
        if export_traces:
            self.export_traces(path=os.path.join(path, 'traces', 'syn'),
                               prefix='traces/syn')
            self.export_traces(path=os.path.join(path, 'traces', 'adj'),
                               prefix='traces/adj')

    def apply_hess(self, path=''):
        """
        High level solver interface

        Computes action of Hessian on a given model vector.
        Note: A gradient evaluation must have already been carried out.

        :type path: str
        :param path: directory to which output files are exported
        """
        unix.cd(self.cwd)
        self.import_model(path)
        unix.mkdir('traces/lcg')
        self.forward('traces/lcg')
        preprocess.prepare_apply_hess(self.cwd)
        self.adjoint()
        self.export_kernels(path)

    def forward(self):
        """
        Low level solver interface

        Calls forward solver

        !!! Must be implemented by subclass !!!
        """
        # must be implemented by subclass
        raise NotImplementedError

    def adjoint(self):
        """
        Low level solver interface

        Calls adjoint solver

        !!! Must be implemented by subclass !!!
        """
        # must be implemented by subclass
        raise NotImplementedError

    @property
    def io(self):
        """
        Solver IO module set by User.

        Located in seisflows.plugins.solver_io
        """
        return getattr(solver_io, PAR.SOLVERIO)

    def load(self, path, parameters=[], prefix='', suffix=''):
        """ 
        Solver I/O: Loads SPECFEM2D/3D models or kernels

        :type path: str
        :param path: directory from which model is read
        :type parameters: list
        :param parameters: material parameters to be read
            (if empty, defaults to self.parameters)
        :type prefix: str
        :param prefix: optional filename prefix
        :type suffix: str
        :param suffix: optional filename suffix, eg '_kernel'
        :rtype: dict
        :return: model or kernels indexed by material parameter and
            processor rank, ie dict[parameter][iproc]
        """
        load_dict = Container()
        for iproc in range(self.mesh_properties.nproc):
            for key in parameters or self.parameters:
                load_dict[key] += self.io.read_slice(
                    path=path, parameters=f"{prefix}{key}{suffix}", iproc=iproc)

        return load_dict

    def save(self, save_dict, path, parameters=['vp','vs','rho'], prefix='',
             suffix=''):
        """ 
        Solver I/O: Saves SPECFEM2D/3D models or kernels

        :type save_dict: dict or Container
        :param save_dict: model stored as a dictionary or Container
        :type path: str
        :param path: directory from which model is read
        :type parameters: list
        :param parameters: list of material parameters to be read
        :type prefix: str
        :param prefix: optional filename prefix
        :type suffix: str
        :param suffix: optional filename suffix, eg '_kernel'
        """
        unix.mkdir(path)

        # Fill in any missing parameters
        missing_keys = diff(parameters, save_dict.keys())
        for iproc in range(self.mesh_properties.nproc):
            for key in missing_keys:
                save_dict[key] += self.io.read_slice(
                    path=PATH.MODEL_INIT, parameters=f"{prefix}{key}{suffix}",
                    iproc=iproc)

        # Write slices to disk
        for iproc in range(self.mesh_properties.nproc):
            for key in parameters:
                self.io.write_slice(data=save_dict[key][iproc], path=path,
                                    parameters=f"{prefix}{key}{suffix}",
                                    iproc=iproc)

    def merge(self, model, parameters=[]):
        """
        Convert dictionary representation `model` to vector representation `m`

        :type model: dict
        :param model: model to be converted
        :type parameters: list
        :param parameters: optional list of parameters,
            defaults to `self.parameters`
        :rtype: np.ndarray
        :return: model as a vector
        """
        m = np.array([])
        for key in parameters or self.parameters:
            for iproc in range(self.mesh_properties.nproc):
                m = np.append(m, model[key][iproc])

        return m

    def split(self, m, parameters=[]):
        """
        Converts vector representation `m` to dictionary representation `model`

        :type m: np.ndarray
        :param m: model to be converted
        :type parameters: list
        :param parameters: optional list of parameters,
            defaults to `self.parameters`
        :rtype: dict
        :return: model as a dictionary
        """
        nproc = self.mesh_properties.nproc
        ngll = self.mesh_properties.ngll
        model = Container()

        for idim, key in enumerate(parameters or self.parameters):
            model[key] = []
            for iproc in range(nproc):
                imin = sum(ngll) * idim + sum(ngll[:iproc])
                imax = sum(ngll) * idim + sum(ngll[:iproc + 1])
                model[key] += [m[imin:imax]]

        return model

    def combine(self, input_path, output_path, parameters=[]):
        """
        Postprocessing wrapper: xcombine_sem
        Sums individual source contributions.

        :type input_path: str
        :param input_path: path to data
        :type output_path: str
        :param output_path: path to export the outputs of xcombine_sem
        :type parameters: list
        :param parameters: optional list of parameters,
            defaults to `self.parameters`
        """
        if not exists(output_path):
            unix.mkdir(output_path)

        unix.cd(self.cwd)
        # Write the source names into the kernel paths file for SEM
        with open('kernel_paths', 'w') as file:
            file.writelines(
                [os.path.join(input_path, f"{name}\n")
                 for name in self.source_names]
            )

        # Call on bin/xcombine_sem
        for name in parameters or self.parameters:
            # Example call:
            # mpiexec ./bin/xcombine_sem alpha_kernel kernel_paths output
            call_solver(mpiexec=system.mpiexec(),
                        executable=" ".join([f"{PATH.SPECFEM_BIN}/xcombine_sem",
                                             f"{name}_kernel", "kernel_paths",
                                             output_path])
                        )

    def smooth(self, input_path, output_path, parameters=[], span_h=0.,
               span_v=0., output='solver.log'):
        """
        Postprocessing wrapper: xsmooth_sem
        Smooths kernels by convolving them with a Gaussian.

        Note:
            paths require a trailing `/` character when calling xsmooth_sem

        :type input_path: str
        :param input_path: path to data
        :type output_path: str
        :param output_path: path to export the outputs of xcombine_sem
        :type parameters: list
        :param parameters: optional list of parameters,
            defaults to `self.parameters`
        :type span_h: float
        :param span_h: horizontal smoothing length in meters
        :type span_v: float
        :param span_v: vertical smoothing length in meters
        :type output: str
        :param output: file to output stdout to
        """
        if not exists(output_path):
            unix.mkdir(output_path)

        # Apply smoothing operator
        unix.cd(self.cwd)
        for name in parameters or self.parameters:
            if PAR.VERBOSE:
                print(f"Smoothing {name}")
            call_solver(mpiexec=system.mpiexec(),
                        executable=" ".join([f"{PATH.SPECFEM_BIN}/xsmooth_sem",
                                             str(span_h), str(span_v),
                                             f"{name}_kernel",
                                             os.path.join(input_path, ""),
                                             os.path.join(output_path, ""),
                                             ".false"]),
                        output=output)

        # Rename output files
        files = glob(os.path.join(output_path, '*'))
        unix.rename(old='_smooth', new='', names=files)

    def combine_vol_data_vtk(self):
        """
        Postprocessing wrapper for xcombine_vol_data_vtk

        !!! must be implemented by subclass !!!
        """
        pass

    def import_model(self, path):
        """
        File transfer utility. Import the model into the workflow.

        :type path: str
        :param path: path to model
        """
        model = self.load(os.path.join(path, "model"))
        self.save(model, self.model_databases)

    def import_traces(self, path):
        """
        File transfer utility. Import traces into the workflow.

        :type path: str
        :param path: path to traces
        """
        src = glob(os.path.join(path, 'traces', self.source_name, '*'))
        dst = os.path.join(self.cwd, 'traces', 'obs')
        unix.cp(src, dst)

    def export_model(self, path, parameters=['rho', 'vp', 'vs']):
        """
        File transfer utility. Export the model to disk.

        Performed by master solver.

        :type path: str
        :param path: path to save model
        :type parameters: list
        :param parameters: list of parameters that define the model
        """
        if self.taskid == 0:
            unix.mkdir(path)
            for key in parameters:
                files = glob(os.path.join(self.model_databases,
                                          f"*{key}.bin")
                             )
                unix.cp(files, path)

    def export_kernels(self, path):
        """
        File transfer utility. Export kernels to disk

        :type path: str
        :param path: path to save kernels
        """
        unix.cd(self.kernel_databases)

        # Work around conflicting name conventions
        self.rename_kernels()

        src = glob('*_kernel.bin')
        dst = os.path.join(path, 'kernels', self.source_name)
        unix.mkdir(dst)
        unix.mv(src, dst)

    def export_residuals(self, path):
        """
        File transfer utility. Export residuals to disk.

        :type path: str
        :param path: path to save residuals
        """
        unix.mkdir(os.path.join(path, 'residuals'))

        src = os.path.join(self.cwd, 'residuals')
        dst = os.path.join(path, 'residuals', self.source_name)
        unix.mv(src, dst)

    def export_traces(self, path, prefix='traces/obs'):
        """
        File transfer utility. Export traces to disk.

        :type path: str
        :param path: path to save traces
        :type prefix: str
        :param prefix: location of traces w.r.t self.cwd
        """
        unix.mkdir(os.path.join(path))

        src = os.path.join(self.cwd, prefix)
        dst = os.path.join(path, self.source_name)
        unix.cp(src, dst)

    def rename_kernels(self):
        """
        Works around conflicting kernel filename conventions by renaming
        `alpha` to `vp` and `beta` to `vs`
        """
        # Rename alpha to vp
        for globfids in ['*proc??????_alpha_kernel.bin',
                         '*proc??????_alpha[hv]_kernel.bin',
                         '*proc??????_reg1_alpha_kernel.bin',
                         '*proc??????_reg1_alpha[hv]_kernel.bin']:
            unix.rename(old='alpha', new='vp', names=glob(globfids))

        # Rename beta to vs
        for globfids in ['*proc??????_beta_kernel.bin',
                         '*proc??????_beta[hv]_kernel.bin',
                         '*proc??????_reg1_beta_kernel.bin',
                         '*proc??????_reg1_beta[hv]_kernel.bin']:
            unix.rename(old='beta', new='vs', names=glob(globfids))

    def rename_data(self, path):
        """
        Works around conflicting data filename conventions
        """
        pass

    # setup utilities
    def initialize_solver_directories(self):
        """
        Setup Utility:
        Creates directory structure expected by SPECFEM3D, copies
        executables, and prepares input files. Executables must be supplied
        by user as there is currently no mechanism for automatically
        compiling from source.

        Directories will act as completely indepedent Specfem run directories.
        This allows for embarassing parallelization while avoiding the need 
        for intra-directory commnunications, at the cost of redundancy and 
        extra files
        """
        unix.mkdir(self.cwd)
        unix.cd(self.cwd)

        # Create directory structure
        unix.mkdir("bin")
        unix.mkdir("DATA")
        unix.mkdir("OUTPUT_FILES")
        unix.mkdir("traces/obs")
        unix.mkdir("traces/syn")
        unix.mkdir("traces/adj")
        unix.mkdir(self.model_databases)
        unix.mkdir(self.kernel_databases)

        # Copy exectuables
        src = glob(os.path.join(PATH.SPECFEM_BIN, "*"))
        dst = os.path.join("bin", "")
        unix.cp(src, dst)

        # Copy all input files except source files
        src = glob(os.path.join(PATH.SPECFEM_DATA, "*"))
        src = [_ for _ in src if not _.startswith(self.source_prefix)]
        dst = os.path.join("DATA", "")
        unix.cp(src, dst)

        # symlink event source specifically
        src = os.path.join(PATH.SPECFEM_DATA, 
                           f"{self.source_prefix}_{self.source_name}")
        dst = os.path.join("DATA", self.source_prefix)
        unix.ln(src, dst)

        if self.taskid == 0: 
            # Symlink taskid_0 as mainsolver in solver directory for convenience
            unix.ln(self.source_name, os.path.join(PATH.SOLVER, "mainsolver"))
        else:
            # Copy the initial model from mainsolver into current directory
            # Avoids the need to run multiple instances of xgenerate_databases
            src = os.path.join(PATH.SOLVER, "mainsolver", "OUTPUT_FILES", 
                               "DATABASES_MPI")
            dst = os.path.join(self.cwd, "OUTPUT_FILES", "DATABASES_MPI")
            unix.cp(src, dst)

        self.check_solver_parameter_files()

    def initialize_adjoint_traces(self):
        """
        Setup utility: Creates the "adjoint traces" expected by SPECFEM

        Note:
            Adjoint traces are initialized by writing zeros for all channels.
            Channels actually in use during an inversion or migration will be
            overwritten with nonzero values later on.
        """
        for filename in self.data_filenames:
            # Read traces
            d = preprocess.reader(path=os.path.join(self.cwd, 'traces', 'obs'),
                                  filename=filename)

            # Initialize by writing zeroes
            for t in d:
                t.data[:] = 0.

            # Write traces
            preprocess.writer(stream=d,
                              path=os.path.join(self.cwd, 'traces', 'adj'),
                              filename=filename)

    def check_mesh_properties(self, path=None):
        """
        Determine if Mesh properties are okay for workflow

        :type path: str
        :param path: path to the mesh file
        """
        # Check the given model path or the initial model
        path = path or PATH.MODEL_INIT
        if not exists(path):
            raise FileNotFoundError(f"Mesh path {path} does not exist")

        # Count slices and grid points
        key = self.parameters[0]
        iproc = 0
        ngll = []
        while True:
            dummy = self.io.read_slice(
                path=path, parameters=key, iproc=iproc)[0]
            ngll += [len(dummy)]
            iproc += 1
            if not exists(f"{path}/proc{iproc:06d}_{key}.bin"):
                break
        nproc = iproc

        # Create coordinate pointers
        coords = Struct()
        for key in ['x', 'y', 'z']:
            coords[key] = partial(self.io.read_slice, self, path, key)

        # Define internal mesh properties
        self._mesh_properties = Struct([['nproc', nproc],
                                        ['ngll', ngll],
                                        ['path', path],
                                        ['coords', coords]]
                                       )

    def check_source_names(self):
        """
        Determines names of sources by applying wildcard rule to
        user-supplied input files

        Note:
            Available sources are sorted and chosen from the start of the list
            until PAR.NTASK
        """
        # Check path
        path = PATH.SPECFEM_DATA
        if not exists(path):
            raise Exception

        # Apply wildcard rule and check for available sources
        wildcard = f"{self.source_prefix}_*"
        globstar = sorted(glob(os.path.join(path, wildcard)))
        if not globstar:
            print(msg.SourceError_SPECFEM.format(path, wildcard))
            sys.exit(-1)

        # Create internal definition of available source names
        names = []
        for path in globstar:
            names += [os.path.basename(path).split('_')[-1]]

        self._source_names = names[:PAR.NTASK]

    def check_solver_parameter_files(self):
        """
        Optional method

        !!! Can be implemented by subclass !!!
        """
        pass

    @property
    def taskid(self):
        """
        It is sometimes useful to overload system.taskid

        :rtype: int
        :return: task id for given solver
        """
        return system.taskid()

    @property
    def source_name(self):
        """
        Returns name of source currently under consideration

        :rtype: str
        :return: given source name for given task id
        """
        return self.source_names[self.taskid]

    @property
    def cwd(self):
        """
        Returns working directory currently in use

        :rtype: str
        :return: current solver working directory
        """
        return os.path.join(PATH.SOLVER, self.source_name)

    @property
    def source_names(self):
        """
        Returns list of source names

        :rtype: list
        :return: list of source names
        """
        if not hasattr(self, '_source_names'):
            self.check_source_names()

        return self._source_names

    @property
    def mesh_properties(self):
        """
        Returns mesh properties

        :rtype: Struct
        :return: Structure of mesh properties
        """
        if not hasattr(self, '_mesh_properties'):
            self.check_mesh_properties()

        return self._mesh_properties

    @property
    def data_filenames(self):
        """
        Template filenames for accessing data

        !!! Must be implemented by subclass !!!
        """
        return NotImplementedError

    @property
    def model_databases(self):
        """
        Template filenames for accessing models

        !!! Must be implemented by subclass !!!
        """
        return NotImplementedError

    @property
    def kernel_databases(self):
        """
        Template filenames for accessing kernels

        !!! Must be implemented by subclass !!!
        """
        return NotImplementedError

    @property
    def source_prefix(self):
        """
        Template filenames for accessing sources

        !!! Must be implemented by subclass !!!
        """
        return NotImplementedError


