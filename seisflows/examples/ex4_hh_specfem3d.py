#!/usr/bin/env python3
"""
                SEISFLOWS SPECFEM3D HOMOGENEOUS HALFSPACE EXAMPLE 4

This example will run an inversion to assess misfit between a homogeneous
halfspace model and the same model with a slight velocity perturbation, using
SPECFEM3D_Cartesian, a grid of surface stations and one or more events. The
default problem is small enough to run on a single compute node (28 cores).

.. note::
    Nothing is downloaded or compiled. Point the example at a SPECFEM3D
    repository that has already been compiled with `--specfem3d_repo`.

.. note::
    The tasks involved include:
    1. Check the SPECFEM3D binaries, set up a SPECFEM3D working directory
    2. Generate starting model from the SPECFEM3D homogeneous halfspace example
    3. Generate target model by perturbing the starting model
    4. Set up a SeisFlows working directory
    5. Run an inversion workflow

.. warning::
    Generating the starting and target models runs the SPECFEM3D mesher and
    database generator with MPI on `--nproc` cores. On a cluster, run this from
    an interactive session on a compute node rather than on the login node.

.. rubric::
    $ seisflows examples run 4 --specfem3d_repo /path/to/specfem3d
"""
import os
import sys
import glob
import difflib
import subprocess
import numpy as np

from seisflows.tools import msg
from seisflows.tools.config import Dict, load_yaml
from seisflows.seisflows import SeisFlows
from seisflows.tools.specfem import getpar, setpar
from seisflows.tools.unix import cd, cp, rm, ln, mv, mkdir
from seisflows.tools.unix import nproc as system_nproc


class SFExample3D:
    """
    A class for running the SPECFEM3D SeisFlows example. Follows the structure
    of the SPECFEM2D examples (see `ex1_homogeneous_halfspace`) so that the
    two read alike.
    """
    # SPECFEM3D executables required by SeisFlows to run this example
    _binaries = ["xmeshfem3D", "xgenerate_databases", "xspecfem3D",
                 "xcombine_sem", "xsmooth_sem_pde", "xcombine_vol_data_vtk"]

    def __init__(self, ntask=None, nsta=None, niter=None, nproc=None,
                 method="run", specfem3d_repo=None, specfem2d_repo=None,
                 mpiexec="mpirun", perturb=None, system=None, par=None,
                 par_yaml=None, **kwargs):
        """
        Set path structure which is used to navigate around the SPECFEM3D
        repository and the example working directory

        :type ntask: int
        :param ntask: number of events to use in inversion. The first is at the
            center of the domain, any others are spaced evenly on a ring
            around it. Defaults to 1
        :type nsta: int
        :param nsta: number of stations, which are placed on a grid of
            equidistant points on the surface of the domain. Defaults to 25
        :type niter: int
        :param niter: number of iterations to run. Defaults to 1
        :type nproc: int
        :param nproc: number of cores used by each simulation, which also sets
            how the mesh is divided between them. Defaults to 28, the size of
            a node on Chinook
        :type method: str
        :param method: method for running the example problem, can be:
            * 'run': setup and run the example problem
            * 'setup': only setup the example problem, do not `submit` job
        :type specfem3d_repo: str
        :param specfem3d_repo: path to the SPECFEM3D repository, which must
            contain binary executables in bin/. If not given, assumes it is
            called 'specfem3d/' in the current working directory
        :type specfem2d_repo: str
        :param specfem2d_repo: allows the repository to be given with the
            `-r` flag of `seisflows examples`, which names it for SPECFEM2D
        :type mpiexec: str
        :param mpiexec: MPI executable used to run the SPECFEM3D mesher and
            database generator that make the models. Defaults to 'mpirun'
        :type perturb: float
        :param perturb: percentage by which the target model's velocities
            differ from the starting model, e.g., -5 for 5% slower. Defaults
            to 5
        :type system: str
        :param system: name of the SeisFlows system module to run with, e.g.,
            'chinook' or 'miyabi'. Defaults to 'workstation'
        :type par: list of str
        :param par: SeisFlows parameters to set once the parameter file is made
            as 'KEY=VALUE' strings, e.g., ['partition=t1small', 'walltime=30']
        :type par_yaml: str
        :param par_yaml: path to a YAML file of SeisFlows parameters to set,
            for settings that are re-used (e.g., one file per cluster).
            Applied before `par`
        """
        self.cwd = os.getcwd()
        self.sem3d_paths, self.workdir_paths = self.define_dir_structures(
            cwd=self.cwd, specfem3d_repo=specfem3d_repo or specfem2d_repo
        )

        # We set defaults here because `seisflows examples` may input these
        # values as NoneType which would override __init__ defaults.
        self.ntask = ntask or 1
        assert(1 <= self.ntask), \
            f"number of tasks/events must be at least 1, not {self.ntask}"
        self.nsta = nsta or 25
        assert(1 <= self.nsta), \
            f"number of stations must be at least 1, not {self.nsta}"
        self.niter = niter or 1
        assert(1 <= self.niter), \
            f"number of iterations must be at least 1, not {self.niter}"
        self.nproc = nproc or 28
        assert(self.nproc <= system_nproc()), (
            f"your system has a maximum {system_nproc()} processors, which is "
            f"less than the requested value of {self.nproc}. Run from a compute"
            f" node, or adjust with `--nproc`"
        )
        self.perturb = 5. if perturb is None else float(perturb)
        assert(self.perturb != 0.), "`perturb` must not be 0"

        self.mpiexec = mpiexec
        self.system = system
        self.par = list(par or [])
        self.par_yaml = par_yaml

        # This bool information is provided by the User running 'setup' or 'run'
        self.run_example = bool(method == "run")

        # The mesh is divided between processors, and keeps elements about the
        # size of the original example's (3.7km)
        self.nproc_xi, self.nproc_eta = self.divide_mesh(self.nproc)
        self.nex_xi = self.nproc_xi * max(1, int(36 / self.nproc_xi + .5))
        self.nex_eta = self.nproc_eta * max(1, int(36 / self.nproc_eta + .5))

        # Command line tool to use $ seisflows <cmd> from inside Python
        # Zero out sys.argv to ensure that no arguments are given to the CLI
        sys.argv = [sys.argv[0]]
        self.sf = SeisFlows()

        # Set the main SeisFlows modules prior to running the configure() cmd.
        self._modules = {
            "workflow": "inversion",
            "solver": "specfem3d",
        }
        if self.system:
            self._modules["system"] = self.system

        # SeisFlows parameters are set as an attribute so that other examples
        # can overwrite or override
        self._parameters = {
            "ntask": self.ntask,
            "materials": "elastic",  # how velocity model parameterized
            "update_density": False,  # update density or keep constant
            "generate_data": True,  # run simulations through model_tru
            "syn_data_format": "ascii",  # how to output synthetic seismograms
            "obs_data_format": "ascii",
            "unit_output": "disp",
            "nproc": self.nproc,  # number of cores to use for MPI tasks
            "tasktime": 30,  # minutes allowed for each simulation
            "start": 1,  # first iteration
            "end": self.niter,  # final iteration
            "step_count_max": 5,
            "step_len_init": 0.05,  # as a percentage of initial model
            "components": "XYZ",  # SPECFEM3D's Cartesian components
            "attenuation": False,
            "plot_waveforms": False,
            "misfit": "traveltime",  # cross-correlation phase measure
            "adjoint": "traveltime",  # cross-correlation phase measure
            "smooth_type": "pde",  # the default, gaussian, is very slow
            "smooth_h": 10000.,  # smoothing lengths in meters
            "smooth_v": 10000.,
            "path_specfem_bin": self.workdir_paths.bin,
            "path_specfem_data": self.workdir_paths.data,
            "path_model_init": self.workdir_paths.model_init,
            "path_model_true": self.workdir_paths.model_true,
        }
        # A cluster's system module has its own MPI executable
        if not self.system:
            self._parameters["mpiexec"] = self.mpiexec

    def print_dialogue(self):
        """
        Print help/system dialogue message that explains the setup of th
        this workflow
        """
        print(msg.ascii_logo_small)
        print(msg.cli(
            f"This is a [SPECFEM3D] [{(self.system or 'workstation').upper()}]"
            f" example, which will run an inversion to assess misfit between "
            f"two homogeneous halfspace models with a {self.perturb:+g}% "
            f"velocity difference. [{self.ntask} events, {self.nsta} "
            f"stations, {self.niter} iterations, {self.nproc} cores "
            f"({self.nproc_xi}x{self.nproc_eta} mesh slices)]. "
            f"The tasks involved include: ",
            items=["1. Check the compiled SPECFEM3D binaries",
                   "2. [Setup] a SPECFEM3D working directory",
                   "3. [Setup] starting model from the SPECFEM3D example",
                   "4. [Setup] target model w/ perturbed starting model",
                   "5. [Setup] a SeisFlows working directory",
                   "6. [Run] the inversion workflow"],
            header="seisflows example 4",
            border="=")
        )

    @staticmethod
    def divide_mesh(nproc):
        """
        Divide the mesh between `nproc` processors as evenly as possible along
        its two horizontal directions, e.g., 28 -> 7 x 4, 36 -> 6 x 6

        :type nproc: int
        :param nproc: number of processors
        :rtype: (int, int)
        :return: number of slices along xi, along eta
        """
        nproc_eta = int(np.sqrt(nproc))
        while nproc % nproc_eta:
            nproc_eta -= 1
        return nproc // nproc_eta, nproc_eta

    @staticmethod
    def define_dir_structures(cwd, specfem3d_repo):
        """
        Define the example directory structure, which will contain abridged
        versions of the SPECFEM3D working directory.

        :type cwd: str
        :param cwd: current working directory
        :type specfem3d_repo: str
        :param specfem3d_repo: location of the SPECFEM3D repository
        """
        if specfem3d_repo is None:
            specfem3d_repo = os.path.join(cwd, "specfem3d")
        specfem3d_repo = os.path.abspath(specfem3d_repo)

        sem3d = {
            "repo": specfem3d_repo,
            "bin": os.path.join(specfem3d_repo, "bin"),
        }

        # The example is called differently depending on the version
        examples_dir = os.path.join(specfem3d_repo, "EXAMPLES", "applications")
        hex8 = "homogeneous_halfspace_HEX8_elastic_absorbing_Stacey_5sides"
        for name in ["homogeneous_halfspace", hex8]:
            sem3d["example"] = os.path.join(examples_dir, name)
            if os.path.exists(sem3d["example"]):
                break

        # This defines a working directory structure which we will create
        working_directory = os.path.join(cwd, "specfem3d_workdir")
        workdir = {
            "workdir": working_directory,
            "bin": os.path.join(working_directory, "bin"),
            "data": os.path.join(working_directory, "DATA"),
            "output": os.path.join(working_directory, "OUTPUT_FILES"),
            "model_init": os.path.join(working_directory, "OUTPUT_FILES_INIT"),
            "model_true": os.path.join(working_directory, "OUTPUT_FILES_TRUE"),
        }
        return Dict(sem3d), Dict(workdir)

    def check_specfem3d(self):
        """
        Check that the SPECFEM3D repository has the example and the compiled
        binaries that SeisFlows needs. Nothing is compiled by this example
        """
        if not os.path.exists(self.sem3d_paths.example):
            print(msg.cli(
                f"Cannot find the SPECFEM3D homogeneous halfspace example in "
                f"'{self.sem3d_paths.repo}'. Point to a SPECFEM3D repository "
                f"with `--specfem3d_repo`",
                header="missing specfem3d", border="="))
            sys.exit(-1)

        missing = [exc for exc in self._binaries if not
                   os.path.exists(os.path.join(self.sem3d_paths.bin, exc))]
        if missing:
            print(msg.cli(
                f"SPECFEM3D binaries are missing from "
                f"'{self.sem3d_paths.bin}'. Compile SPECFEM3D (e.g., "
                f"`./configure --with-mpi; make all` in the repository) and "
                f"try again. Missing:", items=missing,
                header="missing specfem3d binaries", border="="))
            sys.exit(-1)

        try:
            subprocess.run(f"which {self.mpiexec}", check=True, shell=True,
                           stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            print(msg.cli(
                f"MPI executable '{self.mpiexec}' not found on system. Please "
                f"check that you have MPI loaded. If '{self.mpiexec}' is not "
                f"how you invoke MPI, use the flag `--mpiexec` to change it",
                header="missing mpi executable", border="="))
            sys.exit(-1)

    def create_specfem3d_working_directory(self):
        """
        Create the working directory where we will generate our initial and
        final models using the SPECFEM3D homogeneous halfspace example. Sets up
        the mesh, the processors and the length of the simulation.
        """
        example = self.sem3d_paths.example

        # Incase this example has written before, remove dir. that were created
        rm(self.workdir_paths.workdir)
        mkdir(self.workdir_paths.workdir)
        mkdir(self.workdir_paths.bin)

        # Copy the binary executables, and DATA from the SPECFEM3D example
        cp([os.path.join(self.sem3d_paths.bin, exc) for exc in self._binaries],
           self.workdir_paths.bin)
        mkdir(self.workdir_paths.data)
        cp(os.path.join(example, "DATA", "Par_file"), self.workdir_paths.data)
        cp(os.path.join(example, "meshfem3D_files"),
           os.path.join(self.workdir_paths.data, "meshfem3D_files"))

        print("> Setting the SPECFEM3D Par_file and mesh for this example")
        par_file = os.path.join(self.workdir_paths.data, "Par_file")
        setpar("NPROC", self.nproc, par_file)
        setpar("NSTEP", 2000, par_file)  # 100s

        mesh_par_file = os.path.join(self.workdir_paths.data, "meshfem3D_files",
                                     "Mesh_Par_file")
        setpar("NEX_XI", self.nex_xi, mesh_par_file)
        setpar("NEX_ETA", self.nex_eta, mesh_par_file)
        setpar("NPROC_XI", self.nproc_xi, mesh_par_file)
        setpar("NPROC_ETA", self.nproc_eta, mesh_par_file)
        # The region is one block of material from the top to bottom of the mesh
        self._set_mesh_table_line(mesh_par_file, "NREGIONS",
                                  ["1", str(self.nex_xi), "1",
                                   str(self.nex_eta), "1", "16", "1"])

        # Remember the example's material, the target model is scaled from it
        with open(mesh_par_file, "r") as f:
            lines = f.readlines()
        self._material = lines[self._mesh_table_line(lines, "NMATERIALS")
                               ].split()

        self.write_sources()
        self.write_stations()

    @staticmethod
    def _mesh_table_line(lines, key):
        """
        Mesh_Par_file tables have their values on the first line that is not a
        comment after the parameter, e.g., the material and region definitions

        :type lines: list of str
        :param lines: lines of the Mesh_Par_file
        :type key: str
        :param key: parameter after which the table line follows
        :rtype: int
        :return: index of the line that holds the table values
        """
        start = [i for i, line in enumerate(lines) if
                 line.split("=")[0].strip() == key][0]
        for i in range(start + 1, len(lines)):
            if lines[i].strip() and not lines[i].strip().startswith("#"):
                return i

    def _set_mesh_table_line(self, mesh_par_file, key, values):
        """
        Set the values of a table in the Mesh_Par_file, see `_mesh_table_line`

        :type mesh_par_file: str
        :param mesh_par_file: path to the Mesh_Par_file
        :type key: str
        :param key: parameter after which the table line follows
        :type values: list of str
        :param values: the new contents of the line
        """
        with open(mesh_par_file, "r") as f:
            lines = f.readlines()
        lines[self._mesh_table_line(lines, key)] = "  ".join(values) + "\n"
        with open(mesh_par_file, "w") as f:
            f.writelines(lines)

    def _mesh_extent(self):
        """
        The horizontal extent of the mesh in meters

        :rtype: (float, float, float, float)
        :return: minimum and maximum along latitude, then along longitude
        """
        mesh_par_file = os.path.join(self.workdir_paths.data, "meshfem3D_files",
                                     "Mesh_Par_file")
        return tuple(float(getpar(key, mesh_par_file, delim="=")[1])
                     for key in ["LATITUDE_MIN", "LATITUDE_MAX",
                                 "LONGITUDE_MIN", "LONGITUDE_MAX"])

    def write_sources(self):
        """
        Write a CMTSOLUTION_??? file for each event, based on the example's
        source. The first is at the center of the domain and the rest are on a
        ring around it, all with the same mechanism and depth.
        """
        lat_min, lat_max, lon_min, lon_max = self._mesh_extent()
        lat_c, lon_c = (lat_min + lat_max) / 2, (lon_min + lon_max) / 2
        radius = 0.2 * min(lat_max - lat_min, lon_max - lon_min)

        with open(os.path.join(self.sem3d_paths.example, "DATA", "CMTSOLUTION"),
                  "r") as f:
            template = f.readlines()

        print(f"> Writing {self.ntask} event(s) to the SPECFEM3D DATA/")
        for i in range(self.ntask):
            if i == 0:
                lat, lon = lat_c, lon_c
            else:
                angle = 2 * np.pi * (i - 1) / (self.ntask - 1)
                lat, lon = (lat_c + radius * np.sin(angle),
                            lon_c + radius * np.cos(angle))
            name = f"{i + 1:0>3}"

            lines = []
            for line in template:
                if line.startswith("event name"):
                    line = f"event name:       homog_test_{name}\n"
                elif line.startswith("latorUTM"):
                    line = f"latorUTM:       {lat:.1f}\n"
                elif line.startswith("longorUTM"):
                    line = f"longorUTM:      {lon:.1f}\n"
                lines.append(line)
            with open(os.path.join(self.workdir_paths.data,
                                   f"CMTSOLUTION_{name}"), "w") as f:
                f.writelines(lines)

    def write_stations(self):
        """
        Write a STATIONS file with `nsta` stations on a grid of equidistant
        points at the surface. The grid is as square as `nsta` allows (e.g.,
        25 is 5x5) and fills the domain except for a margin along its edges,
        where the absorbing boundaries are. If `nsta` does not fill a grid, the
        last row is incomplete.
        """
        lat_min, lat_max, lon_min, lon_max = self._mesh_extent()
        margin = 0.15  # fraction of the domain
        n_lat = int(np.ceil(np.sqrt(self.nsta)))
        n_lon = int(np.ceil(self.nsta / n_lat))

        def grid(vmin, vmax, n):
            if n == 1:
                return [(vmin + vmax) / 2]
            return np.linspace(vmin + margin * (vmax - vmin),
                               vmax - margin * (vmax - vmin), n)

        lats = grid(lat_min, lat_max, n_lat)
        lons = grid(lon_min, lon_max, n_lon)
        print(f"> Writing {self.nsta} stations on a {n_lat}x{n_lon} grid")
        with open(os.path.join(self.workdir_paths.data, "STATIONS"), "w") as f:
            for i in range(self.nsta):
                f.write(f"S{i + 1:0>3} DB {lats[i // n_lon]:.2f} "
                        f"{lons[i % n_lon]:.2f} 0.0 0.0\n")

    def run_specfem3d_binaries(self, choice, scale=1.):
        """
        Runs the mesher and then the database generator, which write the model
        that SeisFlows takes as the starting or target model, and then moves
        the model into place. Runs with MPI if there is more than one process.

        :type choice: str
        :param choice: 'INIT' or 'TRUE' for the starting or target model
        :type scale: float
        :param scale: factor by which to multiply the velocities of the model
        """
        cd(self.workdir_paths.workdir)

        # The mesher defines the model, so this is where its velocities change
        values = list(self._material)
        values[2] = f"{float(values[2]) * scale:.2f}"  # vp
        values[3] = f"{float(values[3]) * scale:.2f}"  # vs
        self._set_mesh_table_line(
            os.path.join(self.workdir_paths.data, "meshfem3D_files",
                         "Mesh_Par_file"), "NMATERIALS", values
        )
        print(f"> Model {choice}: rho, vp, vs = {values[1]}, {values[2]}, "
              f"{values[3]}")

        rm(self.workdir_paths.output)
        mkdir(os.path.join(self.workdir_paths.output, "DATABASES_MPI"))

        # SPECFEM3D looks for one source, SeisFlows chooses between them later
        source = os.path.join(self.workdir_paths.data, "CMTSOLUTION")
        rm(source)
        ln(f"{source}_001", source)

        mpicmd = f"{self.mpiexec} -n {self.nproc} " if self.nproc > 1 else ""
        # The mesher and database generator, and the file each writes its
        # own messages to, which say more than the MPI messages do
        for exc, log, output in [
                ("xmeshfem3D", "mesher", "output_meshfem3D.txt"),
                ("xgenerate_databases", "database",
                 "output_generate_databases.txt")]:
            log = os.path.join("OUTPUT_FILES", f"{log}.log.txt")
            cmd = f"{mpicmd}bin/{exc} > {log} 2>&1"
            print(f"Running SPECFEM3D with command: {cmd}")
            try:
                subprocess.run(cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                # The program's own message is before MPI's (-----) message
                tail = []
                with open(log, "r") as f:
                    for line in f.read().strip().splitlines():
                        if line.startswith("-----"):
                            break
                        tail.append(line)
                output = os.path.join("OUTPUT_FILES", output)
                if os.path.exists(output):
                    with open(output, "r") as f:
                        tail += f.read().strip().splitlines()[-5:]
                print(msg.cli(f"{exc} has failed, see '{log}' and '{output}'. "
                              f"The end of the output is:",
                              items=tail[-5:] + [str(e)],
                              header="specfem3d error", border="="))
                sys.exit(-1)
        rm(source)

        # SeisFlows expects the model files, and the mesh (for visualization)
        model_dir = self.workdir_paths[f"model_{choice.lower()}"]
        mkdir(model_dir)
        for name in ["vp", "vs", "rho", "x", "y", "z"]:
            cp(glob.glob(os.path.join(self.workdir_paths.output,
                                      "DATABASES_MPI", f"proc*_{name}.bin")),
               model_dir)
        rm(self.workdir_paths.output)

    def setup_seisflows_working_directory(self):
        """
        Create and set the SeisFlows parameter file, making sure all required
        parameters are set correctly for this example problem. The User's
        parameters are applied last, and may override any of the example's.
        """
        cd(self.cwd)

        self.sf.init(force=True)  # Force will delete existing parameter file
        for key, val in self._modules.items():
            self.sf.par(key, val)

        self.sf.configure()

        overrides = {}
        if self.par_yaml:
            overrides.update(load_yaml(self.par_yaml))
        for item in self.par:
            assert("=" in item), f"`--par` expects KEY=VALUE, not '{item}'"
            key, val = item.split("=", 1)
            overrides[key] = val

        # Setting a parameter that is not in the file does nothing, so be sure
        known = {key.lower(): key for key in load_yaml("parameters.yaml")}
        for key in list(self._parameters) + list(overrides):
            if key.lower() not in known:
                print(msg.cli(
                    f"'{key}' is not a parameter of the SeisFlows modules in "
                    f"this example. Did you mean: "
                    f"{difflib.get_close_matches(key.lower(), known) or 'none'}"
                    f". Parameters depend on the system, see `--system`",
                    header="unknown parameter", border="="))
                sys.exit(-1)

        for key, val in {**self._parameters, **overrides}.items():
            self.sf.par(key, val)

    def finalize_specfem3d_par_file(self):
        """
        Final changes to the SPECFEM3D Par_file before running SeisFlows.
        Par_file will be used to control all the child specfem3d directories.
        Need to tell them to read models from .bin files
        """
        print("> Finalizing SPECFEM3D Par_file for SeisFlows example")
        par_file = os.path.join(self.workdir_paths.data, "Par_file")
        setpar("MODEL", "gll", par_file)

        # The mesher is run again by the workflow. The model comes from the
        # files now, but leave the starting model in the mesh files too
        self._set_mesh_table_line(
            os.path.join(self.workdir_paths.data, "meshfem3D_files",
                         "Mesh_Par_file"), "NMATERIALS", self._material
        )

    def run_sf_example(self):
        """
        Use subprocess to run the SeisFlows example we just set up
        """
        print(msg.cli("RUNNING SEISFLOWS EXAMPLE WORKFLOW", border="="))
        cd(self.cwd)
        try:
            subprocess.run("seisflows submit", check=True, shell=True)
        except subprocess.CalledProcessError as e:
            print(msg.cli("EXAMPLE FAILED", items=[str(e)], border="="))
            sys.exit(-1)
        print(msg.cli("EXAMPLE COMPLETED SUCCESFULLY", border="="))

    def main(self):
        """
        Setup the example and then optionally run the actual seisflows workflow
        """
        print(msg.cli("EXAMPLE SETUP", border="="))

        # Step 1: Check SPECFEM3D and make a working directory
        self.check_specfem3d()
        # A new parameter file would pick up the state of an earlier workflow
        # (e.g., its completed tasks and its optimization), not start afresh
        old = [fid for fid in ["sfstate.txt", "scratch", "output"] if
               os.path.exists(os.path.join(self.cwd, fid))]
        if old:
            print(msg.cli(
                f"This directory has the files of an earlier workflow, which "
                f"a new example would continue instead of starting over. Run "
                f"the example in an empty directory, or remove them first "
                f"(`seisflows clean`):", items=old,
                header="existing workflow", border="="))
            sys.exit(-1)
        self.create_specfem3d_working_directory()
        # Step 2: Generate MODEL_INIT and MODEL_TRUE
        print(msg.cli("GENERATING INITIAL MODEL", border="="))
        self.run_specfem3d_binaries(choice="INIT")
        print(msg.cli("GENERATING TRUE/TARGET MODEL", border="="))
        self.run_specfem3d_binaries(choice="TRUE", scale=1 + self.perturb / 100)
        # Step 3: Prepare the parameter files
        self.finalize_specfem3d_par_file()
        self.setup_seisflows_working_directory()
        print(msg.cli("COMPLETE EXAMPLE SETUP", border="="))
        # Step 4: Run the workflwo
        if self.run_example:
            self.run_sf_example()
