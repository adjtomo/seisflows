#!/usr/bin/env python3
"""
A Model class for interfacing with SPECFEM velocity models, used to read
models into memory, manipulate and write to disk. Also contains plotting
functions for SPECFEM2D models
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from glob import glob
from seisflows import logger
from seisflows.tools import unix
from seisflows.tools.config import Dict
from seisflows.tools.math import poissons_ratio
from seisflows.tools.graphics import plot_2d_image
from seisflows.tools.specfem import read_fortran_binary, write_fortran_binary


class Model:
    """
    A container for reading, storing and manipulating model/gradient/kernel
    parameters from SPECFEM2D/3D/3D_GLOBE.
    Stores metadata information alongside model data allowing models to be
    converted to back and forth between vector representations required by the
    optimization library.
    Also contains utility functions to read/write itself so that models can be
    saved alongside their metadata.
    """
    # Dictate the parameters that Model can handle, which does not cover all
    # files created by SPECFEM, which includes things like 'ibool', 'info' etc.
    acceptable_parameters = ["vp", "vs", "rho",
                             "vpv", "vph", "vsv", "vsh", "eta"]
    # Add kernel tag to all acceptable parameters for adjoint simulation results
    acceptable_parameters.extend([f"{_}_kernel" for _ in acceptable_parameters])
    # Add source mask for SPECFEM3D_ source masking
    acceptable_parameters.append("mask_source")
    # Edit acceptable parameters for 3DGLOBE, which must include region name
    for parameter in acceptable_parameters[:]:
        for region in ["1", "2", "3"]:
            acceptable_parameters.append(f"reg{region}_{parameter}")
    

    def __init__(self, path=os.getcwd(), fmt="", parameters=None, regions="123", 
                 flavor=None):
        """
        Model only needs path to model to determine model parameters. Format
        `fmt` can be provided by the user or guessed based on available file
        formats

        .. note::
            The `vector` representation is based completely on the `model`
            attribute. In order to update the model based on vector manipulation
            you must use the update() function.

        :type path: str
        :param path: path to SPECFEM model/kernel/gradient files
        :type fmt: str
        :param fmt: expected format of the files (e.g., '.bin'), if None, will
            attempt to guess based on the file extensions found in `path`
            Available formats are: .bin, .dat
        :type parameters: list
        :param parameters: the list of parameters to consider when defining
            the 'model', which is what will be updated during an inversion. 
        :type regions: str
        :param regions: only for SPECFEM3D_GLOBE, which regions of the chunk 
            to consider in your 'model'. Valid regions are 1, 2 and 3. If you 
            want all regions, set as '123'. If you only want region 1, set as 
            '1'. Order insensitive.
        :type flavor: str
        :param flavor: optional, tell Model what version of SPECFEM was used
            to generate the model, acceptable values are ['2D', '3D', '3DGLOBE']
            If None, will try to guess based on file matching
        """
        self.path = path
        self.fmt = fmt
        self.flavor = flavor
        self.model = None
        if regions:
            self.regions = sorted([f"reg{i}" for i in regions])
        else:
            self.regions = None
        self.coordinates = None
        self._parameters = parameters
        self._ngll = None
        self._nproc = None
    
        # Check that User-provided (optional) flavor matches acceptable values
        if self.flavor is not None:
            acceptable_flavors = ["2D", "3D", "3DGLOBE"]
            assert(self.flavor in acceptable_flavors), \
                f"User-defined `flavor' must be in {acceptable_flavors}"

        # Load an existing model if a valid path is given
        if self.path and os.path.exists(path):
            # Read existing model from a previously saved .npz file
            if os.path.splitext(path)[-1] == ".npz":
                self.model, self.coordinates, self._ngll, self.fmt = \
                    self.load(file=self.path)
                _first_key = list(self.model.keys())[0]
                self._nproc = len(self.model[_first_key])
            # Read a SPECFEM model from its native output files
            else:
                # Dynamically guess things about the model based on files given
                if not self.fmt:
                    self.fmt = self._guess_file_format()
                if not self.flavor:
                    self.flavor = self._guess_specfem_flavor()
    
                # Gather internal representation of the model for manipulation
                self._nproc, self.available_parameters = \
                    self._get_nproc_parameters()
                self.model = self.read(parameters=parameters)

                # Coordinates are only useful for SPECFEM2D models
                if self.flavor == "2D":
                    self.coordinates = self.read_coordinates_specfem2d()

            # .sorted() enforces parameter order every time, otherwise things
            # can get screwy if keys returns different each time
            self._parameters = sorted(self.model.keys())
        else:
            logger.warning("no/invalid `path` given, initializing empty Model")

    def fnfmt(self, i="*", val="*", ext="*"):
        """
        Expected SPECFEM filename format with some checks to ensure that
        wildcards and numbers are accepted. An example filename is:
        'proc000001_vs.bin'

        :type i: int or str
        :param i: processor number or wildcard. If given as an integer, will be
            converted to a 6 digit zero-leading value to match SPECFEM format
        :type val: str
        :param val: parameter value (e.g., 'vs') or wildcard
        :type ext: str
        :param ext: the file format (e.g., '.bin'). If NOT preceded by a '.'
            will have one prepended
        :rtype: str
        :return: filename formatter for use in model manipulation
        """
        if not ext.startswith("."):
            ext = f".{ext}"

        if isinstance(i, int):
            filename_format = f"proc{i:0>6}_{val}{ext}"
        else:
            filename_format = f"proc{i}_{val}{ext}"

        return filename_format

    @property
    def parameters(self):
        """
        Returns a list of parameters which defines the model.
        """
        if not self._parameters:
            self._parameters = sorted(list(self.model.keys()))
        return self._parameters

    @property
    def ngll(self):
        """
        Provide the number of GLL (Gauss Lobatto Legendre) points per processor
        chunk. Access hidden attribute `_ngll` in the case that the model
        is very large, we only want to count the GLL points once.

        :rtype: list of float
        :return: each float represents the number of GLL points for the chunk
            number that corresponds to its index
        """
        if not self._ngll:
            self._update_ngll_from_model()
        return self._ngll

    def _update_ngll_from_model(self):
        """Convenience function to count NGLL points as length of data arrays
        for each parameter processor chunk"""
        self._ngll = []
        for proc in self.model[self.parameters[0]]:
            self._ngll.append(len(proc))

    @property
    def nproc(self):
        """
        Returns the number of processors that define the model/gradient/kernel.

        :rtype: int
        :return: number of processors that define the model
        """
        if not self._nproc:
            self._nproc = len(self.model[self.parameters[0]])
        return self._nproc

    @property
    def vector(self):
        """
        Conveience property to access the merge() function which creates a
        linear vector defining all model parameters

        :rtype: np.array
        :return: a linear vector of all model parameters
        """
        try:
            return self.merge()
        except TypeError as e:
            raise TypeError("Model cannot merge files into continous "
                            "vector") from e

    def copy(self):
        """Returns a deep copy of self so that models can be transferred"""
        return deepcopy(self)

    def read(self, parameters=None):
        """
        Utility function to load in SPECFEM models/kernels/gradients saved in
        various formats. Will try to guess format of model. Assumes that models
        are saved using the following filename format:

        proc{num}_{val}.{format} where `num` is usually a 6 digit number
        representing the processor number (e.g., 000000), `val` is the parameter
        value of the model/kernel/gradient and `format` is the format of the
        file (e.g., bin)

        :type parameters: list of str
        :param parameters: unique parameters to load model for, if None will load
            all available parameters found in `path`
        :rtype: Dict of np arrays
        :return: Dictionary where keys correspond to model parameters and values
            are vectors (np.arrays) representing the model
        """
        if parameters is None:
            parameters = self.available_parameters
        else:
            assert (set(parameters).union(set(self.available_parameters))), (
                f"user-chosen parameters not in available: "
                f"{self.available_parameters}"
            )

        # Pick the correct read function based on the file format
        load_fx = {".bin": self._read_model_fortran_binary,
                   ".dat": self._read_model_ascii,
                   ".adios": self._read_model_adios  # TODO Check if this okay
                   }[self.fmt]

        # Create a dictionary object containing all parameters and their models
        parameter_dict = Dict({key: [] for key in parameters})
        for parameter in parameters:
            parameter_dict[parameter] = load_fx(parameter=parameter)

        return parameter_dict

    def read_coordinates_specfem2d(self):
        """
        Attempt to read coordinate files from the given model definition.
        This is only really useful for SPECFEM2D, where we can plot the
        model, kernel and gradient using matplotlib.

        .. warning
            This will NOT work for SPECEFM3D. When you get to 3D, the
            coordinates don't match up one-to-one with the model values so you
            need external viewing software (e.g., ParaView) to plot.

        :rtype: Dict
        :return: a dictioanary with the X and Z coordinates read in from
            a SPECFEM2D model, if applicable
        """
        coordinates = {"x": [], "z": []}
        if self.fmt == ".bin":
            coordinates["x"] = self._read_model_fortran_binary(parameter="x")
            coordinates["z"] = self._read_model_fortran_binary(parameter="z")
        elif self.fmt == ".dat":
            fids = glob(os.path.join(self.path,
                                     self.fnfmt(val="*", ext=".dat")))
            for fid in sorted(fids):
                coordinates["x"].append(np.loadtxt(fid).T[:, 0])
                coordinates["z"].append(np.loadtxt(fid).T[:, 0])

        # If nothing is found even though we expected files to be there
        if not list(coordinates["x"]) or not list(coordinates["z"]):
            logger.warning("no coordinates found for assumed SPECFEM2D model, "
                           "will not be able to plot figures")
            return None

        # Internal check for parameter validity by checking length of coord
        # against length of model. If they're not the same, then it's not
        # useful to store coordinates because they can't be used for plotting
        assert (len(coordinates["x"]) == len(coordinates["z"])), \
            f"coordinate arrays do not match in length"
        # Assuming all model parameters have the same length
        assert (len(coordinates["x"]) ==
                len(self.model[list(self.model.keys())[0]])), \
            f"length of coordinates array does not match model length"

        return coordinates

    def merge(self, parameter=None):
        """
        Convert dictionary representation `model` to vector representation `m`
        where all parameters and processors are stored as a single 1D vector.
        This vector representation is used by the optimization library during
        model perturbation.

        :type parameter: str
        :param parameter: single parameter to retrieve model vector from,
            otherwise returns all parameters merged into single vector
        :rtype: np.array
        :return: vector representation of the model
        """
        m = np.array([])
        if parameter is None:
            parameters = self.parameters
        else:
            parameters = [parameter]

        for parameter in parameters:
            for iproc in range(self.nproc):
                m = np.append(m, self.model[parameter][iproc])

        return m

    def write(self, path, fmt=None):
        """
        Save a SPECFEM model/gradient/kernel vector loaded into memory back to
        disk in the appropriate format expected by SPECFEM
        """
        unix.mkdir(path)
        if fmt is None:
            assert (self.fmt is not None), f"must specifiy model format: `fmt`"
            fmt = self.fmt

        # Pick the correct read function based on the file format
        save_fx = {".bin": self._write_model_fortran_binary,
                   # ".dat": _write_model_ascii,
                   # ".adios": _write_model_adios   # TODO Check if  right
                   }[fmt]

        save_fx(path=path)

    def split(self, vector=None):
        """
        Converts internal vector representation `m` to dictionary representation
        `model`. Does this by separating the vector based on how it was
        constructed, parameter-wise and processor-wise

        :type vector: np.array
        :param vector: allow Model to split an input vector. If none given,
            will split the internal vector representation
        :rtype: Dict of np.array
        :return: dictionary of model parameters split up by number of processors
        """
        if vector is None:
            vector = self.vector

        model = Dict({key: [] for key in self.parameters})
        for idim, key in enumerate(self.parameters):
            for iproc in range(self.nproc):
                imin = sum(self.ngll) * idim + sum(self.ngll[:iproc])
                imax = sum(self.ngll) * idim + sum(self.ngll[:iproc + 1])
                model[key].extend([vector[imin:imax]])

            model[key] = np.array(model[key])
        return model

    def print_stats(self):
        """
        Print out all model parameters (min, mean, max) in the logger. Useful
        for checking if models are different, if a model is being satisfactorily 
        udpated, or for debugging purposes.
        """
        # Tell the User min and max values of the updated model
        for key, vals in self.model.items():
            min_val = np.hstack(vals).min()
            max_val = np.hstack(vals).max()
            mean_val = np.hstack(vals).mean()
            # Choose formatter based on the magnitude of the value
            if min_val < 1 or max_val > 1E4:
                parts = (f"{key}: min={min_val:.3E}; mean={mean_val:.3E}; "
                         f"max={max_val:.3E}")
            else:
                parts = (f"{key}: min={min_val:.3f}; mean={mean_val:.3f}; "
                         f"max={max_val:.3f}")
            logger.info(parts)

    def check(self, min_pr=-1., max_pr=0.5):
        """
        Checks parameters in the model. If Vs and Vp present, checks poissons
        ratio. Checks for negative velocity values. 

        :type min_pr: float
        :para min_pr: minimum allowable Poisson's ratio, if applicable
        :type max_pr: float
        :param max_pr: maximum allowable Poisson's ratio, if applicable
        :raises AssertionError: 
            - if the input model has no values for any of its parameters
            - if the model contains any NaN values 
        """
        # Checks to make sure the model is filled out, otherwise the following
        # checks will fail unexpectedly
        for key, val in self.model.items():
            # Make sure there are values in the model (not empty)
            assert(val.any()), (
                 f"SPECFEM_{self.flavor} model '{key}' has no values, please "
                 f"check your input model `path_model_init` and the chosen "
                 f"`material` which controls the expected parameters"
                 )
            # Make sure none of the values are NaNs
            assert(not np.isnan(val).any()), (
                 f"SPECFEM_{self.flavor} model '{key}' contains NaN values and "
                 f"should not, please check your model construction"
                 )

        # Check the physicality of the parameters
        if self.flavor in ["2D", "3D"]:
            self._check_2d3d_parameters(min_pr, max_pr)
        elif self.flavor == "3DGLOBE":
            self._check_3dglobe_parameters(min_pr, max_pr)

    def _check_2d3d_parameters(self, min_pr=-1., max_pr=0.5):
        """
        Checks parameters for SPECFEM2D and SPECFEM3D derived models
        """
        if "vs" in self.parameters and "vp" in self.parameters:
            pr = poissons_ratio(vp=self.merge(parameter="vp"),
                                vs=self.merge(parameter="vs"))
            if pr.min() < 0:
                logger.warning("minimum poisson's ratio is negative")
            if pr.max() < min_pr:
                logger.warning(f"maximum poisson's ratio out of bounds: "
                               f"{pr.max():.2f} > {max_pr}")
            if pr.min() > max_pr:
                logger.warning(f"minimum poisson's ratio out of bounds: "
                               f"{pr.min():.2f} < {min_pr}")

        if "vs" in self.model and np.hstack(self.model.vs).min() < 0:
            logger.warning(f"Vs minimum is negative {self.model.vs.min()}")

        if "vp" in self.model and np.hstack(self.model.vp).min() < 0:
            logger.warning(f"Vp minimum is negative {self.model.vp.min()}")

    def _check_3dglobe_parameters(self, min_pr=-1., max_pr=0.5):
        """
        Checks parameters for SPECFEM3D_GLOBE derived models
        """
        # Check for negative velocities
        for tag in ["vsv", "vsh", "vph", "vpv", "vp", "vs"]:
            for reg in self.regions:
                par = f"{reg}_{tag}"
                if par in self.model and \
                        np.hstack(self.model[par]).min() < 0:
                    logger.warning(f"{par} minimum for is negative "
                                   f"{self.model['par'].min()}")

        # Check Poisson's ratio for all (an)isotropic velocity combinations
        for vs_par in ["vs", "vsv", "vsh"]:
            for vp_par in ["vp", "vpv", "vph"]:
                for reg in self.regions:
                    vs_par = f"{reg}_{vs_par}"  # e.g., 'reg1_vsv'
                    vp_par = f"{reg}_{vp_par}"
                    if vs_par in self.parameters and vp_par in self.parameters:
                        pr = poissons_ratio(vp=self.merge(parameter=vp_par),
                                            vs=self.merge(parameter=vs_par))
                        if pr.min() < 0:
                            logger.warning(f"minimum {vp_par}, {vs_par} "
                                           f"poisson's ratio is negative")
                        if pr.max() < min_pr:
                            logger.warning(f"maximum {vp_par}, {vs_par} "
                                           f"poisson's ratio out of bounds: "
                                           f"{pr.max():.2f} > {max_pr}")
                        if pr.min() > max_pr:
                            logger.warning(f"minimum {vp_par}, {vs_par} "
                                           f"poisson's ratio out of bounds: "
                                           f"{pr.min():.2f} < {min_pr}")

    def save(self, path):
        """
        Save instance attributes (model, vector, metadata) to disk as an
        .npz array so that it can be loaded in at a later time for future use
        """
        model = self.split()
        if self.coordinates:
            # Incase we have model parameters called 'x' or 'z', rename for save
            model["x_coord"] = self.coordinates["x"]
            model["z_coord"] = self.coordinates["z"]

        np.savez(file=path, fmt=self.fmt, **model)

    def _load2d3d(self, file):
        """
        Load in a previously saved .npz file containing model information
        and re-create a Model instance matching the one that was `save`d

        :type file: str
        :param file: .npz file to load data from. Must have been created by
            Model.save()
        :rtype: tuple (Dict, list, str)
        :return: (Model Dictionary, ngll points for each slice, file format)
        """
        model = Dict()
        coords = Dict()
        ngll = []
        data = np.load(file=file)
        for i, key in enumerate(data.files):
            if key == "fmt":
                continue
            # Special case where we are using SPECFEM2D and carry around coords
            elif "coord" in key:
                coords[key[0]] = data[key]  # drop '_coord' suffix from `save`
            else:
                model[key] = data[key]
                # Assign the number of GLL points per slice. Only needs to happen
                # once because all model values should have same points/slice
                if not ngll:
                    for array in model[key]:
                        ngll.append(len(array))

        return model, coords, ngll, str(data["fmt"])

    def load(self, file):
        """
        Load in a previously saved .npz file containing model information
        and re-create a Model instance matching the one that was `save`d

        :type file: str
        :param file: .npz file to load data from. Must have been created by
            Model.save()
        :rtype: tuple (Dict, list, str)
        :return: (Model Dictionary, ngll points for each slice, file format)
        """
        model = Dict()
        coords = Dict()
        ngll = []
        data = np.load(file=file)
        for i, key in enumerate(data.files):
            if key == "fmt":
                continue
            # Special case where we are using SPECFEM2D and carry around coords
            elif "coord" in key:
                coords[key[0]] = data[key]  # drop '_coord' suffix from `save`
            else:
                model[key] = data[key]
                # Assign the number of GLL points per slice. Only needs to happen
                # once because all model values should have same points/slice
                if not ngll:
                    for array in model[key]:
                        ngll.append(len(array))

        return model, coords, ngll, str(data["fmt"])

    def update(self, model=None, vector=None):
        """
        Update internal model/vector defitions. Because these two quantities
        are tied to one another, updating one will update the other. This
        function simply makes that easier.
        """
        if model is not None:
            self.model = model
        elif vector is not None:
            self.model = self.split(vector=vector)

    def plot2d(self, parameter, cmap=None, show=True, title="", save=None):
        """
        Plot internal model parameters as a 2D image plot.

        .. warning::
            This is only available for SPECFEM2D models. SPECFEM3D model
            coordinates do not match the model vectors (because the grids are
            irregular) and cannot be visualized like this.

        :type parameter: str
        :param parameter: chosen internal parameter value to plot.
        :type cmap: str
        :param cmap: colormap which match available matplotlib colormap.
            If None, will choose default colormap based on parameter choice.
        :type show: bool
        :param show: show the figure after plotting
        :type title: str
        :param title: optional title to prepend to some values that are
            provided by default. Useful for adding information about iteration,
            step count etc.
        :type save: str
        :param save: if not None, full path to figure to save the output image
        """
        assert (parameter in self._parameters), \
            f"chosen `parameter` must be in {self._parameters}"

        assert (self.coordinates is not None), (
            f"`plot2d` function requires model coordinates which are only "
            f"available for solver SPECFEM2D"
        )

        # Choose default colormap based on parameter values
        if cmap is None:
            if "kernel" in parameter:
                cmap = "seismic_r"
                zero_midpoint = True
            else:
                cmap = "Spectral"
                zero_midpoint = False

        # 'Merge' the coordinate matrices to get a vector representation
        x, z = np.array([]), np.array([])
        for iproc in range(self.nproc):
            x = np.append(x, self.coordinates["x"][iproc])
            z = np.append(z, self.coordinates["z"][iproc])
        data = self.merge(parameter=parameter)

        f, p, cbar = plot_2d_image(x=x, z=z, data=data, cmap=cmap,
                                   zero_midpoint=zero_midpoint)

        # Set some figure labels based on information we know here
        ax = plt.gca()
        ax.set_xlabel("X [m]")
        ax.set_ylabel("Z [m]")
        # Allow User to title the figure, e.g., with iteration, step count etc.
        _title = (f"{parameter.title()}_min = {data.min()};\n"
                  f"{parameter.title()}_max = {data.max()};\n"
                  f"{parameter.title()}_mean = {data.mean()}")
        if title:
            _title = f"{title}\n{_title}"
        ax.set_title(_title)
        cbar.ax.set_ylabel(parameter.title(), rotation=270, labelpad=15)

        if save:
            plt.savefig(save)
        if show:
            plt.show()

    def _get_nproc_parameters(self):
        """
        Get the number of processors and the available parameters from a list of
        output SPECFEM model files.

        :rtype: tuple (int, list)
        :return: (number of processors, list of available parameters in dir)
        """
        fids = glob(os.path.join(self.path, self.fnfmt(val="*", ext=self.fmt)))
        fids = [os.path.basename(_) for _ in fids]  # drop full path
        fids = [os.path.splitext(_)[0] for _ in fids]  # drop extension

        if self.fmt == ".bin":
            avail_par = list(set(["_".join(_.split("_")[1:]) for _ in fids]))
            # Remove any parameters not accepted by Model
            avail_par = list(set(avail_par).intersection(
                                        set(self.acceptable_parameters)
                                        ))
            # Count the number of files for matching parameters only (do once)
            # Globe version requires the region number in the wild card
            nproc = len(glob(os.path.join(
                self.path, self.fnfmt(val=avail_par[0], ext=self.fmt)))
            )
        elif self.fmt == ".dat":
            # e.g., 'proc000000_rho_vp_vs'
            _, *avail_par = fids[0].split("_")
            # Remove any parameters not accepted by Model
            avail_par = list(set(avail_par).intersection(
                                        set(self.acceptable_parameters)
                                        ))
            nproc = len(fids) + 1
        else:
            raise NotImplementedError(f"{self.fmt} is not yet supported by "
                                      f"SeisFlows")


        return nproc, list(avail_par)

    def _guess_file_format(self):
        """
        Guess the file format of model/kernel/gradient files if none provided by
        the user. Does so by checking file formats against formats expected from
        SPECFEM2D/3D/3D_GLOBE

        :rtype: str
        :return: file format suffix with a leading '.' e.g., '.bin'
        """
        acceptable_formats = {".bin", ".dat"}

        files = glob(os.path.join(self.path, "*"))
        suffixes = set([os.path.splitext(_)[1] for _ in files])
        fmt = acceptable_formats.intersection(suffixes)
        assert (len(fmt) == 1), (
            f"cannot guess model format, multiple matching acceptable formats "
            f"found: {list(suffixes)}"
        )
        return list(fmt)[0]  # pulling single entry from set

    def _guess_specfem_flavor(self):
        """
        Guess if SPECFEM2D/3D/3D_GLOBE was used to generate the model based
        on the format of the file. Check based on unique available files or 
        properties

        :rtype: str
        :return: SPECFEM flavor, one of ['2D', '3D', '3DGLOBE']
        """
        fullpaths = glob(os.path.join(self.path, f"*{self.fmt}"))
        assert fullpaths, f"cannot find files for flavor guessing"

        # Not the most accurate way of doing this, but serves a purpose
        fids = [os.path.basename(_) for _ in fullpaths]
        fids = [os.path.splitext(_)[0] for _ in fids]
        unique_tags = set(["_".join(_.split("_")[1:]) for _ in fids])
  
        # SPECFEM3D_GLOBE
        # Smash all the tags into a single string and look for 'reg1' (or 
        # whatever region User chooses). Assuming here that SPECFEM2D/3D won't
        # have parameters that contain the phrase 'reg1'
        if self.regions and self.regions[0] in "".join(unique_tags):
            flavor = "3DGLOBE"
        # SPECFEM3D is the only one that has a 'y' parameter, globe code
        # doesn't store coordinate information and 2D only has X and Z coords
        elif "y" in unique_tags:
            flavor = "3D"
        else:
            flavor = "2D"

        return flavor

    def _read_model_fortran_binary(self, parameter):
        """
        Load Fortran binary models into disk. This is the preferred model format
        for SeisFlows <-> SPECFEM interaction

        :type parameter: str
        :param parameter: chosen parameter to load model for
        :rtype: np.array
        :return: vector of model values for given `parameter`
        """
        array = []
        fids = glob(os.path.join(
            self.path, self.fnfmt(val=parameter, ext=".bin"))
        )
        for fid in sorted(fids):  # make sure were going in numerical order
            array.append(read_fortran_binary(fid))

        # !!! Causes a visible deprecation warning from NumPy but setting
        # !!! array type as 'object' causes problems with pickling and
        # !!! merging arrays
        array = np.array(array)

        return array

    def _read_model_adios(self, parameter):
        """
        Load ADIOS models into disk

        :type parameter: str
        :param parameter: chosen parameter to load model for
        :rtype: np.array
        :return: vector of model values for given `parameter`
        """
        raise NotImplementedError("ADIOS file formats are not currently "
                                  "implemented into SeisFlows")

    def _read_model_ascii(self, parameter):
        """
        Load ASCII SPECFEM2D models into disk. ASCII models are generally saved
        all in a single file with all parameters together as a N column ASCII 
        file where columns 1 and 2 are the coordinates of the mesh, and the 
        remainder columns are data corresponding to the filenames
        e.g., proc000000_rho_vp_vs.dat, rho is column 3, vp is 4 etc.

        :type parameter: str
        :param parameter: chosen parameter to load model for
        :rtype: np.array
        :return: vector of model values for given `parameter`
        """
        fids = glob(os.path.join(self.path, self.fnfmt(val="*", ext=".dat")))
        _, *available_parameters = fids[0].split("_")
        assert (parameter in available_parameters), (
            f"{parameter} not available for ASCII model"
        )
        # +2 because first 2 columns are the X and Z coordinates in the mesh
        array = []
        column_idx = available_parameters.index(parameter) + 2
        for fid in sorted(fids):
            array.append(np.loadtxt(fid).T[:, column_idx])

        return np.array(array)

    def _write_model_fortran_binary(self, path):
        """
        Save a SPECFEM model back to Fortran binary format.
        Data are written as single precision floating point numbers

        .. note::
            FORTRAN unformatted binaries are bounded by an INT*4 byte count.
            This function mimics that behavior by tacking on the boundary data
            as 'int32' at the top and bottom of the data array.
            https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnc4/index.html
        """
        for parameter in self.parameters:
            for i, data in enumerate(self.model[parameter]):
                filename = self.fnfmt(i=i, val=parameter, ext=".bin")
                filepath = os.path.join(path, filename)
                write_fortran_binary(arr=data, filename=filepath)

