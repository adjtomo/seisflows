#!/usr/bin/env python3
"""
A class for interfacing with SPECFEM velocity models, used to manipulate SPECFEM 
binary models. Processing is done in parallel to keep things efficient

.. TODO::
    - Add support for SPECFEM2D ASCII models
    - Determine if we need a way to filter by parameter, currently we just 
      treat the entire model as a vector of values, which is fine but sometimes
      we may want to only deal with a single parameter?
    - Determine how to efficiently take a dot product?
"""
import os
import numpy as np
from concurrent.futures import ProcessPoolExecutor, wait
from glob import glob
from seisflows.tools import unix
from seisflows.tools.specfem import read_fortran_binary, write_fortran_binary


class Model:
    """
    A container for manipulating models/kernels in SPECFEM2D/3D/3D_GLOBE
    """
    # Dictate the parameters that Model can handle, which does not cover all
    # files created by SPECFEM, which includes things like 'ibool', 'info' etc.
    acceptable_parameters = [
            # Isotropic wavespeeds (kernels alpha and beta are renamed vp, vs)
            "vp", "vs",
            # Transverse Isotropic Wavespeeds
            "vpv", "vph", "vsv", "vsh", "eta",
            # Density, attenuation
            "rho", "rhop", "kappa", "mu",
            # 21 Parameter general anisotropy c_ij
            "c11", "c12", "c13", "c14", "c15", "c16", 
            "c22", "c23", "c24", "c25", "c26", 
            "c33", "c34", "c35", "c36", 
            "c44", "c45", "c46", 
            "c55", "c56", 
            "c66", 
            ]
    # Add kernel tag to all acceptable parameters for adjoint simulation results
    acceptable_parameters.extend([f"{_}_kernel" for _ in acceptable_parameters])
    # Add source mask for SPECFEM3D_ source masking
    acceptable_parameters.append("mask_source")
    # Edit acceptable parameters for 3DGLOBE, which must include region name
    for parameter in acceptable_parameters[:]:
        for region in ["1", "2", "3"]:
            acceptable_parameters.append(f"reg{region}_{parameter}")
    # Define acceptable actions for Model.apply()
    acceptable_actions = ["*", "/", "+", "-"]


    def __init__(self, path, fmt="", parameters=None, regions=None, 
                 parallel=True):
        """
        Model only needs path to model to determine model parameters. Format
        `fmt` can be provided by the user or guessed based on available file
        formats

        :type path: str
        :param path: required, path to SPECFEM model/kernel/gradient files
        :type fmt: str
        :param fmt: expected format of the files (e.g., '.bin'). By default the 
            Model class will attempt to guess based on the file extensions found 
            in `path` (`fmt`==NoneType). Available formats are: 
            - .bin: Fortran binary format
            - .dat: ASCII-based data format used for SPECFEM2D
        :type parameters: list
        :param parameters: optional list of parameters to consider when defining
            the 'model', which is what will be updated during an inversion. If 
            not given, all parameters found in `path` that match 
            `acceptable_parameters` will be used to define the model.
        :type regions: str
        :param regions: SPECFEM3D_GLOBE ONLY! Defines which regions of the chunk 
            to consider in your 'model'. Valid regions are 1, 2 and 3. If you 
            want all regions, set as '123'. If you only want region 1, set as 
            '1'. Order insensitive.
        :type parallel: bool
        :param parallel: run tasks in parallel using Concurrent.futures
        :raises AssertionError: if the `path` does not exist or the `path`
            does not contain expected Model files
        """
        # If the User gives a `path`, we EXPECT model files 
        assert(os.path.exists(path)), \
            f"Model `path` does not exist, cannot init model: '{path}'"
        
        self.path = os.path.abspath(path)
        self.parameters = parameters
        self.parallel = parallel
        self.regions = regions
        # Used to build filenames for Globe which requires 'reg' prefix
        if self.regions is not None:
            self.regions = sorted([f"reg{i}" for i in regions])

        # Guess model format and flavor or have User define it
        self.fmt = fmt or self._guess_file_format()

        # Requires that `fmt` be defined before running
        self.nproc = self._get_nproc()

        # A few checks to make sure things are acceptable
        for par in self.parameters:
            assert(par in self.acceptable_parameters), \
                (f"parameter `{par}` does not match `acceptable_parameters`")

        # Check that required files exist within the `path` for each parameter
        for par in self.parameters:
            assert(glob(os.path.join(path, self._get_filename(val=par, 
                                                             ext=self.fmt)))), \
                f"Model `path` contains no match for parameter {par}"
            
        # Generate filename list. Ensure it is sorted so that other models with
        # the same parameters will have the same indexing
        self.filenames = []
        for par in sorted(self.parameters):
            for iproc in range(self.nproc):
                filename = self._get_filename(i=iproc, val=par, ext=self.fmt)
                self.filenames.append(os.path.join(self.path, filename))
    
    def read(self, filename):
        """
        Read model file based on the particular format specified at init

        :type fid: str
        :param fid: path to file to be read in
        :rtype: np.array
        :return: model values as a numpy vector
        """
        if self.fmt == ".bin":
            arr = read_fortran_binary(filename)
        elif self.fmt == ".dat":
            raise NotImplementedError
        else:
            raise NotImplementedError

        return arr

    def write(self, arr, filename):
        """
        Save a SPECFEM model/gradient/kernel vector loaded into memory back to
        disk in the appropriate format expected by SPECFEM

        :type arr: numpy.array
        :param arr: array to be written out
        :type filename: str
        :param filename: full path and filename for the file that will be 
            written out containing the data in `arr`
        """
        if not os.path.exists(os.path.dirname(filename)):
            unix.mkdir(os.path.dirname(filename))

        if self.fmt == ".bin":
            write_fortran_binary(arr=arr, filename=filename)
        else:
            raise NotImplementedError
    
    def get(self, what):
        """
        Find some value in the model based on `what` you want

        :type what: str
        :param what: the value to find in the model, one of:
            - 'max': maximum value in the model
            - 'min': minimum value in the model
            - 'sum': sum of all values in the model
        """
        if self.parallel:
            # !!! CHECK THIS I don't think it preserves order or even works
            with ProcessPoolExecutor(max_workers=unix.nproc()) as executor:
                futures = [
                    executor.submit(self._get_single, fid, what)
                                    for fid in self.filenames
                                    ]
                wait(futures)     
            vals = [f.result() for f in futures]              
        else:
            vals = []
            for fid in self.filenames:
                vals.append(self._get_single(fid, what))

        # Combine all the values returned from the indvidual files
        if what == "max":
            return np.max(vals)
        elif what == "min":
            return np.min(vals)
        elif what == "sum":  # !!! CHECK THIS
            return np.sum(vals)

    def _get_single(self, filename, what):
        """
        Get a single value from a single file in the model, used by `get` to
        """    
        if what == "max":
            return np.max(self.read(filename))
        elif what == "min":
            return np.min(self.read(filename))
        elif what == "mean":
            return np.mean(self.read(filename))
        elif what == "sum":
            return np.sum(self.read(filename))
        else:
            raise ValueError(f"Unknown action {what} for Model.get()")

    def dot(self, other):
        """
        Dot product of model with `other` model. Performs operations in serial
        or parallel and returns a scalar output
        """
        pass

    def apply(self, actions, values=None, other=None, export_to=None):
        """
        Main model function which manipulates the model with a given list of 
        `actions` and associated `values`. `actions` can also be performed
        with `other` Models, such as dot products or element-wise processes. 
        Lists of actions and values allow chaining of actions in one process
        so that everything can be done in RAM rather than requiring multiple
        reads and writes to disk.

        .. note::

            If using `other`, the FIRST action in `actions` will be applied to
            the other model, and all other `actions` will be applied to `values`    
            such that len(actions) == len(values) + 1. For example:

            Model.apply(actions=["*", "+"], values=[1], other=OtherModel)

            will multiply Model and OtherModel and then add 1 to the result.

        .. note::

            If NOT using `other`, then `actions` must match `values` in length.
            For example:

            Model.apply(actions=["multiply", "add"], values=[2, 1])

            Will multiply the model by 2 and then add 1 to the result.

        .. note::

            Subtraction and division are one-way, i.e., A - B or A / B only. If
            you want to go the other way, then call this function with the 
            other Model object.

        :type other: seisflows.tools.specfem_model.SpecfemModel
        :param other: the `other` model to loop over. Can have different 
            parameters (e.g., model and kernel) but must match format and
            nproc in order to loop with model
        :type action: list of str
        :param action: action to perform using the `actions` functions. 
            If using `other` model, then the first action will be applied and 
            can be one of the following:
            - '*': multiply the model by the other model values
            - '/': divide the model by the `other` model values
            - '+': add the `other` model values to the model
            - '-': subtract the `other` model values from the model
            If not using `other` model, then the actions can only be `multiply` 
            and `add`
        :type values: list of float
        :param values: list of values to use in the `actions`. Must match the
            length of `actions`. 
        :type export_to: str
        :param export_to: path to export the model values to disk. The filenames
            of the exported files will match the input filenames of this Model,
            not the `other` model. If not given, will export to the same
            directory as the input model `self.path`, which means it overwrites
            the current model
        """
        # Default values is empty list so it can be compared
        if values is None:
            values = []
            
        # If applying with another model, ensure matching Model characteristics
        if other:
            assert(self.fmt == other.fmt), (
                f"Model format {self.fmt} does not match other model "
                f"{other.fmt}")
            assert(self.nproc == other.nproc), (
                f"Model nproc {self.nproc} do not match other model "
                f"{other.nproc}")
            assert(len(actions) == len(values) + 1), (
                f"actions must be 1 longer than values if using `other` model, "
                f"because first action is applied to `other` model"
            )
        else:
            assert(len(actions) == len(values)), (
                f"actions must match values in length"
            )
        
        # Check that actions are acceptable
        for action in actions:
            assert(action in self.acceptable_actions), (
                f"action must be in {self.acceptable_actions}, not {action}"
            )

        # Set export location, defaults to overwriting current model 
        if export_to is None:
            export_to = self.path
            
        # If no `other` model is given, then we only deal with the current model
        if other is None:
            if self.parallel:
                with ProcessPoolExecutor(max_workers=unix.nproc()) as executor:
                    futures = [
                        executor.submit(self._apply_single, fid, actions, 
                                        values, export_to) 
                                        for fid in self.filenames
                                        ]
                    wait(futures)                   
            else:
                for fid in self.filenames:
                    self._apply_single(fid, actions, values, export_to)
                    
        # If `other` model is given, then we loop over the two models together
        else:
            if self.parallel:
                with ProcessPoolExecutor(max_workers=unix.nproc()) as executor:
                    futures = [
                        executor.submit(self._apply_with_single, fid, other_fid, 
                                        actions, values, export_to) 
                                        for fid, other_fid in zip(self.filenames, 
                                                                other.filenames)
                                        ]
                    wait(futures)                   
            else:
                for fid, other_fid in zip(self.filenames, other.filenames):
                    self._apply_with_single(fid, other_fid, actions, values, 
                                            export_to)

    def _apply_single(self, fid, actions, values, export_to):
        """
        Apply an action to a single parameter and processor, this should be 
        parallelized by a main calling function and not called individually

        .. note::
        
            - Multiplying with NumPy is the same speed as np.multiply
              https://stackoverflow.com/questions/49493482/\
                  numpy-np-multiply-vs-operator
    
        :type fid: str
        :param fid: full path to filename to read and maninpulate
        :type action: str
        :param action: action to perform using the `actions` functions
            - '*': multiply the model by `value`. If you want to divide
              then `multiply` by `1/value`
            - '+': add `value` to the model. If you want to subtract
                then `add` by `-1*value`
            
        :type value: float
        :param value: value to use in the `action`
        :type export_to: str
        :param export_to: path to the directory to export the action'ed file.
            filename will match the input filename of `fid`. If None then
            will export to the same path as the input model `self.path`, which
            is an inplace action
        """
        arr = self.read(filename=fid)

        for action, value in zip(actions, values):
            if action == "*": 
                arr *= value
            elif action == "+":
                arr += value    
            # > ADD ADDITIONAL ACTIONS HERE

        fid_out = os.path.basename(fid)
        self.write(arr=arr, filename=os.path.join(export_to, fid_out))
                                           
    def _apply_with_single(self, fid, other_fid, actions, values, 
                           export_to=None):
        """
        Apply an action to a single parameter and processor, should be 
        parallelized by a main calling function and not called individually
    
        :type fid: str
        :param fid: full path to filename to read and maninpulate
        :type action: str
        :param action: action to perform using the `actions` functions
        :type value: float
        :param value: value to use in the `action`
        :type export_to: str
        :param export_to: path to the directory to export the action'ed file.
            filename will match the input filename of `fid`
        """
        arr = self.read(filename=fid)
        other_arr = self.read(filename=other_fid)

        other_action, *actions = actions  # first action applied to other model

        # Apply the first action to both models 
        if other_action == "*": 
            arr *= other_arr
        elif other_action == "/": 
            arr /= other_arr
        elif other_action == "+":
            arr += other_arr    
        elif other_action == "-":
            arr -= other_arr

        # This is repeated from `_apply_single` but it's short enough.
        # This will get skipped if there are no other `actions`` to apply
        for action, value in zip(actions, values):
            if action == "*": 
                arr *= value
            elif action == "+":
                arr += value    
            # > ADD ADDITIONAL ACTIONS HERE

        fid_out = os.path.basename(fid)
        self.write(arr=arr, filename=os.path.join(export_to, fid_out))

    def _get_filename(self, i="*", val="*", ext="*"):
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
                        
    def _get_nproc(self):
        """
        Get the number of processors and the available parameters from a list of
        output SPECFEM model files.

        :rtype: tuple (int, list)
        :return: (number of processors, sorted list of available parameters)
        """
        if self.fmt == ".bin":
            # Count the number of files for matching parameters only (do once)
            # Globe version requires the region number in the wild card
            par = self.parameters[0]
            nproc = len(glob(
                os.path.join(self.path, self._get_filename(val=par, 
                                                          ext=self.fmt)))
            )
        elif self.fmt == ".dat":
            # e.g., files contain all parameters e.g., 'proc000000_rho_vp_vs'
            fids = glob(os.path.join(self.path, 
                                     self._get_filename(val="*", ext=self.fmt))
                                     )
            nproc = len(fids) + 1  # indexing starts at 0
        else:
            raise NotImplementedError(f"{self.fmt} is not yet supported by "
                                      f"SeisFlows")

        return nproc
    
    def _get_available_parameters(self):
        """
        Get a list of available parameters in the `path`. This might be useful
        for searching for specific model parameters.

        .. note::

            Currently not used, but may be helpful in the future
        """
        fids = glob(os.path.join(self.path, 
                                 self._get_filename(val="*", ext=self.fmt))
                                 )
        fids = [os.path.basename(_) for _ in fids]  # drop full path
        fids = [os.path.splitext(_)[0] for _ in fids]  # drop extension

        if self.fmt == ".bin":
            avail_par = list(set(["_".join(_.split("_")[1:]) for _ in fids]))
            # Remove any parameters not accepted by Model
            avail_par = list(set(avail_par).intersection(
                                        set(self.acceptable_parameters)
                                        ))
        elif self.fmt == ".dat":
            # e.g., 'proc000000_rho_vp_vs'
            _, *avail_par = fids[0].split("_")
            # Remove any parameters not accepted by Model
            avail_par = list(set(avail_par).intersection(
                                        set(self.acceptable_parameters)
                                        ))
        else:
            raise NotImplementedError(f"{self.fmt} is not yet supported by "
                                      f"SeisFlows")
        
        return sorted(list(avail_par))

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

    # SPECFEM2D-specific functions below
    # def _read_coordinates_specfem2d(self):
    #     """
    #     Attempt to read coordinate files from the given model definition.
    #     This is only really useful for SPECFEM2D, where we can plot the
    #     model, kernel and gradient using matplotlib.

    #     .. warning
    #         This will NOT work for SPECEFM3D. When you get to 3D, the
    #         coordinates don't match up one-to-one with the model values so you
    #         need external viewing software (e.g., ParaView) to plot.

    #     :rtype: Dict
    #     :return: a dictioanary with the X and Z coordinates read in from
    #         a SPECFEM2D model, if applicable
    #     """
    #     coordinates = {"x": [], "z": []}
    #     if self.fmt == ".bin":
    #         coordinates["x"] = self._read_model_fortran_binary(parameter="x")
    #         coordinates["z"] = self._read_model_fortran_binary(parameter="z")
    #     elif self.fmt == ".dat":
    #         fids = glob(os.path.join(self.path,
    #                                  self._get_filename(val="*", ext=".dat")))
    #         for fid in sorted(fids):
    #             coordinates["x"].append(np.loadtxt(fid).T[:, 0])
    #             coordinates["z"].append(np.loadtxt(fid).T[:, 0])

    #     # If no coordinates then move on, sometimes we don't save the 
    #     # coordinates if we're just using the model for updates. But without
    #     # coordinates the default plot functions will not work
    #     if not list(coordinates["x"]) or not list(coordinates["z"]):
    #         return None

    #     # Internal check for parameter validity by checking length of coord
    #     # against length of model. If they're not the same, then it's not
    #     # useful to store coordinates because they can't be used for plotting
    #     assert (len(coordinates["x"]) == len(coordinates["z"])), \
    #         f"coordinate arrays do not match in length"
    #     # Assuming all model parameters have the same length
    #     assert (len(coordinates["x"]) ==
    #             len(self.model[list(self.model.keys())[0]])), \
    #         f"length of coordinates array does not match model length"

    #     return coordinates

    # def plot2d(self, parameter, cmap=None, show=True, title="", save=None):
    #     """
    #     TODO Consider moving this into a plotting utility
    #     Plot internal model parameters as a 2D image plot.

    #     .. warning::
    #         This is only available for SPECFEM2D models. SPECFEM3D model
    #         coordinates do not match the model vectors (because the grids are
    #         irregular) and cannot be visualized like this.

    #     :type parameter: str
    #     :param parameter: chosen internal parameter value to plot.
    #     :type cmap: str
    #     :param cmap: colormap which match available matplotlib colormap.
    #         If None, will choose default colormap based on parameter choice.
    #     :type show: bool
    #     :param show: show the figure after plotting
    #     :type title: str
    #     :param title: optional title to prepend to some values that are
    #         provided by default. Useful for adding information about iteration,
    #         step count etc.
    #     :type save: str
    #     :param save: if not None, full path to figure to save the output image
    #     """
    #     assert (parameter in self._parameters), \
    #         f"chosen `parameter` must be in {self._parameters}"

    #     assert (self.coordinates is not None), (
    #         f"`plot2d` function requires model coordinates which are only "
    #         f"available for solver SPECFEM2D"
    #     )

    #     # Choose default colormap based on parameter values
    #     if cmap is None:
    #         if "kernel" in parameter:
    #             cmap = "seismic_r"
    #             zero_midpoint = True
    #         else:
    #             cmap = "Spectral"
    #             zero_midpoint = False

    #     # 'Merge' the coordinate matrices to get a vector representation
    #     x, z = np.array([]), np.array([])
    #     for iproc in range(self.nproc):
    #         x = np.append(x, self.coordinates["x"][iproc])
    #         z = np.append(z, self.coordinates["z"][iproc])
    #     data = self.merge(parameter=parameter)

    #     f, p, cbar = plot_2d_image(x=x, z=z, data=data, cmap=cmap,
    #                                zero_midpoint=zero_midpoint)

    #     # Set some figure labels based on information we know here
    #     ax = plt.gca()
    #     ax.set_xlabel("X [m]")
    #     ax.set_ylabel("Z [m]")
    #     # Allow User to title the figure, e.g., with iteration, step count etc.
    #     _title = (f"{parameter.title()}_min = {data.min()};\n"
    #               f"{parameter.title()}_max = {data.max()};\n"
    #               f"{parameter.title()}_mean = {data.mean()}")
    #     if title:
    #         _title = f"{title}\n{_title}"
    #     ax.set_title(_title)
    #     cbar.ax.set_ylabel(parameter.title(), rotation=270, labelpad=15)

    #     if save:
    #         plt.savefig(save)
    #     if show:
    #         plt.show()

 
    def check(self, arr, min_pr=-1., max_pr=0.5):
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
        assert(arr.size != 0), (
                 f"Model has no values, please "
                 f"check your input model `path_model_init` and the chosen "
                 f"`material` which controls the expected parameters"
                 )
        # Make sure none of the values are NaNs
        assert(not np.isnan(arr)), (
                f"Model contains NaN values and "
                f"should not, please check your model construction"
                )
            
        # # SPECFEM 2D and 3D parameters
        # if "vs" in self.parameters and "vp" in self.parameters:
        #     pr = poissons_ratio(vp=self.merge(parameter="vp"),
        #                         vs=self.merge(parameter="vs"))
        #     if pr.min() < 0:
        #         logger.warning("minimum poisson's ratio is negative")
        #     if pr.max() < min_pr:
        #         logger.warning(f"maximum poisson's ratio out of bounds: "
        #                        f"{pr.max():.2f} > {max_pr}")
        #     if pr.min() > max_pr:
        #         logger.warning(f"minimum poisson's ratio out of bounds: "
        #                        f"{pr.min():.2f} < {min_pr}")

        # if "vs" in self.model and np.hstack(self.model.vs).min() < 0:
        #     logger.warning(f"Vs minimum is negative {self.model.vs.min()}")

        # if "vp" in self.model and np.hstack(self.model.vp).min() < 0:
        #     logger.warning(f"Vp minimum is negative {self.model.vp.min()}")


        # # SPECFEM3D_GLOBE: Check Poisson's ratio for all (an)isotropic velocity 
        # check = True
        # for par in ["vs", "vsv", "vsh", "vp", "vpv", "vph"]:
        #     if par not in self.parameters:
        #         check = False
        # if check:
        #     for tag in ["vsv", "vsh", "vph", "vpv", "vp", "vs"]:
        #         for reg in self.regions:
        #             par = f"{reg}_{tag}"
        #             if par in self.model and \
        #                     np.hstack(self.model[par]).min() < 0:
        #                 logger.warning(f"{par} minimum for is negative "
        #                             f"{self.model['par'].min()}")
        #     for vs_par in ["vs", "vsv", "vsh"]:
        #         for vp_par in ["vp", "vpv", "vph"]:
        #             for reg in self.regions:
        #                 vs_par = f"{reg}_{vs_par}"  # e.g., 'reg1_vsv'
        #                 vp_par = f"{reg}_{vp_par}"
        #                 if vs_par in self.parameters and \
        #                     vp_par in self.parameters:
        #                     pr = poissons_ratio(vp=self.merge(parameter=vp_par),
        #                                         vs=self.merge(parameter=vs_par))
        #                     if pr.min() < 0:
        #                         logger.warning(f"minimum {vp_par}, {vs_par} "
        #                                     f"poisson's ratio is negative")
        #                     if pr.max() < min_pr:
        #                         logger.warning(f"maximum {vp_par}, {vs_par} "
        #                                     f"poisson's ratio out of bounds: "
        #                                     f"{pr.max():.2f} > {max_pr}")
        #                     if pr.min() > max_pr:
        #                         logger.warning(f"minimum {vp_par}, {vs_par} "
        #                                     f"poisson's ratio out of bounds: "
        #                                     f"{pr.min():.2f} < {min_pr}")


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
        fids = glob(os.path.join(self.path, 
                                 self._get_filename()(val="*", ext=".dat")))
        _, *available_parameters = fids[0].split("_")
        assert (parameter in available_parameters), (
            f"{parameter} not available for ASCII model"
        )
        # +2 because first 2 columns are the X and Z coordinates in the mesh
        array = []
        column_idx = available_parameters.index(parameter) + 2
        for fid in sorted(fids):
            array.append(np.loadtxt(fid).T[:, column_idx])

        return np.array(array, dtype="object")
