"""
Utilities to interact with, manipulate or call on the external solver, 
i.e., SPECFEM2D/3D/3D_GLOBE
"""
import os
import numpy as np
from glob import glob
from seisflows import logger
from seisflows.tools.core import Dict
from seisflows.tools import unix
from seisflows.tools.math import poissons_ratio


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
    acceptable_parameters.extend([f"{_}_kernel" for _ in acceptable_parameters])

    def __init__(self, path, fmt=None, parameters=None, load=False):
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
        :type load: bool
        :param load: load the model into disk as dictionary and vector
            representations. If False, will guess how to import data based on
            file extensions or lack thereof
        """
        assert os.path.exists(path), f"specfem model path {path} does not exist"
        self.path = path
        self.parameters = parameters

        # Load an existing model
        if os.path.splitext(path)[-1] == ".npz" or load:
            self.model, self.ngll, self.fmt = self.load(file=self.path)
            _first_key = list(self.model.keys())[0]
            self.nproc = len(self.model[_first_key])
        # Read a model from files
        else:
            if fmt is None:
                self.fmt = self._guess_file_format()
            else:
                self.fmt = fmt
            self.nproc, self.available_parameters = self._get_nproc_parameters()
            self.model, self.ngll = self.read(parameters=parameters)

        # .sorted() enforces parameter order every time, otherwise things can
        # get screwy if keys returns different each time
        self.parameters = sorted(self.model.keys())

    @staticmethod
    def fnfmt(i="*", val="*", ext="*"):
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
    def vector(self):
        """conveience property to access the merge() property which creates a
        linear vector defining all model parameters"""
        return self.merge()

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
            assert(set(parameters).union(set(self.available_parameters))), (
                f"user-chosen parameters not in available: "
                f"{self.available_parameters}"
            )

        # Pick the correct read function based on the file format
        load_fx = {".bin": self._read_model_fortran_binary,
                   ".dat": self._read_model_ascii,
                   ".adios": self._read_model_adios   # TODO Check if this okay
                   }[self.fmt]

        # Create a dictionary object containing all parameters and their models
        parameter_dict = Dict({key: [] for key in parameters})
        for parameter in parameters:
            parameter_dict[parameter] = load_fx(parameter=parameter)

        # Save some metadata to be able to manipulate model slices freely
        ngll = []
        for array in parameter_dict[parameters[0]]:
            ngll.append(len(array))

        return parameter_dict, ngll

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
            assert(self.fmt is not None), f"must specifiy model format: `fmt`"
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

    def check(self, min_pr=-1., max_pr=0.5):
        """
        Checks parameters in the model. If Vs and Vp present, checks poissons
        ratio. Checks for negative velocity values. And prints out model
        min/max values
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

        if "vs" in self.model and self.model.vs.min() < 0:
            logger.warning(f"Vs minimum is negative {self.model.vs.min()}")

        if "vp" in self.model and self.model.vp.min() < 0:
            logger.warning(f"Vp minimum is negative {self.model.vp.min()}")

        # Tell the User min and max values of the updated model
        for key, vals in self.model.items():
            # Choose formatter based on the magnitude of the value
            if vals.min() < 1 or (vals.max() > 1E4):
                parts = "{minval:.2E} <= {key} <= {maxval:.2E}"
            else:
                parts = "{minval:.2f} <= {key} <= {maxval:.2f}"
            logger.info(parts.format(minval=vals.min(), key=key,
                                     maxval=vals.max()))

    def save(self, path):
        """
        Save instance attributes (model, vector, metadata) to disk as an
        .npz array so that it can be loaded in at a later time for future use
        """
        model = self.split()
        np.savez(file=path, fmt=self.fmt, **model)

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
        ngll = []
        data = np.load(file=file)
        for i, key in enumerate(data.files):
            if key == "fmt":
                continue
            model[key] = data[key]
            if not ngll:
                for array in model[key]:
                    ngll.append(len(array))

        return model, ngll, str(data["fmt"])

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

        return self

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
            nproc = len(glob(os.path.join(
                self.path, self.fnfmt(val=avail_par[0], ext=self.fmt)))
                        )
        elif self.fmt == ".dat":
            # e.g., 'proc000000_rho_vp_vs'
            _, *avail_par = fids[0].split("_")
            nproc = len(fids) + 1
        else:
            raise NotImplementedError(f"{self.fmt} is not yet supported by "
                                      f"SeisFlows")

        # Remove any parameters not accepted by Model
        avail_par = set(avail_par).intersection(set(self.acceptable_parameters))

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

    def _read_model_fortran_binary(self, parameter):
        """
        Load Fortran binary models into disk. This is the preferred model format
        for SeisFlows <-> SPECFEM interaction

        :type parameter: str
        :param parameter: chosen parameter to load model for
        :rtype: np.array
        :return: vector of model values for given `parameter`
        """
        def _read(filename):
            """Read a single slice (e.g., proc000000_vs.bin) binary file"""
            nbytes = os.path.getsize(filename)
            with open(filename, 'rb') as file:
                # read size of record
                file.seek(0)
                n = np.fromfile(file, dtype='int32', count=1)[0]

                if n == nbytes-8:
                    file.seek(4)
                    data = np.fromfile(file, dtype='float32')
                    return data[:-1]
                else:
                    file.seek(0)
                    data = np.fromfile(file, dtype='float32')
                    return data

        array = []
        fids = glob(os.path.join(
            self.path, self.fnfmt(val=parameter, ext=".bin"))
        )
        for fid in sorted(fids):  # make sure were going in numerical order
            array.append(_read(fid))

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
        all in a single file with all parameters together as a N column ASCII file
        where columns 1 and 2 are the coordinates of the mesh, and the remainder
        columns are data corresponding to the filenames
        e.g., proc000000_rho_vp_vs.dat, rho is column 3, vp is 4 etc.

        :type parameter: str
        :param parameter: chosen parameter to load model for
        :rtype: np.array
        :return: vector of model values for given `parameter`
        """
        fids = glob(os.path.join(self.path, self.fnfmt(val="*", ext=".dat")))
        _, *available_parameters = fids[0].split("_")
        assert(parameter in available_parameters), (
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
                buffer = np.array([4 * len(data)], dtype="int32")
                data = data.astype("float32")

                with open(filepath, 'wb') as f:
                    buffer.tofile(f)
                    data.tofile(f)
                    buffer.tofile(f)


def getpar(key, file, delim="=", match_partial=False):
    """
    Reads and returns parameters from a SPECFEM or SeisFlows parameter file
    Assumes the parameter file is formatted in the following way:

    # comment comment comment
    {key}      {delim} VAL

    :type key: str
    :param key: case-insensitive key to match in par_file. must be EXACT match
    :type file: str
    :param file: The SPECFEM Par_file to match against
    :type delim: str
    :param delim: delimiter between parameters and values within the file.
        default is '=', which matches for SPECFEM2D and SPECFEM3D_Cartesian
    :type match_partial: bool
    :param match_partial: allow partial key matches, e.g., allow key='tit' to
        return value for 'title'. Defaults to False as this can have
        unintended consequences
    :rtype: tuple (str, str, int)
    :return: a tuple of the key, value and line number (indexed from 0).
        The key will match exactly how it looks in the Par_file
        The value will be returned as a string, regardless of its expected type
        IF no matches found, returns (None, None, None)
    """
    lines = open(file, "r").readlines()

    for i, line in enumerate(lines):
        # Find the first occurence, CASE-INSENSITIVE search, strip whitespace
        # To allow for nested parameters
        if line.strip().upper().startswith(key.upper()):
            try:
                key_out, val = line.strip().split(delim)
            except ValueError as e:
                raise ValueError(
                    f"Could not split line with delimiter '{delim}'. Error "
                    f"message: {e}"
                )
            # Strip white space now so we can match keys
            key_out = key_out.strip()
            # Unless explcitely authorized, don't allow partial matches,
            # which str.find() will return on. Can have unintended consequences
            # e.g., 'tit' will return 'title'
            if (not match_partial) and (key_out.upper() != key.upper()):
                continue
            # Drop any trailing line comments
            if "#" in val:
                val = val.split("#")[0]
            # One last strip to remove any whitespace
            val = val.strip()
            # Address the fact that SPECFEM Par_file sometimes lists values as
            # formatted strings, e.g., 38.0d-2
            try:
                if len(val.split("d")) == 2:
                    num, exp = val.split("d")
                    val = str(float(num) * 10 ** int(exp))
            except ValueError as e:
                # This will break on anything other than than the above type str
                pass
            break
    else:
        raise KeyError(f"Could not find matching key '{key}' in file: {file}")
    return key_out, val, i


def setpar(key, val, file, delim="=", match_partial=False):
    """
    Overwrites parameter value to a SPECFEM Par_file.

    :type key: str
    :param key: case-insensitive key to match in par_file. must be EXACT match
    :type val: str
    :param val: value to OVERWRITE to the given key
    :type file: str
    :param file: The SPECFEM Par_file to match against
    :type delim: str
    :param delim: delimiter between parameters and values within the file.
        default is '=', which matches for SPECFEM2D and SPECFEM3D_Cartesian
    :type match_partial: bool
    :param match_partial: allow partial key matches, e.g., allow key='tit' to
        return value for 'title'. Defaults to False as this can have
        unintended consequences
    """
    key_out, val_out, i = getpar(key, file, delim, match_partial)
    if key_out is None:
        return

    lines = open(file, "r").readlines()

    # Replace value in place
    if val_out != "":
        lines[i] = lines[i].replace(val_out, str(val))
    else:
        # Special case where the initial parameter is empty so we just replace
        # the newline formatter at the end
        lines[i] = lines[i].replace("\n", f" {val}\n")
    with open(file, "w") as f:
        f.writelines(lines)


def getpar_vel_model(file):
    """
    SPECFEM2D doesn't follow standard formatting when defining its internal
    velocity models so we need a special function to address this specifically.
    Velocity models are ASSUMED to be formatted in the following way in the
    SPECFEM2D Par_file (with any number of comment lines in between)

    nbmodels                        = 4
    1 1 2700.d0 3000.d0 1732.051d0 0 0 9999 9999 0 0 0 0 0 0
    2 1 2500.d0 2700.d0 0 0 0 9999 9999 0 0 0 0 0 0
    3 1 2200.d0 2500.d0 1443.375d0 0 0 9999 9999 0 0 0 0 0 0
    4 1 2200.d0 2200.d0 1343.375d0 0 0 9999 9999 0 0 0 0 0 0
    TOMOGRAPHY_FILE                 = ./DATA/tomo_file.xyz

    :type file: str
    :param file: The SPECFEM Par_file to match against
    :rtype: list of str
    :return: list of all the layers of the velocity model as strings
    """
    _, _, i_start = getpar("nbmodels", file)
    _, _, i_end = getpar("tomography_file", file)

    # i_start + 1 to skip over the 'nbmodels' parameter
    lines = open(file, "r").readlines()[i_start + 1:i_end]
    vel_model = []
    for line in lines:
        # Skip comments, empty lines, newlines
        for not_allowed in [" ", "#", "\n"]:
            if line.startswith(not_allowed):
                break
        else:
            vel_model.append(line.strip())
    return vel_model


def setpar_vel_model(file, model):
    """
    Set velocity model values in a SPECFEM2D Par_file, see getpar_vel_model
    for more information.

    Deletes the old model from the Par_file, writes the new model in the same
    place, and then changes the value of 'nbmodels'

    :type file: str
    :param file: The SPECFEM Par_file to match against
    :type model: list of str
    :param model: input model
    :rtype: list of str
    :return: list of all the layers of the velocity model as strings, e.g.:
        model = ["1 1 2700.d0 3000.d0 1732.051d0 0 0 9999 9999 0 0 0 0 0 0",
                 "2 1 2500.d0 2700.d0 0 0 0 9999 9999 0 0 0 0 0 0"]
    """
    _, nbmodels, i_start = getpar("nbmodels", file)
    i_start += 1  # increase by one to start AFTER nbmodels line
    _, _, i_end = getpar("tomography_file", file)

    lines = open(file, "r").readlines()
    model_lines = []
    # i_start + 1 to skip over the 'nbmodels' parameter
    for i, line in enumerate(lines[i_start:i_end]):
        # Skip comments, empty lines, newlines
        for not_allowed in [" ", "#", "\n"]:
            if line.startswith(not_allowed):
                break
        else:
            # We will use these indices to delete the original model
            model_lines.append(i)

    # Make sure that our indices are relative to the list and not enumeration
    model_lines = [_ + i_start for _ in model_lines]

    # one-liner to get rid of the original model
    lines = [i for j, i in enumerate(lines) if j not in model_lines]

    # Throw a new line onto the last line of the model to get proper formatting
    model[-1] = model[-1] + "\n"

    # Drop in the new model one line at a time
    for i, val in enumerate(model):
        lines.insert(i + model_lines[0], f"{val.strip()}\n")

    # Overwrite with new lines containing updated velocity model
    with open(file, "w") as f:
        f.writelines(lines)

    # Set nbmodels to the correct value
    setpar(key="nbmodels", val=len(model), file=file)


def _read(filename):
    """
    Legacy code: Reads Fortran style binary data into numpy array.

    .. note::
        Has been rewritten into the Model class but left here if useful
    """
    nbytes = os.path.getsize(filename)
    with open(filename, 'rb') as file:
        # read size of record
        file.seek(0)
        n = np.fromfile(file, dtype='int32', count=1)[0]

        if n == nbytes-8:
            file.seek(4)
            data = np.fromfile(file, dtype='float32')
            return data[:-1]
        else:
            file.seek(0)
            data = np.fromfile(file, dtype='float32')
            return data


def _write(v, filename):
    """
    Legacy code: Writes Fortran style binary files
    Data are written as single precision floating point numbers

    .. note::
        Has been rewritten into the Model class but left here if useful

    .. note::
        FORTRAN unformatted binaries are bounded by an INT*4 byte count. This
        function mimics that behavior by tacking on the boundary data.
        https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnc4/index.html
    """
    n = np.array([4 * len(v)], dtype='int32')
    v = np.array(v, dtype='float32')

    with open(filename, 'wb') as file:
        n.tofile(file)
        v.tofile(file)
        n.tofile(file)
