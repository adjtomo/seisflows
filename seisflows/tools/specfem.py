"""
Utilities to interact with, manipulate or call on the external solver, 
i.e., SPECFEM2D/3D/3D_GLOBE
"""
import os
import numpy as np
import shutil
from glob import glob
from seisflows.tools import msg
from seisflows.tools.config import Dict

# Used for accessing loc. information from SPECFEM2D/3D/3D_GLOBE source files
SOURCE_KEYS = {
    "SOURCE": {"lat": "xs", "lon": "zs", "delim": "=", },
    "FORCESOLUTION_3D": {"lat": "latorUTM", "lon": "longorUTM", "delim": ":"},
    "CMTSOLUTION_3D": {"lat": "latorUTM", "lon": "longorUTM", "delim": ":"},
    "FORCESOLUTION_3DGLOBE": {"lat": "latitude", "lon": "longitude",
                              "delim": ":"},
    "CMTSOLUTION_3DGLOBE": {"lat": "latitude", "lon": "longitude",
                            "delim": ":"},
}


def convert_stations_to_sources(stations_file, source_file, source_type,
                                output_dir="./"):
    """
    Used for ambient noise adjoint tomography inversions where each station
    is treated like a virtual source. This requires generating source files
    for each station in a station file.

    .. warning::

        In SPECFEM3D_GLOBE, the FORCESOLUTION first line requires a specific 
        format. I haven't tested this thoroughly but it either cannot exceed
        a char limit, or cannot have characters like spaces or hyphens. To be 
        safe, something like FORCE{N} will work, where N should be the source
        number. This function enforces this just incase.

    :type stations_file: str
    :param stations_file: full path to SPECFEM STATIONS file which should be
        formatted 'STATION NETWORK LATITUDE LONGITUDE ELEVATION BURIAL',
        elevant and burial will not be used
    :type source_file: str
    :param source_file: path to 
    :type source_type: str
    :param source_type: tells SeisFlows what type of file we are using, which
        in turn defines the specific keys and delimiters to use when editing 
        the source file. 

        - SOURCE: SPECFEM2D source file
        - FORCESOLUTION_3D: SPECFEM3D Cartesian FORCESOLUTION file 
        - FORCESOLUTION_3DGLOBE: SPECFEM3D_GLOBE FORCESOLUTION file
    """
    assert(source_type in SOURCE_KEYS.keys()), \
        f"`source_type` must be in {SOURCE_KEYS.keys()}"
    keys = SOURCE_KEYS[source_type]
    _fid = source_type.split("_")[0]  # e.g., FORCESOLUTION

    stations = np.loadtxt(stations_file, dtype="str")
    for i, sta in enumerate(stations):
        station, network, latitude, longitude, *_ = sta

        # Copy the original source file to a new file
        new_source = os.path.join(output_dir, f"{_fid}_{network}_{station}")
        shutil.copy(source_file, new_source)

        # Set the new location based on the station location
        setpar(key=keys["lat"], val=latitude, file=new_source,
               delim=keys["delim"])
        setpar(key=keys["lon"], val=longitude, file=new_source,
               delim=keys["delim"])

        # Enforce first line comment, see warning in docstring for explanation
        lines = open(new_source, "r").readlines()
        with open(new_source, "w") as f:
            f.write(f"FORCE{i:0>3}\n")
            f.writelines(lines[1:])


def get_station_locations(stations_file):
    """
    Read the SPECFEM STATIONS file to get metadata information about stations.
    This functionality is required in a few preprocessing or utility functions

    :type stations_file: str
    :param stations_file: full path to STATIONS file that will be read. We
        assume the structure of the STATIONS file to be a text file with the
        first 4 columns defining 'station ID, network ID, latitude, longitude'
    :rtype: Dict
    :return: keys are `network_station` and values are dictionaries contianing
        'lat' and 'lon' for latitude and longitude values
    """
    sta_dict = {}
    stations = np.loadtxt(stations_file, dtype=str)
    for i, sta in enumerate(stations):
        station, network, latitude, longitude, *_ = sta
        sta_dict[f"{network}_{station}"] = {"latitude": float(latitude),
                                            "longitude": float(longitude)
                                            }
    return Dict(sta_dict)


def get_source_locations(path_to_sources, source_prefix):
    """
    Read in SOURCE/CMTSOLUTION/FORCESOLUTION files to get lat/lon/depth or x/y/z
    values which can be used to determine relative locations between sources
    and receivers. Used by the preprocessing module.


    :type path_to_sources: str
    :param path_to_sources: full path to all source files which should start
        with `source_prefix`. Will run a wildcard glob search on this path
        to look for ALL available source files
    :type source_prefix: str
    :param source_prefix: source file name prefix used for wildcard searching.
        This function will look for glob(`path_to_sources`/`source_prefix`_*),
        and tag each of the sources with whatever tag comes within wildcard *
    :rtype: Dict
    :return: keys are `source tag` and values are dictionaries contianing
        'lat' and 'lon' for latitude and longitude values. If None is returned,
        then function could not determine what keys to use for reading source.
    """
    src_dict = {}
    src_files = glob(os.path.join(path_to_sources, f"{source_prefix}_*"))

    # Guess the source key based on the source prefix since it's difficult
    # to get the exact solver information into this function without a lot
    # of extra bookkeeping
    for source_key in SOURCE_KEYS.keys():
        if source_prefix in source_key:
            try:
                keys = SOURCE_KEYS[source_key]
                getpar(key=keys["lat"], file=src_files[0], delim=keys["delim"])
                source_key_actual = source_key
                break
            except KeyError:
                continue
    else:
        return None

    for src_file in glob(os.path.join(path_to_sources, f"{source_prefix}_*")):
        src_name = os.path.basename(src_file).split(f"{source_prefix}_")[-1]
        key = SOURCE_KEYS[source_key_actual]
        latitude = getpar(key=key["lat"], file=src_file, delim=key["delim"])[1]
        longitude = getpar(key=key["lon"], file=src_file, delim=key["delim"])[1]
        src_dict[src_name] = {"latitude": float(latitude),
                              "longitude": float(longitude)
                              }

    return Dict(src_dict)


def check_source_names(path_specfem_data, source_prefix, ntask=None):
    """
    Determines names of sources by applying wildcard rule to user-supplied
    input files. Source names are only provided up to PAR.NTASK and are
    returned in alphabetical order.

    .. note::
        SeisFlows expects sources to be stored in the DATA/ directory with a
        prefix and a source name, e.g., {source_prefix}_{source_name} which
        would evaluate to something like CMTSOLUTION_001

    :type path_specfem_data: str
    :param path_specfem_data: path to a
    :type source_prefix: str
    :param source_prefix: type of SPECFEM input source, e.g., CMTSOLUTION
    :type ntask: int
    :parma ntask: if provided, curtails the list of sources up to `ntask`. If
        None, returns all files found matching the wildcard
    :rtype: list
    :return: alphabetically ordered list of source names up to PAR.NTASK
    """
    wildcard = f"{source_prefix}_*"
    fids = sorted(glob(os.path.join(path_specfem_data, wildcard)))
    if not fids:
        print(msg.cli("No matching source files when searching PATH for the "
                      "given WILDCARD", header="source error", border="=",
                      items=[f"PATH: {path_specfem_data}",
                             f"WILDCARD: {wildcard}"]))
        return
    if ntask is not None:
        assert(len(fids) >= ntask), (
            f"Number of requested tasks/events {ntask} exceeds number "
            f"of available sources {len(fids)}"
        )
        fids = fids[:ntask]

    # Create internal definition of sources names by stripping prefix
    names = ["_".join(os.path.basename(fid).split("_")[1:]) for fid in fids]

    return names


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

    # Replace value in place by splitting on the delimiter and replacing the 
    # first instance of the old value (avoids replacing comment values or 
    # matching values inside the key)
    if val_out != "":
        key, val_and_comment = lines[i].split(delim)
        val_and_comment = val_and_comment.replace(val_out, str(val), 1)
        lines[i] = delim.join([key, val_and_comment])
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


def read_fortran_binary(filename):
    """
    Reads Fortran-style unformatted binary data into numpy array.

    .. note::
        The FORTRAN runtime system embeds the record boundaries in the data by
        inserting an INTEGER*4 byte count at the beginning and end of each
        unformatted sequential record during an unformatted sequential WRITE.
        see: https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnc4/index.html

    :type filename: str
    :param filename: full path to the Fortran unformatted binary file to read
    :rtype: np.array
    :return: numpy array with data with data read in as type Float32
    """
    nbytes = os.path.getsize(filename)
    with open(filename, "rb") as file:
        # read size of record
        file.seek(0)
        n = np.fromfile(file, dtype="int32", count=1)[0]

        if n == nbytes - 8:
            file.seek(4)
            data = np.fromfile(file, dtype="float32")
            return data[:-1]
        else:
            file.seek(0)
            data = np.fromfile(file, dtype="float32")
            return data


def write_fortran_binary(arr, filename):
    """
    Writes Fortran style binary files. Data are written as single precision
    floating point numbers.

    .. note::
        FORTRAN unformatted binaries are bounded by an INT*4 byte count. This
        function mimics that behavior by tacking on the boundary data.
        https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnc4/index.html

    :type arr: np.array
    :param arr: data array to write as Fortran binary
    :type filename: str
    :param filename: full path to file that should be written in format
        unformatted Fortran binary
    """
    buffer = np.array([4 * len(arr)], dtype="int32")
    data = np.array(arr, dtype="float32")

    with open(filename, "wb") as file:
        buffer.tofile(file)
        data.tofile(file)
        buffer.tofile(file)
