"""
Utilities to interact with, manipulate or call on the external solver,
i.e., SPECFEM2D/3D/3D_GLOBE
"""
import os
import numpy as np
import shutil
from glob import glob
from obspy.geodetics import gps2dist_azimuth

from seisflows import logger
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


def get_src_rcv_lookup_table(path_to_data, source_prefix="CMTSOLUTION",
                             stations_file="STATIONS"):
    """
    Generate a lookup table that gives relative distance, azimuth etc. for
    each source and receiver used in the workflow. Source and station
    locations will be gathered from metadata files stored in the SPECFEM
    data directory.

    :type path_to_data: str
    :param path_to_data: full path to the SPECFEM DATA/ directory which should
        contain the source file (e.g., CMTSOLUTION, FORCESOLUTION), and the
        STATIONS file whose name is defined by `station_file`.
    :type source_prefix: str
    :param source_prefix: prefix of all the source files that will be wildcard
        searched for using the following wildcard: {source_prefix}_*. e.g., if
        all your files are CMTSOLUTIONS (CMTSOLUTION_001, CMTSOLUTION_002), then
        `source_prefix` should be 'CMTSOLUTION'
    :type stations_file: str
    :param stations_file: the name of the STATIONS file in the SPECFEM DATA
        directory `path_to_data`. By default this is STATIONS
    """
    # Determine source locations 'latitude' and 'longitude'
    src_dict = get_source_locations(path_to_sources=path_to_data,
                                    source_prefix=source_prefix)

    # Station dictionary of locations 'latitude', 'longitude'
    rcv_dict = get_station_locations(
        stations_file=os.path.join(path_to_data, stations_file)
    )
    dict_out = Dict()
    for src_name, src_vals in src_dict.items():
        dict_out[src_name] = Dict()
        for rcv_name, rcv_vals in rcv_dict.items():
            # Calculate the azimuth from North of the source (theta) and azimuth
            # from North of the receiver (theta prime). See Fig. 1 from
            # Wang et al. (2019) for diagrammatic explanation.
            dist_m, az, baz = gps2dist_azimuth(lat1=src_vals["latitude"],
                                               lon1=src_vals["longitude"],
                                               lat2=rcv_vals["latitude"],
                                               lon2=rcv_vals["longitude"]
                                               )
            # Theta is the azimuth from north of the source, and theta prime
            # is the azimuth from north of the receiver. Theta != Theta' for
            # a spherical Earth, but they will be close.
            theta = np.deg2rad(az)
            theta_p = np.deg2rad((baz - 180) % 360)

            dict_out[src_name][rcv_name] = Dict(dist_m=dist_m, az=az,
                                                baz=baz, theta=theta,
                                                theta_p=theta_p)

    return dict_out


def return_matching_waveform_files(obs_path, syn_path, obs_fmt="ASCII",
                                   syn_fmt="ASCII", components=None):
    """
    Gather a list of filenames of matching waveform IDs which are sorted and
    match for given paths to 'observed' and 'synthetic' waveform data. This is
    useful because often the available 'observed' data will not match all of
    the 'synthetic' data that is generated during a simulation, so this function
    helps determine which files/stations can be accessed for misfit
    quantification.

    .. note::

        Waveform files are expected to be in the format
        NN.SSS.CCc.* (N=network, S=station, C=channel, c=component;
        following SPECFEM ASCII formatting). They will be matched on
        `NN.SSS.c` (dropping channel naming because SEED convention may have
        different channel naming). For example, synthetic name
        'AA.S001.BXZ.semd' will be converted to 'AA.S001.Z', and matching
        observation 'AA.S001.HHZ.SAC' will be converted to 'AA.S001.Z'.
        These two will be matched.

    :type obs_path: str
    :param obs_path: the path to observed data waveforms for a given source.
        In SeisFlows this will typically point to somewhere like:
         'scratch/solver/<SOURCE_NAME>/traces/obs/'
    :type syn_path: str
    :param syn_path: the path to synthetic data waveforms for a given source.
        In SeisFlows this will typically point to somewhere like:
         'scratch/solver/<SOURCE_NAME>/traces/syn/'
    :type obs_fmt: str
    :param obs_fmt: expected file format of the observed waveforms. Used for
        safety checks that only one expected file format will be read. Defaults
        to 'ASCII'
    :type syn_fmt: str
    :param syn_fmt: expected file format of the synthetic waveforms. Used for
        safety checks that only one expected file format will be read. Defaults
        to 'ASCII'
    :type components: list
    :param components: optional list of components to ignore preprocessing
        traces that do not have matching components. The adjoint sources for
        these components will be 0. E.g., ['Z', 'N']. If None, all available
        components will be considered.
    :rtype: list of tuples
    :return: [(observed filename, synthetic filename)]. tuples will contain
        filenames for matching stations + component for obs and syn
    """
    observed = sorted(os.listdir(obs_path))
    synthetic = sorted(os.listdir(syn_path))

    logger.debug(f"found {len(observed)} obs and {len(synthetic)} syn "
                 f"waveforms for given paths")

    assert(len(observed) != 0 and len(synthetic) != 0), \
        f"cannot quantify misfit, missing observed or synthetic traces"

    # Verify all traces have acceptable formats that are expected by SeisFlows
    for key, files, fmt in zip(["observed", "synthetic"],
                               [observed, synthetic],
                               [obs_fmt, syn_fmt]):
        # Check all available file extensions to ensure that there is only one
        exts = list(set([os.path.splitext(x)[-1] for x in files]))
        assert(len(exts) == 1), (
            f"{key} files have too many available file formats ({exts}), and "
            f"only 1 was expected ({fmt})"
        )
        ext = exts[0].upper()  # e.g., .ASCII

        # Check if the expected file format matches the provided files
        if fmt.upper() == "ASCII":
            # SPECFEM3D and 3D_GLOBE name their ASCII files slightly differently
            ext_ok = (ext == ".ASCII") or (".SEM" in ext)
        elif fmt.upper() == "SU":
            raise NotImplementedError
        else:
            ext_ok = bool(ext == f".{fmt}")
        assert ext_ok, f"{key} data unexpected file format {ext} != {fmt}"

    # Format path/to/NN.SSS.CCc* -> NN.SSS.c (see docstring note for details)
    match_obs, match_syn = [], []
    for files, lists in zip([observed, synthetic], [match_obs, match_syn]):
        # Drop full path incase these are given as absolute paths
        full_fids = [os.path.basename(_) for _ in files]
        # Split into expected format NN.SSS.CCc, drop extension
        fids = [_.split(".")[:3] for _ in full_fids]
        for fid in fids:
            net, sta, cha = fid
            comp = cha[-1]
            # Ignore non-selected components if User requests
            if components and comp not in components:
                lists.append("")
                continue
            # NN.SSS.c
            lists.append(f"{net}.{sta}.{comp}")

    # only return traces that have both observed and synthetic file match
    match_traces = \
            sorted(list(set(match_syn).intersection(set(match_obs))))
    logger.info(f"{len(match_traces)} traces matching between `obs` and `syn`")

    # Final check to makes sure that we still have data that can be compared
    assert(len(match_traces) != 0), (
        f"there are no traces with both observed and synthetic files for "
        f"the given paths. Please verify that waveform files within paths "
        f"have the format 'NN.SSS.CCc*', which match on variables 'N', 'S', "
        f"and 'c'"
    )

    # Generate the list of full path waveform fids for matching obs + syn
    obs_paths, syn_paths = [], []
    for short_fid in match_traces:
        # Find the corresponding full path name based on the short match vrs
        obs_fid = observed[match_obs.index(short_fid)]
        obs_paths.append(os.path.join(obs_path, obs_fid))

        syn_fid = synthetic[match_syn.index(short_fid)]
        syn_paths.append(os.path.join(syn_path, syn_fid))

    return obs_paths, syn_paths


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


def rename_as_adjoint_source(fid, fmt):
    """
    Rename SPECFEM synthetic waveform filenames consistent with how SPECFEM
    expects adjoint sources to be named. Usually this just means adding
    a '.adj' to the end of the filename.

    :type fid: str
    :param fid: file path of synthetic waveform to rename as adjoint source
    :type fmt: str
    :param fmt: expected format of the input synthetic waveform, because
        different file formats have different filename structure. Available
        are 'SU' (seismic unix) and 'ASCII'. Case-insensitive
    :rtype: str
    :return: renamed file that matches expected SPECFEM filename format
        for adjoint sources
    """
    # Safety check to make sure were not trying to rename adjoint source files
    if fid.endswith(".adj"):
        logger.warning(f"file to be renamed already ends with .adj: {fid}")
        return fid

    if fmt.upper() == "SU":
        fid = f"{fid}.adj"
    elif fmt.upper() == "ASCII":
        # Differentiate between SPECFEM3D and 3D_GLOBE file naming
        # SPECFEM3D: NN.SSSS.CCC.sem?
        # SPECFEM3D_GLOBE: NN.SSSS.CCC.sem.ascii
        ext = os.path.splitext(fid)[-1]
        # SPECFEM3D
        if ".sem" in ext:
            fid = fid.replace(ext, ".adj")
        # GLOBE (!!! Hardcoded to only work with ASCII format)
        elif ext == ".ascii":
            root, ext1 = os.path.splitext(fid)  # .ascii
            root, ext2 = os.path.splitext(root)  # .sem
            fid = fid.replace(f"{ext2}{ext1}", ".adj")

    return fid


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


def getpar(key, file, delim="=", match_partial=False, comment="#",
           _fmt_dbl=True):
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
    :type comment: str
    :param comment: character used to delimit comments in the file. Defaults to
        '#' for the SPECFEM Par_file, but things like the FORCESOLUTION use '!'
    :type _fmt_dbl: bool
    :param _fmt_dbl: the SPECFEM files use FORTRAN double precision notation
        to define floats (e.g., 2.5d0 == 2.5, or 1.2e2 == 1.2*10^2). Usually it
        is preferable to convert this notation directly to a Python float, so
        this is set True by default. However, in cases where we are doing string
        replacement, like when using `setpar`, we do not want to format the
        double precision values, so this should be set False.
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
                # Split on the first occurrence, assuming key has not delim
                key_out, val = line.strip().split(delim, 1)
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
            if comment in val:
                val = val.split(comment)[0]
            # One last strip to remove any whitespace
            val = val.strip()
            # Address the fact that SPECFEM Par_file sometimes lists values as
            # formatted strings, e.g., 38.0d-2
            try:
                if len(val.split("d")) == 2 and _fmt_dbl:
                    num, exp = val.split("d")
                    val = str(float(num) * 10 ** int(exp))
            except ValueError as e:
                # This will break on anything other than than the above type str
                pass
            break
    else:
        raise KeyError(f"Could not find matching key '{key}' in file: {file}")
    return key_out, val, i


def setpar(key, val, file, delim="=", **kwargs):
    """
    Overwrites parameter value to a SPECFEM Par_file. Kwargs passed to `getpar`

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
    key_out, val_out, i = getpar(key, file, delim, _fmt_dbl=False, **kwargs)
    if key_out is None:
        return

    lines = open(file, "r").readlines()

    # Replace value in place by splitting on the delimiter and replacing the
    # first instance of the old value (avoids replacing comment values or
    # matching values inside the key)
    if val_out != "":
        key, val_and_comment = lines[i].split(delim)
        val_and_comment = val_and_comment.replace(val_out, str(val), 1)
        assert(str(val) in val_and_comment), (
                f"unexpected error in parameter value replacement when trying "
                f"to replace '{val_out}' with '{val}' in '{val_and_comment}'"
                )
        lines[i] = delim.join([key, val_and_comment])
    else:
        # Special case where the initial parameter is empty so we just replace
        # the newline formatter at the end
        lines[i] = lines[i].replace("\n", f" {val}\n")
    with open(file, "w") as f:
        f.writelines(lines)


def _getidx_vel_model(lines):
    """
    Get the line indices of a velocity model, which can be used to retrieve 
    or replace the model values in a SPECFEM2D paramter file. Used by 
    `getpar_vel_model` and `setpar_vel_model`

    :type lines: list
    :param lines: list of strings read from the par_file
    :rtype idxs: list
    :param idxs: list of integer indices of the velocity model lines
    """
    idxs = []
    for l, line in enumerate(lines):
        # Skip over all other parameters, comments, newlines etc.
        if "=" in line:
            continue
        elif line.startswith(" "):
            continue
        elif line.startswith("#"):
            continue
        elif line.startswith("\n"):
            continue

        # Matches formatting expected by velocity model and starts with integer
        # Should be enough to avoid matching other strings
        lines = line.strip().split()
        if len(lines) == 15 and lines[0].isdigit():
            idxs.append(l)

    return idxs


def getpar_vel_model(file, strip=False):
    """
    SPECFEM2D doesn't follow standard key = val formatting when defining its 
    internal velocity models so we need a special function to address this 
    specifically.
    
    Velocity models are ASSUMED to be formatted in the following way

    1 1 2700.d0 3000.d0 1732.051d0 0 0 9999 9999 0 0 0 0 0 0

    That is, 15 entries separated by spaces. We use that to find all relevant 
    lines of the model.
    
    :type file: str
    :param file: The SPECFEM Par_file to match against
    :type strip: bool
    :param strip: strip newline '\n' from each of the model lines
    :rtype: list of str
    :return: list of all the layers of the velocity model as strings
    """
    # i_start + 1 to skip over the 'nbmodels' parameter
    lines = open(file, "r").readlines()
    idxs = _getidx_vel_model(lines)
    vel_model = []
    for idx in idxs:
        if strip:
            line = lines[idx].strip()
        else:
            line = lines[idx]
        vel_model.append(line)

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
    lines = open(file, "r").readlines()
    model_lines = _getidx_vel_model(open(file, "r").readlines())

    # one-liner to get rid of the original model from the file
    lines = [i for j, i in enumerate(lines) if j not in model_lines]
    model_idx_start = model_lines[0]

    # Throw a new line onto the last line of the model to get proper formatting
    model[-1] = model[-1] + "\n"

    # Drop in the new model one line at a time
    for i, val in enumerate(model):
        lines.insert(model_idx_start + i, f"{val.strip()}\n")

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
