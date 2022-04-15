"""
Utilities to interact with, manipulate or call on the external solver, 
i.e., SPECFEM2D/3D/3D_GLOBE
"""
import os
import sys
import numpy as np
import subprocess

from collections import defaultdict
from seisflows3.tools import msg
from seisflows3.tools.math import poissons_ratio
from seisflows3.tools.wrappers import iterable


class Minmax(defaultdict):
    """
    Keeps track of min, max values of model or kernel
    """
    def __init__(self):
        super(Minmax, self).__init__(lambda: [+np.inf, -np.inf])

    def update(self, keys, vals):
        for key, val in _zip(keys, vals):
            if min(val) < self.dict[key][0]:
                self.dict[key][0] = min(val)
            if max(val) > self.dict[key][1]:
                self.dict[key][1] = max(val)

    def __call__(self, key):
        return self.dict[key]


class Container(defaultdict):
    """
    Dictionary-like object for holding models or kernels
    """
    def __init__(self):
        super(Container, self).__init__(lambda: [])
        self.minmax = Minmax()


def call_solver(mpiexec, executable, output="solver.log"):
    """
    Calls MPI solver executable to run solver binaries, used by individual
    processes to run the solver on system. If the external solver returns a 
    non-zero exit code (failure), this function will return a negative boolean.

    :type mpiexec: str
    :param mpiexec: call to mpi. If None (e.g., serial run, defaults to ./)
    :type executable: str
    :param executable: executable function to call
    :type output: str
    :param output: where to redirect stdout
    """
    # mpiexec is None when running in serial mode, so e.g., ./xmeshfem2D
    if mpiexec is None:
        exc_cmd = f"./{executable}"
    # Otherwise mpiexec is system dependent (e.g., srun, mpirun)
    else:
        exc_cmd = f"{mpiexec} {executable}"

    try:
        # Write solver stdout (log files) to text file
        f = open(output, "w")
        subprocess.run(exc_cmd, shell=True, check=True, stdout=f)
    except (subprocess.CalledProcessError, OSError) as e:
        print(msg.cli("The external numerical solver has returned a nonzero "
                      "exit code (failure). Consider stopping any currently "
                      "running jobs to avoid wasted computational resources. "
                      f"Check 'scratch/solver/mainsolver/{output}' for the "
                      f"solvers stdout log message. "
                      f"The failing command and error message are: ",
                      items=[f"exc: {exc_cmd}", f"err: {e}"],
                      header="external solver error",
                      border="=")
              )
        sys.exit(-1)
    finally:
        f.close()


def getpar(key, file, delim="=", match_partial=False):
    """
    Reads and returns parameters from a SPECFEM or SeisFlows3 parameter file
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


def check_poissons_ratio(vp, vs, min_val=-1., max_val=0.5):
    """
    Check Poisson's ratio based on Vp and Vs model vectors. Exit SeisFlows3 if
    Poisson's ratio is outside `min_val` or `max_val` which by default are
    set internally by SPECFEM. Otherwise return the value

    :type vp: np.array
    :param vp: P-wave velocity model vector
    :type vs: np.array
    :param vp: S-wave velocity model vector
    :type min_val: float
    :param min_val: minimum model-wide acceptable value for poissons ratio
    :type max_val: float
    :param max_val: maximum model-wide acceptable value for poissons ratio
    :return:
    """
    poissons = poissons_ratio(vp=vp, vs=vs)
    if (poissons.min() < min_val) or (poissons.max() > max_val):
        print(msg.cli(f"The Poisson's ratio of the given model is out of "
                      f"bounds with respect to the defined range "
                      f"({min_val}, {max_val}). "
                      f"The model bounds were found to be:",
                      items=["{pmin:.2f} < PR < {pmax:.2f}"], border="=",
                      header="Poisson's Ratio Error")
              )
        sys.exit(-1)
    return poissons


def _split(string, sep):
    """
    Utility function to split a string by a given separation character or str

    :type string: str
    :param string: string to split
    :type sep: str
    :param sep: substring to split by
    """
    n = string.find(sep)
    if n >= 0:
        return string[:n], string[n + len(sep):]
    else:
        return string, ''


def _merge(*parts):
    """
    Utility function to merge various strings together with no breaks
    """
    return ' '.join(parts)


def _zip(keys, vals):
    """
    Zip together keys and vals

    :type keys: dict_keys
    :param keys: keys
    :type vals: dict_values
    :param vals: values
    """
    return zip(iterable(keys), iterable(vals))
