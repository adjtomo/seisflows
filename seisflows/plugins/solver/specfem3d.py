#!/usr/bin/env python3
"""
Plugins used for the numerical solver SPECFEM3D_CARTESIAN
"""
import os
from seisflows.tools.wrappers import findpath
from seisflows.tools.specfem import setpar


def write_sources(par, h, path="."):
    """
    Writes FORCESOLUTION source information to text file

    :type par: dict
    :param par: seisflows.PAR
    :type h:
    :param h:
    :type path: str
    :param path: path to write sources to
    """
    file = os.path.join(findpath("seisflows.plugins"), "specfem3d",
                        "FORCESOLUTION")
    with open(file, "r") as f:
        lines = f.readlines()

    file = "DATA/FORCESOURCE"
    with open(file, "w") as f:
        f.writelines(lines)

    # adjust coordinates
    setpar("xs", h.sx[0], file)
    setpar("zs", h.sz[0], file)
    setpar("ts", h.ts, file)

    # adjust wavelet
    setpar("f0", par["F0"], file)


def write_receivers(h):
    """
    Writes receiver information to text file

    :type h:
    :param h:
    """
    file = "DATA/STATIONS"
    lines = []

    # loop over receivers
    for ir in range(h.nr):
        line = ""
        line += "S%06d" % ir + " "
        line += "AA" + " "
        line += "%11.5e" % h.rx[ir] + " "
        line += "%11.5e" % h.rz[ir] + " "
        line += "%3.1f" % 0. + " "
        line += "%3.1f" % 0. + "\n"
        lines.extend(line)

    with open(file, "w") as f:
        f.writelines(lines)


