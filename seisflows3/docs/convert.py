#! /usr/bin/env python
"""
Convert empty IPython notebook to a Sphinx .rst doc page in HTML

Copied and edited from Pyflex documentation
"""
import io
import os
import sys
import glob

from nbformat import read, write


def convert_nb(nbname):
    """
    Open up a Jupyter notebook, execute the code inside and save to a new
    notebook, and then convert the executed notebook to a .rst file for docs.
    """
    # Remove file extension
    nbname = os.path.splitext(nbname)[0]
    filename = f"{nbname}.ipynb"

    # Run nbconvert to execute a notebook, allowing errors through and saving
    # the new, executed notebook in place
    os.system(f"jupyter nbconvert --to notebook "
              f"--execute {filename} --inplace")
    os.system("rm -rf ./index_files")

    # Run nbconvert to convert executed notebook to RST format for docs page
    os.system(f"jupyter nbconvert --to rst {filename}")


if __name__ == "__main__":
    if sys.argv[1:]:
        for nbname in sys.argv[1:]:
            convert_nb(nbname)
    else:
        for nbname in glob.glob("*ipynb"):
            convert_nb(nbname)
