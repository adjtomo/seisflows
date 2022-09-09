#! /usr/bin/env python
"""
Run an IPython notebook (execute all cells) and then convert the notebook to a 
Sphinx .rst doc page in HTML

.. rubric::
    python convert.py  # to convert all .ipynb files
    OR
    python convert.py {notebook}.ipynb  # to convert a given notebook

.. note::
    You will not need this if using the nbsphinx extension. 

Copied and edited from Pyflex documentation
"""
import os
import sys
import shutil
import glob


def convert_nb(nbname):
    """
    Open up a Jupyter notebook, execute the code inside and save to a new
    notebook, and then convert the executed notebook to a .rst file for docs.
    """
    nbname = os.path.splitext(nbname)[0]
    filename = f"{nbname}.ipynb"
    rst_filename = filename.replace("ipynb", "rst")

    # Run nbconvert to convert executed notebook to RST format for docs page
    os.system(f"jupyter nbconvert --to rst {filename}")

    # Overwrite any existing .rst file with newly created one
    rst_path = f"../{rst_filename}"
    if os.path.exists(rst_path):
        os.remove(rst_path)
    os.rename(rst_filename, rst_path)


def adjust_nb(nbname):
    """
    Since the notebooks are not stored in the same directory, their images
    will point at the wrong directory. Just make sure that we adjust the
    path and then move the created images to the correct location. Probably
    an easier way to do this but for this scale this should work.
    """
    nbname = os.path.splitext(nbname)[0]
    filename = f"{nbname}.ipynb"
    rst_filename = filename.replace("ipynb", "rst")
    rst_path = f"../{rst_filename}"

    # Move the created images files to a new directory
    nb_files_path_old = f"{nbname}_files"
    nb_files_path_new = f"../images/{nb_files_path_old}"

    if os.path.exists(nb_files_path_old):
        if os.path.exists(nb_files_path_new):
            shutil.rmtree(nb_files_path_new)
        os.rename(nb_files_path_old, nb_files_path_new)

    with open(rst_path, "r") as f:
        rst_file = f.readlines()

    # Check if any images in the converted notebook
    for i, line in enumerate(rst_file[:]):
        if ".. image::" in line:
            # If so, replace with correct image location
            new_line = line.replace(".. image:: ", ".. image:: images/")
            rst_file[i] = new_line

    with open(rst_path, "w") as f:
        f.writelines(rst_file)


if __name__ == "__main__":
    if sys.argv[1:]:
        for nbname in sys.argv[1:]:
            convert_nb(nbname)
            adjust_nb(nbname)
    else:
        for nbname in glob.glob("*ipynb"):
            convert_nb(nbname)
            adjust_nb(nbname)
