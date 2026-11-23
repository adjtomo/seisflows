#!/usr/bin/env python3
"""
Find and remove the large files of a SeisFlows working directory that can be
made again by running more simulations, and keep the rest: the parameter file,
state file, logs, optimization vectors, residuals and the directory structure
(including the run directory of every event) that a workflow needs to carry on.

Files are put into categories so that a User can see where the space is, and
choose what to remove. A file is only ever in one category, whichever of the
categories are chosen, so that leaving one out does not remove its files under
the name of another. Only regular files below the scratch directory are ever
removed, never a directory or a link, or anything a link points to.
"""
import os
from glob import glob
from datetime import datetime

from seisflows.tools.specfem import getpar

# Name: (removed unless left out, what it is). In the order that a file is
# given a category, so a forward array is not a 'database'
CATEGORIES = {
    "forward_arrays": (True, "forward wavefields and absorbing boundaries "
                             "saved for the adjoint simulations. Usually "
                             "the largest"),
    "vtk": (True, "VTK files of the databases of each run directory"),
    "databases": (True, "model and mesh databases in the run directories "
                        "(DATABASES_MPI), made by the mesher or from a model"),
    "kernels": (True, "event kernels, gradient and misfit kernel in "
                      "eval_grad"),
    "models": (True, "copies of the models in eval_func and eval_grad, the "
                     "optimization vectors that they come from are kept"),
    "traces": (True, "synthetic waveforms and adjoint sources in the run "
                     "directories"),
    "obs": (False, "observed waveforms in the run directories"),
    "solver_logs": (False, "timestamps and other small output of the "
                           "solver"),
}

# Where SPECFEM puts its databases, relative to a run directory. 'DATA' is
# where SPECFEM2D reads models from. The 'LOCAL_PATH' of the Par_file is added
DATABASE_DIRS = ["DATA", "OUTPUT_FILES",
                 os.path.join("OUTPUT_FILES", "DATABASES_MPI")]

# File name patterns in the database directories, in order of category
DATABASE_PATTERNS = {
    "forward_arrays": ["proc??????_save_forward_arrays*",
                       "proc??????_absorb_*", "proc??????_reg?_absorb_*",
                       "absorb_*_*"],
    "vtk": ["proc??????_*.vt?"],
    "databases": ["proc??????_*"],
}

# File name patterns relative to a run directory
RUN_PATTERNS = {
    "traces": [os.path.join("traces", "syn", "*"),
               os.path.join("traces", "adj", "*"),
               os.path.join("SEM", "*"),
               os.path.join("OUTPUT_FILES", "*.sem?"),
               os.path.join("OUTPUT_FILES", "*.su")],
    "obs": [os.path.join("traces", "obs", "*")],
    "solver_logs": [os.path.join("OUTPUT_FILES", "timestamp*"),
                    os.path.join("OUTPUT_FILES", "output_list_s*.txt"),
                    os.path.join("OUTPUT_FILES", "sr.vtk"),
                    os.path.join("OUTPUT_FILES", "starttimeloop.txt")],
}


def select_categories(include=None, exclude=None):
    """
    Choose the categories to remove: the ones that are removed unless they are
    left out, and the ones that are added, less the ones that are left out.

    :type include: list of str
    :param include: categories to remove as well as the usual ones
    :type exclude: list of str
    :param exclude: categories to leave alone
    :rtype: list of str
    :return: names of the categories to remove
    :raises ValueError: if a name is not a category
    """
    include, exclude = list(include or []), list(exclude or [])
    for name in include + exclude:
        if name not in CATEGORIES:
            raise ValueError(f"'{name}' is not a category, choose from: "
                             f"{', '.join(CATEGORIES)}")
    return [name for name, (usual, _) in CATEGORIES.items()
            if (usual or name in include) and name not in exclude]


def find_files(scratch, solver=None, eval_func=None, eval_grad=None,
               categories=None):
    """
    The files in each category, found from the directories of a working
    directory. Only regular files that are really below `scratch` are found.
    Nothing is removed.

    :type scratch: str
    :param scratch: scratch directory, which nothing outside of is touched
    :type solver: str
    :param solver: directory of the run directories of the events, defaults
        to `scratch/solver`
    :type eval_func: str
    :param eval_func: line search evaluation directory, defaults to
        `scratch/eval_func`
    :type eval_grad: str
    :param eval_grad: gradient evaluation directory, defaults to
        `scratch/eval_grad`
    :type categories: list of str
    :param categories: the categories to return, defaults to the usual ones.
        All categories are looked for to give each file its category, whatever
        is chosen
    :rtype: dict
    :return: category: sorted list of absolute paths of the files in it
    """
    solver = solver or os.path.join(scratch, "solver")
    eval_func = eval_func or os.path.join(scratch, "eval_func")
    eval_grad = eval_grad or os.path.join(scratch, "eval_grad")
    if categories is None:
        categories = select_categories()

    root = os.path.realpath(scratch)
    found = {name: [] for name in CATEGORIES}
    seen = set()

    def add(category, fid):
        """A regular file that is really inside of scratch, not yet seen"""
        parent = os.path.realpath(os.path.dirname(fid))
        if (fid in seen or os.path.islink(fid) or not os.path.isfile(fid) or
                os.path.commonpath([root, parent]) != root):
            return
        seen.add(fid)
        found[category].append(fid)

    def add_tree(category, path):
        """Every regular file below a directory, which is not a link"""
        if os.path.islink(path):
            return
        for dirpath, _, fids in os.walk(path, followlinks=False):
            for fid in fids:
                add(category, os.path.join(dirpath, fid))

    # The run directory of each event, not the link to the first of them
    run_dirs = sorted(path for path in glob(os.path.join(solver, "*"))
                      if os.path.isdir(path) and not os.path.islink(path))
    for run_dir in run_dirs:
        database_dirs = list(DATABASE_DIRS)
        try:
            local_path = os.path.normpath(getpar(
                "LOCAL_PATH", os.path.join(run_dir, "DATA", "Par_file"))[1])
            if not os.path.isabs(local_path) and local_path not in \
                    database_dirs:
                database_dirs.append(local_path)
        except (Exception, SystemExit):
            pass  # not a SPECFEM3D Par_file

        # The forward arrays of noise inversion, one directory for each force
        for path in glob(os.path.join(run_dir, "*_FWD_ARR")):
            add_tree("forward_arrays", path)
        for category, patterns in DATABASE_PATTERNS.items():
            for database_dir in database_dirs:
                for pattern in patterns:
                    for fid in glob(os.path.join(run_dir, database_dir,
                                                 pattern)):
                        add(category, fid)
        for category, patterns in RUN_PATTERNS.items():
            for pattern in patterns:
                for fid in glob(os.path.join(run_dir, pattern)):
                    add(category, fid)

    for name in ["kernels", "gradient", "misfit_kernel"]:
        add_tree("kernels", os.path.join(eval_grad, name))
    for path in [eval_func, eval_grad]:
        add_tree("models", os.path.join(path, "model"))

    return {name: sorted(found[name]) for name in categories}


def size_of(fids):
    """
    Total size of files

    :type fids: list of str
    :param fids: paths of files
    :rtype: int
    :return: size in bytes
    """
    return sum(os.path.getsize(fid) for fid in fids)


def human_size(nbytes):
    """
    A size for a person to read

    :type nbytes: int
    :param nbytes: size in bytes
    :rtype: str
    :return: e.g., '1.5 GB'
    """
    for unit in ["B", "kB", "MB", "GB"]:
        if nbytes < 1000 or unit == "GB":
            return f"{nbytes:.0f} {unit}" if unit == "B" else \
                f"{nbytes:.1f} {unit}"
        nbytes /= 1000


def format_files(found, workdir, list_files=False):
    """
    Say what is to be removed, for a person to read

    :type found: dict
    :param found: category: files, from `find_files`
    :type workdir: str
    :param workdir: paths are given relative to this
    :type list_files: bool
    :param list_files: if True, name every file and its size, otherwise only
        say how many there are of each category
    :rtype: str
    :return: the text
    """
    lines = []
    for name, fids in found.items():
        lines.append(f"{name}: {len(fids)} files, "
                     f"{human_size(size_of(fids))}  ({CATEGORIES[name][1]})")
        if list_files:
            for fid in fids:
                lines.append(f"    {os.path.relpath(fid, workdir)}  "
                             f"{human_size(os.path.getsize(fid))}")
    nfiles = sum(len(fids) for fids in found.values())
    total = human_size(sum(size_of(fids) for fids in found.values()))
    lines.append(f"total: {nfiles} files, {total}")
    return "\n".join(lines)


def remove_files(found, workdir, record=None):
    """
    Remove the files, and keep a record of what was removed

    :type found: dict
    :param found: category: files, from `find_files`
    :type workdir: str
    :param workdir: paths in the record are given relative to this
    :type record: str
    :param record: path of a text file that a line for each file removed is
        added to, with the category and size
    :rtype: tuple
    :return: (number of files removed, their size in bytes, the files that
        could not be removed)
    """
    nremoved, nbytes, failed = 0, 0, []
    lines = [f"# cleaned {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"]
    for name, fids in found.items():
        for fid in fids:
            try:
                size = os.path.getsize(fid)
                os.remove(fid)
            except OSError:
                failed.append(fid)
                continue
            nremoved += 1
            nbytes += size
            lines.append(f"{name}\t{size}\t{os.path.relpath(fid, workdir)}\n")

    if record and nremoved:
        os.makedirs(os.path.dirname(record), exist_ok=True)
        with open(record, "a") as f:
            f.writelines(lines)

    return nremoved, nbytes, failed


def tree_summary(path):
    """
    How much there is in a file or directory, for saying what is to be removed
    when everything is

    :type path: str
    :param path: path to a file or directory
    :rtype: tuple
    :return: (number of files, size in bytes), links are counted as the file
        that is removed, which is the link
    """
    if os.path.islink(path) or os.path.isfile(path):
        return 1, os.lstat(path).st_size
    nfiles, nbytes = 0, 0
    for dirpath, _, fids in os.walk(path, followlinks=False):
        for fid in fids:
            nfiles += 1
            nbytes += os.lstat(os.path.join(dirpath, fid)).st_size
    return nfiles, nbytes
