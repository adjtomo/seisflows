"""
Unix functions wrapped in Python3. Used to simplify function calling and
provide a uniform look to SeisFlows source code. These Python functions
are meant to mimic commonly-used shell commands.
"""
import os
import time
import random
import shutil
import socket
import subprocess


def _iterable(arg):
    """
    Make an argument iterable. Allows for more generalized inputs to these
    unix-style functions.

    :type arg: anything
    :param arg: an argument to make iterable
    :rtype: list
    :return: iterable argument
    """
    if not isinstance(arg, (list, tuple)):
        return [arg]
    else:
        return arg


def cat(src, dst=None):
    """
    Concatenate files and print to standard output or write to file

    :type src: str
    :param src: file to read
    :type dst: str
    :param dst: file to write contents to
    """
    with open(src, "r") as f:
        contents = f.read()

    if dst is None:
        print(contents)
    else:
        with open(dst, "w") as f:
            f.write(contents)


def cd(path):
    """
    Change directory to `path`

    :type path: str
    :param path: path to change directory to
    """
    os.chdir(path)


def cp(src, dst):
    """
    Copy files

    :type src: str or list or tuple
    :param src: source to copy from
    :type dst: str
    :param dst: destination to copy to
    """
    if isinstance(src, (list, tuple)):
        if len(src) > 1:
            assert os.path.isdir(dst), f"unix.cp 'dst' must be directory: {dst}"
        for sub in src:
            cp(sub, dst)
        return

    if os.path.isdir(dst):
        dst = os.path.join(dst, os.path.basename(src))
        if os.path.isdir(dst):
            for sub in ls(src):
                cp(os.path.join(src, sub), dst)
            return

    if os.path.isfile(src):
        shutil.copy(src, dst)

    elif os.path.isdir(src):
        shutil.copytree(src, dst)


def hostname():
    """
    Check the system's hostname

    :rtype: str
    :return: system hostname
    """
    return socket.gethostname().split('.')[0]


def ln(src, dst):
    """
    Make a symbolic link between files

    .. rubric::
        >>> from seisflows.tools.unix import ln
        >>> ln("example_file", "path/to/sylink/new_filename")
        >>> # OR
        >>> ln("example_file", "path/to/sylink/")

    :type src: str
    :param src: path to file or directory to symlink
    :type dst: str
    :param dst: path or file to symlink `src` to
    """
    src = os.path.abspath(src)
    dst = os.path.abspath(dst)
    # If symlinking directly to a directory, we need to append file name,
    # otherwise we expect the user wants to change the filename with the symlink
    if os.path.isdir(dst):
        dst = os.path.join(dst, os.path.basename(src))
    os.symlink(src, dst)


def ls(path=os.getcwd(), show_all=False):
    """
    List directory contents, option to show hidden files

    :type path: str
    :param path: path to list directory contents
    :type show_all: bool
    :param: if True, shows hidden files (starting with '.'), same as the
        -a or --all in shell ls command
    """
    dirs = os.listdir(path)
    if not show_all:
        for dir_ in dirs:
            if dir_[0].startswith("."):
                dirs.remove(dir_)
    return dirs


def mkdir(dirs):
    """
    Make directory or directories

    .. note::
        Random wait times to prevent overloading disk when multiple subprocesses
        are running mkdir
        
    :type dirs: str or list
    :param dirs: pathnames to make
    """
    time.sleep(2 * random.random())  # interval [0, 2]s

    for dir_ in _iterable(dirs):
        if not os.path.isdir(dir_):
            os.makedirs(dir_)


def mv(src, dst):
    """
    Move contents from `src` to `dst`

    :type src: str or list or tuple
    :param src: input file(s) to move
    :type dst: str
    :param dst: location to move path to
    """
    # Src can be multiple files
    if isinstance(src, (list, tuple)):
        if len(src) > 1:
            assert os.path.isdir(dst), \
                "'dst' must be a directory for multiple input `src` in unix.mv"
        for sub in src:
            mv(sub, dst)
    else:
        if os.path.isdir(dst):
            dst = os.path.join(dst, os.path.basename(src))
        shutil.move(src, dst)


def rename(old, new, names):
    """
    Rename multiple files

    :type old: str
    :param old: expression to replace
    :type new: str
    :param new: replacement expression
    :type names: list
    :param names: files to replace expressions in
    """
    for name in _iterable(names):
        if name.find(old) >= 0:
            os.rename(name, name.replace(old, new))


def rm(path):
    """
    Remove files or directories
    """
    for name in _iterable(path):
        if os.path.isfile(name) or os.path.islink(name):
            os.remove(name)
        elif os.path.isdir(name):
            shutil.rmtree(name)


def select(items, prompt=''):
    """
    Monitor file descriptors, waiting for one or more descriptor to be "ready"
    """
    while True:
        if prompt:
            print(prompt)
        for i, item in enumerate(items):
            print(f"{i+1:2d}) {item}")
        try:
            reply = int(input().strip())
            status = (1 <= reply <= len(items))
        except (ValueError, TypeError, OverflowError):
            status = 0
        if status:
            return items[reply - 1]


def touch(filename, times=None):
    """
    Update timestamps on files

    :type filename: str
    :param filename: file to touch
    :type times: None or (atime, mtime)
    :param times: if None, set time to current time, otherwise
        (accesstime, modifiedtime) need to be set
    """
    with open(filename, "a"):
        os.utime(filename, times)


def which(name):
    """
    Shows the full path of shell commands and executables

    :type name: str
    :param name: name of shell command to check
    """
    def isexe(file):
        return os.path.isfile(file) and os.access(file, os.X_OK)

    dirname, filename = os.path.split(name)
    if dirname:
        if isexe(name):
            return name

    for path in os.environ["PATH"].split(os.pathsep):
        path = path.strip('"')
        fullname = os.path.join(path, name)
        if isexe(fullname):
            return fullname
    else:
        return None


def nproc():
    """
    Get the number of processors available. Same as calling 'nproc' from
    Linux command line.

    TODO replace all instances of nproc() with os.cpu_count()


    :rtype: int
    :return: number of processors
    :raises EnvironmentError: if nproc cannot be determined
    """
    _nproc = os.cpu_count()

    # Method 1 calls 'nproc'. May fail and return '' if 'nproc' not avail.
    if not _nproc:
        _nproc = subprocess.run("nproc", shell=True, text=True,
                                stdout=subprocess.PIPE).stdout.strip()
    # Method 2 checks /proc/cpuinfo
    if not _nproc:
        if os.path.exists("/proc/cpuinfo"):
            processors = 0
            lines = open("/proc/cpuinfo", "r").readlines()
            for line in lines:
                if line.startswith("processor"):
                    processors += 1
            if processors:
                _nproc = processors
    if not _nproc:
        raise EnvironmentError("Could not access 'nproc' information on system")

    return int(_nproc)
