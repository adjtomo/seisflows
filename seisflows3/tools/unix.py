"""
Unix functions wrapped in Python
"""
import os
import random
import shutil
import socket
import time

from seisflows3.tools.tools import iterable


def cat(src, *dst):
    """
    Concatenate files and print on standard output
    """
    with open(src, 'r') as f:
        contents = f.read()

    if not dst:
        print(contents)
    else:
        with open(dst, 'w') as f:
            f.write(contents)


def cd(path):
    """
    Change directory
    """
    os.chdir(path)


def cp(src='', dst=''):
    """
    Copy files

    :type src: str or list or tuple
    :param src: source to copy from
    :type dst: str
    :param dst: destination to copy to
    """
    if isinstance(src, (list, tuple)):
        if len(src) > 1:
            assert os.path.isdir(dst), "unexpected type for unix.cp 'dst'"
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
    Check the hostname
    """
    return socket.gethostname().split('.')[0]


def ln(src, dst):
    """
    Make a symbolic link
    """
    dst = os.path.abspath(dst)
    if os.path.isdir(dst):
        for name in iterable(src):
            s = os.path.abspath(name)
            d = os.path.join(dst, os.path.basename(name))
            os.symlink(s, d)
    else:
        os.symlink(src, dst)


def ls(path):
    """
    List directory contents
    """
    dirs = os.listdir(path)
    for dir_ in dirs:
        if dir_[0] == '.':
            dirs.remove(dir_)
    return dirs


def mkdir(dirs):
    """
    Make directory
    Note: Random wait times to prevent overloading disk
        
    :type dirs: str or list
    :param dirs: pathnames to make
    """
    time.sleep(2 * random.random())
    for dir_ in iterable(dirs):
        if not os.path.isdir(dir_):
            os.makedirs(dir_)


def mv(src='', dst=''):
    """
    Move contents
    """
    if isinstance(src, (list, tuple)):
        if len(src) > 1:
            assert os.path.isdir(dst), "unexpected type for 'dst' in unix.mv"
        for sub in src:
            mv(sub, dst)
        return

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
    for name in iterable(names):
        if name.find(old) >= 0:
            os.rename(name, name.replace(old, new))


def rm(path=''):
    """
    Remove files or directories
    """
    for name in iterable(path):
        if os.path.isfile(name):
            os.remove(name)
        elif os.path.islink(name):
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
    with open(filename, 'a'):
        os.utime(filename, times)


def which(name):
    """
    Shows the full path of shell commands

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

