#
# This is Seisflows
#
# See LICENCE file
#
###############################################################################

# This file contains a bunch of useful functions to manipulate files and
# directories and more generally to mimic unix bash commands

# Import system modules
import os
import random
import shutil
import socket
import subprocess
import sys
import time

# Local imports
from os.path import abspath, basename, isdir, isfile, join
from seisflows.tools.tools import iterable


def cat(src, *dst):
    """ Open a file and print it. If file dst is given the file is printed into
        dst
    """
    f = open(src, 'r')
    contents = f.read()
    f.close()

    if not dst:
        print contents
    else:
        f = open(dst, 'w')
        f.write(contents)
        f.close()


def cd(path):
    """ Change current directory to the one given
    """
    os.chdir(path)


def cp(src='', dst=''):
    """ Copy all (files or directories, as list or tupples) given in src to
        directory (or file) dst 
    """
    if isinstance(src, (list, tuple)):
        if len(src) > 1:
            assert isdir(dst)

        for sub in src:
            cp(sub, dst)
        return

    if isdir(dst):
        dst = join(dst, basename(src))
        if isdir(dst):
            for sub in ls(src):
                cp(join(src, sub), dst)
            return

    if isfile(src):
        shutil.copy(src, dst)

    elif isdir(src):
        shutil.copytree(src, dst)


def hostname():
    return socket.gethostname().split('.')[0]


def ln(src, dst):
    dst = abspath(dst)
    if os.path.isdir(dst):
        for name in iterable(src):
            s = abspath(name)
            d = join(dst, basename(name))
            os.symlink(s, d)
    else:
        os.symlink(src, dst)


def ls(path):
    dirs = os.listdir(path)
    for dir in dirs:
        if dir[0] == '.':
            dirs.remove(dir)
    return dirs


def mkdir(dirs):
    """ Make directories or list of directories
        In parallel environments, random delays are useful prior to prevent
        filesystem overload from multiple simultaneous commands.
        This method is used in large-scale environments like amazon or
        google web services
    """
    time.sleep(1.0 * random.random())
    for dir in iterable(dirs):
        if not os.path.isdir(dir):
            os.makedirs(dir)


def mv(src='', dst=''):
    if isinstance(src, (list, tuple)):
        if len(src) > 1:
            assert isdir(dst)
        for sub in src:
            mv(sub, dst)
        return

    if isdir(dst):
        dst = join(dst, basename(src))

    shutil.move(src, dst)


def rename(old, new, names):
    for name in iterable(names):
        if name.find(old) >= 0:
            os.rename(name, name.replace(old, new))


def rm(path=''):
    for name in iterable(path):
        if os.path.isfile(name):
            os.remove(name)
        elif os.path.islink(name):
            os.remove(name)
        elif os.path.isdir(name):
            shutil.rmtree(name)


def select(items, prompt=''):
    while True:
        if prompt:
            print prompt
        for i, item in enumerate(items):
            print("%2d) %s" % (i + 1, item))
        try:
            reply = int(raw_input().strip())
            status = (1 <= reply <= len(items))
        except (ValueError, TypeError, OverflowError):
            status = 0
        if status:
            return items[reply - 1]


def touch(filename, times=None):
    with open(filename, 'a'):
        os.utime(filename, times)


def which(name):
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
