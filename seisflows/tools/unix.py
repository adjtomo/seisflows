
import os
import shutil
import socket
import subprocess
import sys

from os.path import abspath, basename, isdir, isfile, join


def cat(src, *dst):
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
    os.chdir(path)


def cp(src='', dst='', opt=''):
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
    return socket.gethostname()


def ln(src, dst):
    if os.path.isdir(dst):
        for name in _strlist(src):
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
    for dir in _strlist(dirs):
        if not os.path.isdir(dir):
            os.makedirs(dir)


def mkdir_gpfs(dirs):
    # hack to deal with race condition
    try:
        for dir in _strlist(dirs):
            if not os.path.isdir(dir):
                os.makedirs(dir)
    except:
        pass


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


def pwd():
    return os.getcwd()


def rename(old, new, names):
    for name in names:
        if name.find(old) >= 0:
            os.rename(name, name.replace(old, new))


def rm(path=''):
    for name in _strlist(path):
        if os.path.isfile(name):
            os.remove(name)
        elif os.path.islink(name):
            os.remove(name)
        elif os.path.isdir(name):
            shutil.rmtree(name)


def run(args):
    child = subprocess.Popen(args, shell=1)
    streamdata = child.communicate()[0]
    if child.returncode!=0:
        sys.exit(-1)


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


def whoami():
    import getpass

    return getpass.getuser()


def _strlist(object):
    if not isinstance(object, list):
        return [object]
    else:
        return object
