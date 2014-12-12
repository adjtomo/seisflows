#!/usr/bin/python

from parameters import *

from shared import *


def main(dir):
    # assert len(kernels) is NSRC
    assert exists(join(dir, 'kernels'))
    assert exists('collect_kernels.exe')

    cmd = ['collect_kernels.exe']
    cmd.extend(str(NSRC))
    cmd.extend(os.path.join(path, dir))

    subprocess.call(cmd)

    return v
