#!/usr/bin/env python

from glob import glob

from seisflows.tools import unix

unix.rm('niter')
unix.rm('scratch')
unix.rm(glob('output*'))
unix.rm(glob('*.pyc'))

