#!/usr/bin/python

from seisflows.tools import unix
from seisflows.tools.code import glob

unix.rm('scratch')

unix.rm(glob('output*'))
unix.rm(glob('*.pyc'))

