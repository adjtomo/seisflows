#!/usr/bin/python

from seisflows.tools import unix
from seisflows.tools.code import glob

unix.rm('pickle')
unix.rm('scratch')
unix.rm('nodenum')

unix.rm(glob('job*'))
unix.rm(glob('output*'))
unix.rm(glob('*.pyc'))

