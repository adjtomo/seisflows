#!/bin/env python

import glob

from seisflows.tools import unix

unix.rm('pickle')
unix.rm('scratch')
unix.rm('nodenum')

unix.rm(glob('job*'))
unix.rm(glob('output*'))
unix.rm(glob('*.pyc'))

