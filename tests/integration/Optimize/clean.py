#!/bin/env python

import glob

from seisflows.tools import unix

unix.rm('scratch')

unix.rm(glob('output*'))
unix.rm(glob('*.pyc'))

