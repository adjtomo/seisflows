#!/usr/bin/env python

from glob import glob
import os
import shutil

try:
    shutil.rmtree('scratch')
except:
    pass

print glob('output*')

[os.remove(f) for f in glob('output*')]
[os.remove(f) for f in glob('*.pyc')]

