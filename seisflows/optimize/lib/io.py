
from os.path import abspath, exists, join

import os
import numpy as np

class StepWriter(object):
    """ Utility for writing one or more columns to text file
    """
    def __init__(self, path='.', filename='output.optim'):
        self.iter = 0
        self.filename = self.fullfile(path, filename)

        self.write_header()

    def __call__(self, steplen=None, funcval=None):
        with open(self.filename, 'a') as fileobj:
            if self.iter == 0:
                self.iter += 1
                fmt = '%10d  %10.3e  %10.3e\n'
                fileobj.write(fmt % (self.iter, steplen, funcval))
            elif steplen == 0.:
                self.iter += 1
                fmt = '%10d  %10.3e  %10.3e\n'
                fileobj.write(fmt % (self.iter, steplen, funcval))
            else:
                fmt = 12*' ' + '%10.3e  %10.3e\n'
                fileobj.write(fmt % (steplen, funcval))

    def write_header(self):
        # write header
        headers = []
        headers += ['ITER']
        headers += ['STEPLEN']
        headers += ['MISFIT']

        with open(self.filename, 'a') as fileobj:
            for header in headers:
                fmt = '%%%ds  ' % 10
                fileobj.write('%10s  ' % header)
            fileobj.write('\n')
            for _ in range(len(headers)):
                fileobj.write('%10s  ' % (10*'='))
            fileobj.write('\n')

    def newline(self):
        with open(self.filename, 'a') as fileobj:
                fileobj.write('\n')


    def fullfile(self, path, filename):
        try:
            fullpath = abspath(path)
        except:
            raise IOError
        return join(fullpath, filename)


class Writer(object):
    """Utility for appending values to text files"""

    def __init__(self, path='.'):
        if not exists(path):
            raise Exception
        self.path = join(path, 'NonlinearOptimization')

        os.mkdir(self.path)

        self.__call__('step_count', 0)

    def __call__(self, filename, val):
        fullfile = join(self.path, filename)
        with open(fullfile, 'a') as f:
            f.write('%e\n' % val)

