
import os
import numpy as np

from collections import Mapping
from os.path import abspath, join, exists
from string import find
from seisflows.tools import unix
from seisflows.tools.tools import Struct


class SeisStruct(Struct):
    """ Holds information about data
    """

    def __init__(self, nr=0, nt=0, dt=0., ts=0.,
                 sx=[], sy=[], sz=[],
                 rx=[], ry=[], rz=[],
                 nrec=[], nsrc=[]):
        super(SeisStruct, self).__init__(
            [['nr', nr], ['nt', nt], ['dt', dt], ['ts', ts],
             ['sx', sx], ['sy', sy], ['sz', sz],
             ['rx', rx], ['ry', ry], ['rz', rz],
             ['nrec', nrec], ['nsrc', nsrc]])


def getpar(key, file='DATA/Par_file', sep='=', cast=str):
    """ Reads parameter from SPECFEM parfile
    """
    val = None
    with open(file, 'r') as f:
        # read line by line
        for line in f:
            if find(line, key) == 0:
                # read key
                key, val = _split(line, sep)
                if not key:
                    continue
                # read val
                val, _ = _split(val, '#')
                val.strip()
                break

    if val:
        if cast == float:
            val = val.replace('d', 'e')
        return cast(val)

    else:
        print 'Not found in parameter file: %s\n' % key
        raise Exception


def setpar(key, val, file='DATA/Par_file', path='.', sep='='):
    """ Writes parameter to SPECFEM parfile
    """

    val = str(val)

    # read line by line
    with open(path + '/' + file, 'r') as f:
        lines = []
        for line in f:
            if find(line, key) == 0:
                # read key
                key, _ = _split(line, sep)
                # read comment
                _, comment = _split(line, '#')
                n = len(line) - len(key) - len(val) - len(comment) - 2
                # replace line
                if comment:
                    line = _merge(key, sep, val, ' '*n, '#', comment)
                else:
                    line = _merge(key, sep, str(val), '\n')
            lines.append(line)

    # write file
    _writelines(path + '/' + file, lines)


def Model(keys):
    return dict((key, []) for key in keys)


class Minmax(object):
    def __init__(self, keys):
        self.keys = keys
        self.minvals = dict((key, +np.Inf) for key in keys)
        self.maxvals = dict((key, -np.Inf) for key in keys)

    def items(self):
        return ((key, self.minvals[key], self.maxvals[key]) for key in self.keys)

    def update(self, keys, vals):
        for key,val in zip(keys, vals):
            minval = val.min()
            maxval = val.max()
            minval_all = self.minvals[key]
            maxval_all = self.maxvals[key]
            if minval < minval_all: self.minvals.update({key: minval})
            if maxval > maxval_all: self.maxvals.update({key: maxval})

    def write(self, path, logpath):
        if not logpath:
            return
        filename = join(logpath, 'output.minmax')
        with open(filename, 'a') as f:
            f.write(abspath(path)+'\n')
            for key,minval,maxval in self.items():
                f.write('%-15s %10.3e %10.3e\n' % (key, minval, maxval))
            f.write('\n')


class StepWriter(object):
    """ Utility for writing one or more columns to text file
    """
    def __init__(self, path='./output.optim'):
        self.iter = 0
        self.filename = abspath(path)

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


class Writer(object):
    """Utility for appending values to text files"""

    def __init__(self, path='./output.stat'):
        self.path = abspath(path)
        try:
            os.mkdir(path)
        except:
            raise IOError

        self.__call__('step_count', 0)

    def __call__(self, filename, val):
        fullfile = join(self.path, filename)
        with open(fullfile, 'a') as f:
            f.write('%e\n' % val)



### utility functions

def _split(str, sep):
    n = find(str, sep)
    if n >= 0:
        return str[:n], str[n + len(sep):]
    else:
        return str, ''


def _merge(*parts):
    return ''.join(parts)


def _writelines(file, lines):
    """ Writes text file
    """
    with open(file, 'w') as f:
        f.writelines(lines)

