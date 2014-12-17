
import os as _os
import struct as _struct

from collections import Mapping

import numpy as _np

from seisflows.tools.code import Struct


class BinaryReader(object):
    """Generic binary file reader"""

    def __init__(self, fname, endian='|'):
        # opens binary file
        self.file = open(fname, 'r')
        path = _os.path.abspath(self.file.name)
        self.path = _os.path.dirname(path)
        self.name = _os.path.basename(path)
        self.size = _os.path.getsize(path)
        self.endian = endian

    def __del__(self):
        self.file.close()

    def read(self, fmt, length=1, offset=0):
        # reads binary file
        if offset != 0:
            self.file.seek(offset)

        if fmt is 'bit48':
            return [[]]

        val = []
        fmtlist = self.endian + mychar(fmt)
        size = mysize(fmt)

        for _ in range(length):
            string = self.file.read(size)
            val.append(_struct.unpack(fmtlist, string)[0])

        return val

    def scan(self, fmtlist, origin=0, contiguous=1):
        # reads binary file
        self.file.seek(origin)

        position = 0
        h = Struct()

        for item in fmtlist:
            fmt = item[0]
            length = item[1]
            offset = item[2]
            name = item[3]

            if not contiguous:
                self.file.seek(offset - position, 1)
                position = offset + mysize(fmt)

            # if length is 1:
            #   h[name] = self.read(fmt,length)[0]
            # else:
            #   h[name] = self.read(fmt,length)
            h[name] = self.read(fmt, length)[0]

        return h


class BinaryWriter(object):
    """Generic binary file writer"""

    def __init__(self, fname, endian='|'):
        # open binary file
        self.file = open(fname, 'w')
        path = _os.path.abspath(self.file.name)
        self.path = _os.path.dirname(path)
        self.name = _os.path.basename(path)
        self.size = _os.path.getsize(path)
        self.endian = endian

    def __del__(self):
        self.file.close()

    def write(self, fmt, vals, length=1, offset=0):
        # write binary data
        if offset != 0:
            self.file.seek(offset)

        if fmt is 'bit48':
            return [[]]

        fmtlist = self.endian + mychar(fmt)

        if length == 1:
            vals = [vals]

        for val in vals:
            self.file.write(_struct.pack(fmtlist, val))

    def printf(self, fmts, vals, origin=0, contiguous=1):
        # write binary data
        self.file.seek(origin)
        position = 0

        for i, item in enumerate(fmts):
            fmt = item[0]
            length = item[1]
            offset = item[2]
            val = vals[i]

            if not contiguous:
                self.file.seek(offset - position, 1)
                position = offset + mysize(fmt)

            self.write(fmt, val, length)


class OutputWriter(object):
    def __init__(self, filename, keys):

        try:
            self.filename = _os.path.abspath(filename)
        except:
            raise IOError

        self.nkey = len(keys)
        self.keys = []
        for key in keys:
            self.keys += [key.upper()]

        # format column headers
        line = ''
        for key in self.keys:
            line += '%10s  ' % key

        # write headers to file
        if _os.path.exists(filename):
            fileobj = open(filename, 'a')
        else:
            fileobj = open(filename, 'w')
            fileobj.write(line + '\n')
            fileobj.write((self.nkey*((10*'=') + '  ')) + '\n')
        fileobj.close()

    def __call__(self, *vals):
        fileobj = open(self.filename, 'a')
        nval = len(vals)
        if nval != self.nkey:
            raise Exception
        line = ''
        for val in vals:
            line += self._getline(val)
        fileobj.write(line + '\n')
        fileobj.close()

    def _getline(self, val):
        if val == '':
            return 12*' '
        if not val:
            return 12*' '
        if type(val) is int:
            return '%10d  ' % val
        if type(val) is float:
            return '%10.3e  ' % val
        if type(val) is str:
            return '%10s  ' % val


def loadbin(filename):
    """Reads Fortran style binary data"""
    with open(filename, 'rb') as file:
        # read size of record
        file.seek(0)
        n = _np.fromfile(file, dtype='int32', count=1)[0]

        # read contents of record
        file.seek(4)
        v = _np.fromfile(file, dtype='float32')

    return v[:-1]


def savebin(v, filename):
    """Writes Fortran style binary data"""
    n = _np.array([4*len(v)], dtype='int32')
    v = _np.array(v, dtype='float32')

    with open(filename, 'wb') as file:
        n.tofile(file)
        v.tofile(file)
        n.tofile(file)


def mychar(fmt):
    chars = {'int8': 'b',
             'uint8': 'B',
             'int16': 'h',
             'uint16': 'H',
             'int32': 'i',
             'uint32': 'I',
             'int64': 'q',
             'uint64': 'Q',
             'float': 'f',
             'float32': 'f',
             'double': 'd',
             'char': 's'}
    if fmt in chars:
        return chars[fmt]
    else:
        return fmt

def mysize(fmt):
    return _struct.calcsize(mychar(fmt))

