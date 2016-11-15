
from os import mkdir
from os.path import abspath, basename, dirname, exists, getsize, join
from struct import calcsize, pack, unpack

import numpy as _np

from seisflows.tools.code import Struct


class BinaryReader(object):
    """Generic binary file reader"""
    # TODO: Seems very fragile. More robust implementation?

    def __init__(self, fname, endian='|'):
        """
        For endianness values, see:
         docs.python.org/2/library/struct.html#byte-order-size-and-alignment
        """
        # opens binary file
        self.file = open(fname, 'r')
        path = abspath(self.file.name)
        self.path = dirname(path)
        self.name = basename(path)
        self.size = getsize(path)
        self.endian = endian

    def __del__(self):
        self.file.close()

    def read(self, fmt, length=1, offset=0):
        """
        For fmt values, see:
         docs.python.org/2/library/struct.html#format-characters
        """
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
            val.append(unpack(fmtlist, string)[0])

        return val

    def scan(self, fmtlist, origin=0, contiguous=1):
        """
        Take a list of formats:
            [[fmt, origin, offset, description], [...], ...]
        and returns a dictionary:
            {description: value, ...}
        """
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
        path = abspath(self.file.name)
        self.path = dirname(path)
        self.name = basename(path)
        self.size = getsize(path)
        self.endian = endian

    def __del__(self):
        self.file.close()

    def write(self, fmt, vals, length=1, offset=0):
        """Write binary data to file, all with the same format"""
        # write binary data
        if offset != 0:
            self.file.seek(offset)

        if fmt is 'bit48':
            return [[]]

        fmtlist = self.endian + mychar(fmt)

        if length == 1:
            vals = [vals]

        for val in vals:
            self.file.write(pack(fmtlist, val))

    def printf(self, fmts, vals, origin=0, contiguous=True):
        """Write binary data to file, with variable format"""
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
    return calcsize(mychar(fmt))

