
import os as _os
import struct as _struct

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
        path = _os.path.abspath(self.file.name)
        self.path = _os.path.dirname(path)
        self.name = _os.path.basename(path)
        self.size = _os.path.getsize(path)
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
            val.append(_struct.unpack(fmtlist, string)[0])

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
        path = _os.path.abspath(self.file.name)
        self.path = _os.path.dirname(path)
        self.name = _os.path.basename(path)
        self.size = _os.path.getsize(path)
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
            self.file.write(_struct.pack(fmtlist, val))

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


class OutputWriter(object):
    """ Utility for writing one or more columns to text file
    """
    def __init__(self, filename, width=10):
        try:
            self.filename = _os.path.abspath(filename)
        except:
            raise IOError

        # column width
        self.width = width


    def __call__(self, *vals):
        fileobj = open(self.filename, 'a')
        line = ''
        for val in vals:
            line += self._getline(val)
        fileobj.write(line + '\n')
        fileobj.close()

    def _getline(self, val):
        w = self.width

        if type(val) is int:
            fmt = '%%%dd  ' % w
        elif type(val) is float:
            fmt = '%%%d.%de  ' % (w,w-7)
        elif type(val) is str:
            fmt = '%%%ds  ' % w
        else:
            return (w+2)*' '

        return fmt % val

    def newline(self):
        # writes newline
        fileobj = open(self.filename, 'a')
        fileobj.write('\n')
        fileobj.close()

    def header(self, *keys):
        # writes column headers
        line = ''
        for key in keys:
            fmt = '%%%ds  ' % self.width
            line += fmt % key.upper()

        if _os.path.exists(self.filename):
            fileobj = open(self.filename, 'a')
        else:
            fileobj = open(self.filename, 'w')
            fileobj.write(line + '\n')
            fileobj.write((len(keys)*((self.width*'=') + '  ')) + '\n')
        fileobj.close()



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

