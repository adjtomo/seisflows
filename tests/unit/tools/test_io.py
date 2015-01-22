import unittest

import os
import struct
from tempfile import NamedTemporaryFile

import seisflows.tools.io as tools


class TestBinaryReader(unittest.TestCase):
    def setUp(self):
        # Create a temporary binary file. Make sure it is not deleted upon
        # closing.
        self.tmp_file = NamedTemporaryFile(mode='wb', delete=False)

    def tearDown(self):
        os.remove(self.tmp_file.name)

    def test_read_native_chars(self):
        # Write some binary values
        value = 'abcdef'
        self.tmp_file.file.write(bytearray(value))
        self.tmp_file.close()

        # Read the values back
        reader = tools.BinaryReader(self.tmp_file.name, endian='=')
        r = reader.read('c', len(value), 0)
        self.assertEqual(r, list(value))

    def test_scan(self):
        # Write some binary values
        fmts = [
            ['i', 1, 0, 'Int32Bits'],
            ['h', 1, 4, 'Int16Bits'],
            ['c', 1, 6, 'Character']]
        values = (42, 33, 'a')
        endian = '='  # native
        s = struct.Struct(endian + 'i h c')
        self.tmp_file.write(s.pack(*values))
        self.tmp_file.close()

        # Read the values back
        reader = tools.BinaryReader(self.tmp_file.name, endian='=')
        r = reader.scan(fmts)
        self.assertEqual(r['Int32Bits'], 42)
        self.assertEqual(r['Int16Bits'], 33)
        self.assertEqual(r['Character'], 'a')


class TestBinaryWriter(unittest.TestCase):
    def setUp(self):
        # Create a temporary binary file. Make sure it is not deleted upon
        # closing.
        self.tmp_file = NamedTemporaryFile(mode='wb', delete=False)
        self.tmp_file.close()

    def tearDown(self):
        os.remove(self.tmp_file.name)

    def test_write_single_values(self):
        writer = tools.BinaryWriter(self.tmp_file.name, endian='=')
        writer.write('i', 42)
        writer.write('h', 33)
        writer.write('c', 'a')
        del writer

        reader = tools.BinaryReader(self.tmp_file.name, endian='=')
        self.assertEqual([42], reader.read('i', 1, 0))
        self.assertEqual([33], reader.read('h', 1, 4))
        self.assertEqual(['a'], reader.read('c', 1, 6))

    def test_write_int_array(self):
        writer = tools.BinaryWriter(self.tmp_file.name, endian='=')
        values = [2,4,6,8,10]
        writer.write('i', values, length=len(values))
        del writer

        reader = tools.BinaryReader(self.tmp_file.name, endian='=')
        self.assertEqual(values, reader.read('i', length=len(values)))


