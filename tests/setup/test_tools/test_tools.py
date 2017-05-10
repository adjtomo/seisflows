import unittest

import os
import uuid
from tempfile import NamedTemporaryFile

from seisflows.tools import tools


class TestToolsCode(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_divides(self):
        divisible = [(1, 1), (2, 1), (100, 1), (50, 2), (12, 3)]
        not_divisible = [(1, 0), (5, 0), (5, 4), (13, 2), (3, 6)]
        for (n, m) in divisible:
            self.assertTrue(tools.divides(n, m))
        for (n, m) in not_divisible:
            self.assertFalse(tools.divides(n, m))

    def test_exists(self):
        existent = os.listdir(os.curdir)
        non_existent = [str(uuid.uuid4().get_hex().upper()[0:6])
                        for i in range(10)]
        for name in existent:
            self.assertTrue(tools.exists(name))
        for name in non_existent:
            self.assertFalse(tools.exists(name))

    def test_saveload(self):
        filename = "tmp_pickle_file"
        obj_init = "Something"
        tools.saveobj(filename, obj_init)
        obj_read = tools.loadobj(filename)
        self.assertEqual(obj_read, obj_init)
        os.remove(filename)

    def test_saveload_json(self):
        filename = "tmp_json_file"
        obj_init =[u'foo', {u'bar': [u'baz', None, 1.0, 2]}]
        tools.savejson(filename, obj_init)
        obj_read = tools.loadjson(filename)
        self.assertEqual(obj_read, obj_init)
        os.remove(filename)

    def test_savetxt(self):
        filename = "tmp_savetxt"
        x = 3.14159265359
        tools.savetxt(filename, x)
        with open(filename, 'r') as f:
           y = float(f.readline())
        os.remove(filename)
        self.assertAlmostEqual(x, y, places=6)

    def test_saveloadtxt(self):
        tmp_file = NamedTemporaryFile(mode='wb', delete=False)
        tmp_file.close()

        x = 3.14159265359
        tools.savetxt(tmp_file.name, x)
        self.assertAlmostEqual(x, tools.loadtxt(tmp_file.name), 6)


if __name__ == '__main__':
    unittest.main()
