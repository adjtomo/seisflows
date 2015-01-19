import unittest

import os
import uuid

import seisflows.tools.code as tools


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
        not_names = [None, 0, 1, False, True]
        for name in existent:
            self.assertTrue(tools.exists(name))
        for name in non_existent:
            self.assertFalse(tools.exists(name))
        for name in not_names:
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

    def test_setdiff(self):
        list_ints = range(1, 10)
        list_odds = range(1, 10, 2)
        list_evens = range(2, 10, 2)

        self.assertEqual(set(list_odds), tools.setdiff(list_ints, list_evens))
        self.assertEqual(set(list_evens), tools.setdiff(list_ints, list_odds))
        self.assertEqual(set(list_ints), tools.setdiff(list_ints, []))

    def test_unique(self):
        self.assertEqual([], tools.unique([]))
        self.assertEqual([1,2,3,4], tools.unique([1,2,3,4]))
        self.assertEqual([1], tools.unique([1,1,1,1]))

    def test_savetxt(self):
        filename = "tmp_savetxt"
        x = 3.14159265359
        tools.savetxt(filename, x)
        with open(filename, 'r') as f:
           y = float(f.readline())
        os.remove(filename)
        self.assertAlmostEqual(x, y, places=6)


if __name__ == '__main__':
    unittest.main()