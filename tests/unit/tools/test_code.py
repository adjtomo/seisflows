import unittest

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


if __name__ == '__main__':
    unittest.main()