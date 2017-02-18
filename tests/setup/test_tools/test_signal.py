
import unittest

import numpy as np

from seisflows.tools import signal


class TestSeistoolsSignal(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_mask(self):
        # time scheme
        nt = 100
        dt = 1.e-3

        t1 = dt
        t2 = t1 + (nt-1)*dt

        # mute parameters
        slope = 0.     # units: time/distance
        const = t2/2.  # units: time
        offset = 0.    # units: distance

        w = signal.mask(slope, const, offset, (nt, dt, t1), length=10)

        # check monotonicity
        assert(np.all(np.diff(w) >= 0))


if __name__ == '__main__':
    unittest.main()


