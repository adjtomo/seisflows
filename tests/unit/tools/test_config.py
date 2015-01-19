import unittest

import sys
import uuid

import seisflows.tools.config as tools


class TestConfigObj(unittest.TestCase):
    def setUp(self):
        self.sys_modules_bck = sys.modules.copy()

    def tearDown(self):
        sys.modules = self.sys_modules_bck.copy()

    def test_init(self):
        name = str(uuid.uuid4())
        co = tools.ConfigObj(name)
        self.assertTrue(name in sys.modules)

    def test_registration(self):
        names = [str(uuid.uuid4().get_hex().upper()[0:6]) for i in range(10)]
        values = [str(uuid.uuid4().get_hex().upper()[0:6]) for i in range(10)]
        config = tools.ConfigObj()

        # Before registering object, no keys are present.
        self.assertEqual(config.keys, [])

        for name, value in zip(names, values):
            config.register(name, value)

        # All registered objects should be present.
        self.assertItemsEqual(names, config.keys)

        for k in config:
            # Get the correct value for the present key.
            idx = names.index(k)
            value = values[idx]
            # Make sure that the correct value is associated to the correct key.
            self.assertEqual(value, sys.modules[k])

        for name in names:
            config.unregister(name)

        # Un-registering all objects, yields an empty set.
        self.assertEqual(config.keys, [])

        for name in names:
            # Un-registering removes objects from sys.modules
            self.assertNotIn(name, sys.modules.keys())
