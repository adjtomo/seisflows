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
        names = ['m' + str(uuid.uuid4().get_hex()[0:6]) for i in range(10)]
        values = ['v' + str(uuid.uuid4().get_hex()[0:6]) for i in range(10)]

        # Needs at least one name for initialization.
        config = tools.ConfigObj(names[0])

        # Before registering object, no keys are present.
        # self.assertEqual(config.keys, [])

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
        self.assertEqual(config.keys, set([]))

        for name in names:
            # Un-registering removes objects from sys.modules
            self.assertNotIn(name, sys.modules.keys())

    def test_saveload(self):
        # FIXME: too be completed.
        # Save and load might need some modifications. In particular:
        #     * save appends a ".p" to the name of the file being saved,
        #        while load does not.
        #     * load uses a hardcoded list of modules to look for. It matches
        #       seisflows structures. It will be troublesome if we need to
        #       change it.
        pass


class TestParameterObj(unittest.TestCase):
    def setUp(self):
        self.sys_modules_bck = sys.modules.copy()

    def tearDown(self):
        sys.modules = self.sys_modules_bck.copy()

    def test_init_non_exisiting(self):
        name = 'm' + str(uuid.uuid4().get_hex()[0:6])
        param = tools.ParameterObj(name)

        self.assertEqual(param, sys.modules[name])

    def test_init_exisiting(self):
        name = 'm' + str(uuid.uuid4().get_hex()[0:6])
        # Registering the same object twice
        param1 = tools.ParameterObj(name)

        param2 = tools.ParameterObj(name)

        self.assertEqual(param1, sys.modules[name])
        self.assertEqual(param2, sys.modules[name])
        self.assertEqual(param1, param2)

    def test_attr(self):
        name = 'm' + str(uuid.uuid4().get_hex()[0:6])
        param = tools.ParameterObj(name)

        keys = ['k' + str(uuid.uuid4().get_hex()[0:6]) for i in range(1, 10)]
        values = ['v' + str(uuid.uuid4().get_hex()[0:6]) for i in range(1, 10)]

        for k, v in zip(keys, values):
            param.__setattr__(k, v)
            self.assertEqual(param.__getattr__(k), v)

        # No attributes modification
        for k, v in zip(keys, values):
            self.assertRaises(Exception, param.__setattr__, k, v)

        # No attributes deletion
        for k, v in zip(keys, values):
            self.assertRaises(Exception, param.__delattr__, k)

    def test_update(self):
        name = 'm' + str(uuid.uuid4().get_hex()[0:6])
        param = tools.ParameterObj(name)

        keys = ['k' + str(uuid.uuid4().get_hex()[0:6]) for i in range(1, 10)]
        values = ['v' + str(uuid.uuid4().get_hex()[0:6]) for i in range(1, 10)]
        dic = dict(zip(keys, values))

        param.update(dic)

        # Checking initial dictionary update
        for (k, v) in dic.iteritems():
            self.assertEqual(v, param.__getattr__(k))

        keys = ['k' + str(uuid.uuid4().get_hex()[0:6]) for i in range(1, 10)]
        values = ['v' + str(uuid.uuid4().get_hex()[0:6]) for i in range(1, 10)]
        new_dic = dict(zip(keys, values))

        param.update(new_dic)
        # Checking second dictionary update
        for (k, v) in new_dic.iteritems():
            self.assertEqual(v, param.__getattr__(k))
        # Checking old dictionary items have been erased.
        for (k, v) in dic.iteritems():
            self.assertRaises(KeyError, param.__getattr__, k)

    def test_save(self):
        # Fill a ParamObj with a random dictionary
        name = 'm' + str(uuid.uuid4().get_hex()[0:6])
        param = tools.ParameterObj(name)

        keys = ['k' + str(uuid.uuid4().get_hex()[0:6]) for i in range(1, 10)]
        values = ['v' + str(uuid.uuid4().get_hex()[0:6]) for i in range(1, 10)]
        dic = dict(zip(keys, values))

        param.update(dic)

        # Save the object dictionary to file
        file_name = 'f' + str(uuid.uuid4().get_hex()[0:6])
        param.save(file_name)

        # Read back the file
        from seisflows.tools.code import loadjson
        read_dic = loadjson(file_name)

        self.assertEqual(dic, read_dic)