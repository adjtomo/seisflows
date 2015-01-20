import unittest

import sys
import uuid
import imp
import os
import os.path

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
        os.remove(file_name)


class TestNull(unittest.TestCase):
    def test_init(self):
        self.assertDictEqual(tools.Null().__dict__, {})
        self.assertDictEqual(tools.Null('args').__dict__, {})

    def test_call(self):
        n = tools.Null()
        self.assertEqual(n, n())
        self.assertEqual(n, n(1))

    def test_nonzero(self):
        self.assertFalse(bool(tools.Null()))

    def test_getattr(self):
        n = tools.Null()
        self.assertEqual(n, n.__getattr__('key'))

    def test_setattr(self):
        n = tools.Null()
        self.assertEqual(n, n.__setattr__('key', 'value'))

    def test_delattr(self):
        n = tools.Null()
        self.assertEqual(n, n.__delattr__('key'))


class TestLoadClass(unittest.TestCase):
    def test_noargs(self):
        self.assertEqual(tools.Null, tools.loadclass())

    def test_load_base_module(self):
        # Get a class from one of seisflow module
        cls = tools.loadclass('system', 'serial')
        # Check if we can instanciate the class.
        self.assertIsInstance(cls(), cls),

    def test_load_extension_module(self):
        # Get a class from one of seisflow module
        cls = tools.loadclass('system', 'tiger_sm_job')
        # Check if we can instanciate the class.
        self.assertIsInstance(cls(), cls),


class TestParse(unittest.TestCase):
    def test_no_package(self):
        args = ['m' + str(uuid.uuid4().get_hex()[0:6]) for i in range(1, 10)]
        self.assertEqual(args, tools._parse(tuple(args)))

    def test_with_package(self):
        package = 'p' + str(uuid.uuid4().get_hex()[0:6])
        args = ['m' + str(uuid.uuid4().get_hex()[0:6]) for i in range(1, 10)]
        self.assertEqual([package] + args,
                         tools._parse(tuple(args), package=package))

