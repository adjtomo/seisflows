import json
import os
import pickle
import sys

import numpy as np


class Struct(dict):
    def __init__(self, *args, **kwargs):
        super(Struct, self).__init__(*args, **kwargs)
        self.__dict__ = self


def divides(i, j):
    """Returns true if j divides i"""
    if j is 0:
        return False
    elif i % j:
        return False
    else:
        return True


def exists(name):
    """Wrapper for os.path.exists"""
    if not isinstance(name, basestring):
        return False
    else:
        return os.path.exists(name)


def loadobj(filename):
    """Load object using pickle"""
    with open(filename, 'rb') as file:
        return pickle.load(file)


def saveobj(filename, obj):
    """Save object using pickle"""
    with open(filename, 'wb') as file:
        pickle.dump(obj, file)


def loadjson(filename):
    """Load object using json"""
    with open(filename, 'rb') as file:
        return json.load(file)


def savejson(filename, obj):
    """Save object using json"""
    with open(filename, 'wb') as file:
        json.dump(obj, file, sort_keys=True, indent=4)


def setdiff(list1, list2):
    """Returns the difference of two list in a set.
    :param list1:
    :param list2:
    :return:
    """
    set1 = set(list1)
    set2 = set(list2)
    return set1.difference(set2)


def unique(mylist):
    """Finds unique elements of list"""
    return list(set(mylist))


def loadtxt(filename):
    """Load scalar from text file"""
    return float(np.loadtxt(filename))


def savetxt(filename, v):
    """Save scalar to text file"""
    np.savetxt(filename, [v], '%11.6e')
