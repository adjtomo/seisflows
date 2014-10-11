
import glob as _glob
import json as _json
import os as _os
import pickle as _pickle
import sys as _sys

import numpy as _np


class Struct(dict):
    "On the fly structure"
    def __init__(self,*args,**kwargs):
        super(Struct,self).__init__(*args,**kwargs)
        self.__dict__ = self


def abspath(*args,**kwargs):
    return _os.path.abspath(*args,**kwargs)

def divides(i,j):
    "Returns true if j divides i"
    if j is 0:
        return False
    elif i%j:
        return False
    else:
        return True

def exists(name):
    "Wrapper for os.path.exists"
    if not name:
        return False
    else:
        return _os.path.exists(name)

def glob(*args):
    "Wrapper for glob.glob"
    return _glob.glob(*args)

def irange(*args):
    if len(args)==1:
        return range(1,args[0]+1)
    elif len(args)==2:
        return range(args[0],args[1]+1)

def join(*args):
    "Wrapper for os.path.join"
    return _os.path.join(*args)

def loadobj(filename):
    "Load object using pickle"
    with open(filename,'rb') as file:
        return _pickle.load(file)

def saveobj(filename,obj):
    "Save object using pickle"
    with open(filename,'wb') as file:
        _pickle.dump(obj,file)

def loadjson(filename):
    "Load object using json"
    with open(filename,'rb') as file:
        return _json.load(file)

def savejson(filename,obj):
    "Save object using json"
    with open(filename,'wb') as file:
        _json.dump(obj, file, sort_keys=True, indent=4)

def setdiff(list1,list2):
    set1 = set(list1)
    set2 = set(list2)
    diff = set1.difference(set2)
    return list(diff)

def unique(mylist):
    "Finds unique elements of list"
    return list(set(mylist))


#

def loadtxt(filename):
    return float(_np.loadtxt(filename))

def savetxt(filename,v):
    _np.savetxt(filename,[v],'%11.6e')
