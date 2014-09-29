
import os as _os
import struct as _struct


class Reader(object):

  def __init__(self,fname,endian='|'):
    # opens binary file
    self.file = open(fname,'r')
    path = _os.path.abspath(self.file.name)
    self.path = _os.path.dirname(path)
    self.name = _os.path.basename(path)
    self.size = _os.path.getsize(path)
    self.endian = endian

  def __del__(self):
    self.file.close()

  def read(self,fmt,length=1,offset=0):
    # reads binary file
    if offset != 0:
      self.file.seek(offset)

    if fmt is 'bit48':
      return [[]]

    val = []
    fmtlist = self.endian + mychar(fmt)
    size = mysize(fmt)

    for _ in range(length):
       string = self.file.read(size)
       val.append(_struct.unpack(fmtlist,string)[0])

    return val

  def scan(self,fmtlist,origin=0,contiguous=1):
    # reads binary file
    self.file.seek(origin)

    position = 0
    h = structure()

    for item in fmtlist:
      fmt = item[0]
      length = item[1]
      offset = item[2]
      name = item[3]

      if not contiguous:
        self.file.seek(offset-position,1)
        position = offset+mysize(fmt)

      #if length is 1:
      #  h[name] = self.read(fmt,length)[0]
      #else:
      #  h[name] = self.read(fmt,length)
      h[name] = self.read(fmt,length)[0]

    return h


class Writer(object):
  def __init__(self,fname,endian='|'):
    # open binary file
    self.file = open(fname,'w')
    path = _os.path.abspath(self.file.name)
    self.path = _os.path.dirname(path)
    self.name = _os.path.basename(path)
    self.size = _os.path.getsize(path)
    self.endian = endian

  def __del__(self):
    self.file.close()

  def write(self,fmt,vals,length=1,offset=0):
    # write binary data
    if offset != 0:
      self.file.seek(offset)

    if fmt is 'bit48':
      return [[]]

    fmtlist = self.endian + mychar(fmt)

    if length==1:
       vals = [vals]

    for val in vals:
       self.file.write(_struct.pack(fmtlist,val))

  def printf(self,fmts,vals,origin=0,contiguous=1):
    # write binary data
    self.file.seek(origin)
    position = 0

    for i,item in enumerate(fmts):
      fmt = item[0]
      length = item[1]
      offset = item[2]
      val = vals[i]

      if not contiguous:
        self.file.seek(offset-position,1)
        position = offset+mysize(fmt)

      self.write(fmt,val,length)


class structure(dict):
  def __init__(self,*args,**kwargs):
    super(structure,self).__init__(*args,**kwargs)
    self.__dict__ = self


def mychar(fmt):
  chars = { 'int8'   :'b',
        'uint8'  :'B',
        'int16'  :'h',
        'uint16' :'H',
        'int32'  :'i',
        'uint32' :'I',
        'int64'  :'q',
        'uint64' :'Q',
        'float'  :'f',
        'float32':'f',
        'double' :'d',
        'char'   :'s'}
  if fmt in chars:
    return chars[fmt]
  else:
    return fmt


def mysize(fmt):
  return _struct.calcsize(mychar(fmt))

