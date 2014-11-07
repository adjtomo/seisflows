
### data processing

from signal import sbandpass, smute, swindow

from segy import reader as segyreader
from segy import writer as segywriter

from segy.reader import readsegy, readsu
from segy.writer import writesegy, writesu



### adjoint tomography

import adjoint, misfit



### forward modeling

import specfem2d
import specfem3d
import specfem3d_globe



### visualization

from graphics import splot
from graphics import wplot
