import os
import sys

src_path = os.path.dirname(os.path.realpath(__file__))
src_path += '/../../..'
if src_path not in sys.path:
    sys.path.append(os.path.abspath(src_path))
