import os
import logging

from pkgutil import extend_path

# Extend the search path for the modules which comprise a package. 
# This will add to the package’s __path__ all subdirectories of directories on 
# sys.path named after the package. This is useful if one wants to distribute 
# different parts of a single logical package as multiple directories.
__path__ = extend_path(__path__, __name__)

# Set up the SeisFlows Logging environment
logger = logging.getLogger(__name__)

# The location of this config file, which is the main repository
ROOT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))

# List of module names required by SeisFlows for imports. Order-sensitive
NAMES = ["workflow", "system", "solver", "preprocess", "optimize"]