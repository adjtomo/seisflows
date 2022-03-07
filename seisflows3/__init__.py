import logging
from pkgutil import extend_path


# Extend the search path for the modules which comprise a package. 
# This will add to the package’s __path__ all subdirectories of directories on 
# sys.path named after the package. This is useful if one wants to distribute 
# different parts of a single logical package as multiple directories.
__path__ = extend_path(__path__, __name__)


# Set up the SeisFlows3 Logging environment
logger = logging.getLogger("seisflows3")
logger.setLevel(logging.INFO)  # default level 
ch = logging.StreamHandler()
fmt_str = "[%(asctime)s] %(levelname)s: %(message)s"
formatter = logging.Formatter(fmt_str, datefmt="%Y-%m-%d %H:%M:%S")
ch.setFormatter(formatter)
logger.addHandler(ch)
