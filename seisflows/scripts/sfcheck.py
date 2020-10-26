"""
A suite of functions used to check on the status of a running or to-be-run 
SeisFlows workflow
""" 
import os
import sys
import argparse
from glob import glob

from seisflows.config import Dict, tilde_expand, names
from seisflows.tools.tools import loadyaml, loadobj
from seisflows.scripts.seisflows import parse_null


class Checker:
    """
    Used to check on information from an ongoing SeisFlows workflow and return
    relevant information to the user using print statements.
    """
    def __init__(self, args):
        """
        Ensure all modules are available to the class. Init requires a specific
        order of arguments, that is, function name in 0th position, and all
        arguments passed to that function following.
        """
        # Match internal attributes to the SeisFlows naming convention
        self.PATH = sys.modules["seisflows_paths"]
        self.PAR = sys.modules["seisflows_parameters"]
        self.system = sys.modules["seisflows_system"]
        self.solver = sys.modules["seisflows_solver"]
        self.optimize = sys.modules["seisflows_optimize"]
        self.preprocess = sys.modules["seisflows_preprocess"]
        self.postprocess = sys.modules["seisflows_postprocess"]

        getattr(self, args[0])(*args[1:])

    def model(self, src=None):
        """
        Print out the min/max values from one or all of the currently available
        models. Useful for checking what models are associated with what part of
        the workflow, e.g. evaluate function, evaluate gradient.
        """
        avail = glob(os.path.join(self.PATH.OPTIMIZE, "m_*"))
        srcs = [os.path.basename(_) for _ in avail]
        if src:
            assert(src in srcs), f"{src} not in available models {avail}"
            srcs = [src]
        for tag in srcs:
            m = self.optimize.load(tag)
            self.optimize.check_model_parameters(m, tag)


def main():
    """
    Entry point to the Checker class. Basically need to initate all of SeisFlows
    and then run whatever check function is requested.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("main_args", type=str, nargs="*")
    parser.add_argument("-w", "--workdir", nargs="?", default=os.getcwd())
    parser.add_argument("-p", "--parameter_file", nargs="?",
                        default="parameters.yaml")
    
    args = parser.parse_args()
    if not os.path.exists(args.parameter_file):
        raise FileNotFoundError(f"Parameter file not found: "
                                f"{args.parameter_file}")
    
    # Register parameters
    parameters = loadyaml(args.parameter_file)
    parameters = parse_null(parameters)
    sys.modules["seisflows_parameters"] = Dict(parameters)
        
    # Register paths
    paths = tilde_expand(parameters["PATHS"])
    paths = {key: os.path.abspath(path) for key, path in paths.items()}
    sys.modules["seisflows_paths"] = Dict(paths)

    # Register modules
    for name in names:
        fullfile = os.path.join(args.workdir, "output", f"seisflows_{name}.p")
        sys.modules[f"seisflows_{name}"] = loadobj(fullfile)

    # Initiate defaults, must be run after all modules are registered
    for name in names:
        sys.modules[f"seisflows_{name}"].check() 

    # Run the checker
    if args.main_args:
        checker = Checker(args.main_args)  


if __name__ == "__main__":
    main()
    


    
