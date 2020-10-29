#!/usr/bin/env python
"""
This is the base class seisflows.workflow.migration

This is a main Seisflows class, it controls the main workflow.
"""
import sys

from seisflows.tools import unix
from seisflows.tools.tools import exists
from seisflows.tools.err import ParameterError
from seisflows.workflow.base import Base
from seisflows.config import SeisFlowsPathsParameters


PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
solver = sys.modules['seisflows_solver']
preprocess = sys.modules['seisflows_preprocess']
postprocess = sys.modules['seisflows_postprocess']


class Migration(Base):
    """
    Migration base class.

    Performs the workflow of an inversion up to the postprocessing. In the
    terminology of seismic exploration, implements a 'reverse time migration'.
    """
    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters(super().required)

        # Define the Paths required by this module
        sf.path("SCRATCH", required=True,
                docstr="scratch path to hold temporary data during workflow")

        sf.path("OUTPUT", required=True,
                docstr="directory to save workflow outputs to disk")

        sf.path("LOCAL", required=False, default="null",
                docstr="local path to data available to workflow")

        sf.path("DATA", required=False, default="null",
                docstr="path to data available to workflow")

        sf.path("MODEL_INIT", required=True,
                docstr="location of the initial model to be used for workflow")

        sf.par("SAVEGRADIENT", required=False, default=True, par_type=bool,
               docstr="Save gradient files after each iteration")

        sf.par("SAVEKERNELS", required=False, default=False, par_type=bool,
               docstr="Save event kernel files after each iteration")

        sf.par("SAVETRACES", required=False, default=False, par_type=bool,
               docstr="Save waveform traces after each iteration")

        return sf

    def check(self, validate=True):
        """ 
        Checks parameters and paths
        """
        if validate:
            self.required.validate()
        super().check(validate=False)

        if not exists(PATH.DATA):
            assert "MODEL_TRUE" in PATH, f"DATA or MODEL_TRUE must exist"

    def main(self):
        """ Migrates seismic data
        """
        # set up workflow machinery
        preprocess.setup()
        postprocess.setup()

        # set up solver machinery
        print "Preparing solver..."
        system.run("solver", "setup")

        self.prepare_model()

        # perform migration
        print "Generating synthetics..."
        system.run("solver", "eval_func",
                   path=PATH.SCRATCH,
                   write_residuals=False)

        print "Backprojecting..."
        system.run("solver", "eval_grad",
                   path=PATH.SCRATCH,
                   export_traces=PAR.SAVETRACES)

        system.run_single("postprocess", "process_kernels",
                 path=PATH.SCRATCH+"/"+"kernels",
                 parameters=solver.parameters)

        try:
            system.run_single("postprocess", "process_kernels",
                     path=PATH.SCRATCH+"/"+"kernels",
                     parameters=["rhop"])
        except:
            pass

        if PAR.SAVETRACES:
            self.save_traces()

        if PAR.SAVEKERNELS:
            self.save_kernels()
        else:
            self.save_kernels_sum()

        print "Finished\n"

    def prepare_model(self):
        model = PATH.OUTPUT +"/"+ "model_init"
        assert exists(model)
        unix.cp(model, PATH.SCRATCH +"/"+ "model")

    def save_kernels_sum(self):
        src = PATH.SCRATCH +"/"+ "kernels/sum"
        dst = PATH.OUTPUT +"/"+ "kernels"
        unix.mkdir(dst)
        unix.cp(src, dst)

    def save_kernels(self):
        src = PATH.SCRATCH +"/"+ "kernels"
        dst = PATH.OUTPUT
        unix.mkdir(dst)
        unix.cp(src, dst)

    def save_traces(self):
        src = PATH.SCRATCH +"/"+ "traces"
        dst = PATH.OUTPUT
        unix.cp(src, dst)

