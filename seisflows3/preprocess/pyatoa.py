#!/usr/bin/env python
"""
This is the base class seisflows.preprocess.Pyatoa

This is a main Seisflows class, it controls the preprocessing.
This class uses the Python package Pyatoa to perform preprocessing, and
misfit measurement.

..warning::
    This might break if no residuals are written for a given event
"""
import os
import sys
import pyatoa
import numpy as np
from glob import glob
from seisflows.tools import unix
from seisflows.config import custom_import
from seisflows.config import SeisFlowsPathsParameters
from pyatoa.utils.images import merge_pdfs

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]


class Pyatoa(custom_import("preprocess", "base")):
    """
    Data preprocessing class using the Pyatoa package
    """
    def __init__(self):
        """
        These parameters should not be set by __init__!
        Attributes are just initialized as NoneTypes for clarity and docstrings

        :type data: str
        :param data: directory where data from the preprocessing is stored
        :type figures: str
        :param figures: directory where figures are stored
        """
        self.path_datasets = None
        self.path_figures = None

    @property
    def required(self):
        """
<<<<<<< HEAD
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters()

        # Define the Parameters required by this module
        sf.par("UNIT_OUTPUT", required=True, par_type=str,
               docstr="Data units. Must match the synthetic output of external "
                      "solver. Available: ['DISP': displacement, "
                      "'VEL': velocity, 'ACC': acceleration]")

        sf.par("MIN_PERIOD", required=True, par_type=float,
               docstr="Minimum filter corner in seconds")

        sf.par("MAX_PERIOD", required=True, par_type=float,
               docstr="Maximum filter corner in seconds")

        sf.par("CORNERS", required=False, default=4, par_type=int,
               docstr="Number of filter corners")

        sf.par("CLIENT", required=False, par_type=str,
               docstr="Client name for ObsPy FDSN data gathering")

        sf.par("START_PAD", required=False, default=0, par_type=float,
               docstr="For data gathering; time before origin time to gather. "
                      "START_PAD >= T_0 in SPECFEM constants.h.in. "
                      "Positive values only")

        sf.par("END_PAD", required=True, par_type=float,
               docstr="For data gathering; time after origin time to gather. "
                      "END_PAD >= NT * DT (of Par_file). Positive values only")

        sf.par("ROTATE", required=False, default=False, par_type=bool,
               docstr="Rotate waveform components NEZ -> RTZ")
=======
        # Check the path requirements
        if "PREPROCESS" not in PATH:
            setattr(PATH, "PREPROCESS", 
                    os.path.join(PATH.SCRATCH, "preprocess"))

        if "DATA" not in PATH:
            setattr(PATH, "DATA", None)

        if "RESPONSE" not in PATH:
            setattr(PATH, "RESPONSE", None)

        # Check the existence of required parameters
        required_parameters = ["COMPONENTS", "UNIT_OUTPUT", "MIN_PERIOD",
                               "MAX_PERIOD", "CORNERS", "ROTATE",
                               "ADJ_SRC_TYPE", "FIX_WINDOWS", "PLOT", "FORMAT"
                               ]
        for req in required_parameters:
            if req not in PAR:
                raise ParameterError(PAR, req)

        # Check specific parameter requirements
        if PAR.FORMAT != "ascii":
            raise ValueError("Pyatoa preprocess currently only works with "
                             "the 'ascii' format")

        sf.par("ADJ_SRC_TYPE", required=True, par_type=str,
               docstr="Adjoint source type to use. Available: "
                      "['cc': cross-correlation, 'mt': multitaper, "
                      "wav: waveform']")

        sf.par("PYFLEX_PRESET", required=True, par_type=str,
               docstr="Parameter map for Pyflex config. For available choices, "
                      "see Pyatoa docs page")

        sf.par("FIX_WINDOWS", required=False, default=False,
               par_type="bool or str",
               docstr="Time window evaluation: available: "
                      "[True, False, 'ITER', 'ONCE'] "
                      "True: Same windows for all but i01s00; "
                      "False: New windows at each evaluation; "
                      "'ITER': New windows at first evaluation of each iter; "
                      "'ONCE': New windows at first evaluation of workflow")

        sf.par("PLOT", required=False, default=True, par_type=bool,
               docstr="Plot waveforms and maps as .pdf files")

        sf.par("SNAPSHOT", required=False, default=True, par_type=bool,
               docstr="Copy ASDFDataSets on disk for data redundancy")

        sf.par("LOGGING", required=False, default="DEBUG", par_type=str,
               docstr="Log level. Available: "
                      "['null': no logging, 'warning': warnings only, "
                      "'info': task tracking, 'debug': log everything]")

        # Define the Paths required by this module
        sf.path("PREPROCESS", required=False,
                default=os.path.join(PATH.SCRATCH, "preprocess"),
                docstr="scratch path to store waveform data and figures")

        sf.path("DATA", required=False,
                docstr="Directory to locally stored data")

        return sf

    def check(self, validate=True):
        """ 
        Checks Parameter and Path files, will be run at the start of a Seisflows
        workflow to ensure that things are set appropriately.
        """
        if validate:
            self.required.validate()

        # Check that other modules have set parameters that will be used here
        for required_parameter in ["COMPONENTS", "FORMAT"]:
            assert(required_parameter in PAR), \
                f"Pyatoa requires {required_parameter}"


        if PAR.FORMAT != "ascii":
            raise ValueError("Pyatoa preprocess currently only works with "
                             "the 'ascii' format")

        if PAR.DT * PAR.NT >= PAR.START_PAD + PAR.END_PAD:
            raise ValueError("Pyatoa preprocess parameters START_PAD and "
                             "END_PAD will not provide long enough obs."
                             "traces to match the length of synthetics")

    def setup(self):
        """
        Sets up data preprocessing machinery by establishing an internally
        defined directory structure that will be used to store the outputs 
        of the preprocessing workflow

        Akin to an __init__ class, but to be called externally by the workflow.
        """
        # Late import because preprocess is loaded before optimize
        solver = sys.modules["seisflows_solver"]

        # Inititate a Pyaflowa object to make sure the machinery works
        pyaflowa = pyatoa.Pyaflowa(structure="seisflows", sfpaths=PATH, 
                                   sfpar=PAR)

        # Pull path names from Pyaflowa to keep path structure in one place
        self.path_datasets = pyaflowa.path_structure.datasets
        self.path_figures = pyaflowa.path_structure.figures

    def prepare_eval_grad(self, path, source_name, **kwargs):
        """
        Prepare the gradient evaluation by gathering, preprocessing waveforms, 
        and measuring misfit between observations and synthetics using Pyatoa.
        
        This is a process specific task and intended to be run in parallel

        :type path: str
        :param path: path to the current function evaluation for saving residual
        :type source_name: str
        :param source_name: the event id to be used for tagging and data lookup
        """
        # Late import because preprocess is loaded before optimize,
        # Optimize required to know which iteration/step_count we are at
        optimize = sys.modules["seisflows_optimize"]

        # Inititate the Pyaflowa class which abstracts processing functions
        # Communicate to Pyaflowa the current iteration and step count
        pyaflowa = pyatoa.Pyaflowa(structure="seisflows", sfpaths=PATH, 
                                   sfpar=PAR, iteration=optimize.iter,
                                   step_count=optimize.line_search.step_count)

        # Process all the stations for a given event using Pyaflowa
        misfit = pyaflowa.process_event(source_name, 
                                        fix_windows=PAR.FIX_WINDOWS)

        # Generate the necessary files to continue the inversion
        if misfit:
            # Event misfit defined by Tape et al. (2010)
            self.write_residuals(path=path, scaled_misfit=misfit,
                                 source_name=source_name)
        
        # self.snapshot()

    def finalize(self):
        """
        Run some serial finalization tasks specific to Pyatoa, which will help
        aggregate the collection of output information:
            Aggregate misfit windows using the Inspector class
            Generate PDFS of waveform figures for easy access
        """
        unix.cd(self.path_datasets)
        insp = pyatoa.Inspector(PAR.TITLE, verbose=False)
        insp.discover()
        insp.save() 

        self.make_final_pdfs()

    def write_residuals(self, path, scaled_misfit, source_name):
        """
        Computes residuals and saves them to a text file in the appropriate path

        :type path: str        
        :param path: scratch directory path, e.g. PATH.GRAD or PATH.FUNC
        :type scaled_misfit: float
        :param scaled_misfit: the summation of misfit from each 
            source-receiver pair calculated by prepare_eval_grad()
        :type source_name: str
        :param source_name: name of the source related to the misfit, used
            for file naming
        """
        residuals_dir = os.path.join(path, "residuals")        

        if not os.path.exists(residuals_dir):
            unix.mkdir(residuals_dir)
        
        event_residual = os.path.join(residuals_dir, source_name)        
     
        np.savetxt(event_residual, [scaled_misfit], fmt="%11.6e")

    def sum_residuals(self, files):
        """
        Averages the event misfits and returns the total misfit.
        Total misfit defined by Tape et al. (2010)

        :type files: str
        :param files: list of single-column text files containing residuals
            that will have been generated using prepare_eval_grad()
        :rtype: float
        :return: average misfit
        """
        assert(len(files) == PAR.NTASK), \
            "Number of misfit files does not match the number of events"

        total_misfit = 0
        for filename in files:
            total_misfit += np.sum(np.loadtxt(filename))

        total_misfit /= PAR.NTASK

        return total_misfit

    def snapshot(self):
        """
        Copy all ASDFDataSets in the data directory into a separate snapshot
        directory for redundancy
        """
        if PAR.SNAPSHOT:
            snapshot_dir = os.path.join(self.path_datasets, "snapshot")
            if not os.path.exists(snapshot_dir):
                unix.mkdir(snapshot_dir)
            
            srcs = glob(os.path.join(self.path_datasets, "*.h5"))
            for src in srcs:
                dst = os.path.join(snapshot_dir, os.path.basename(src))
                unix.cp(src, dst)

    def make_final_pdfs(self):
        """
        Utility function to combine all pdfs for a given event, iteration, and
        step count into a single pdf. To reduce on file count and provide easier
        visualization. Removes the original event-based pdfs.
        
        .. warning::
            This is a simple function because it won't account for missed 
            iterations i.e. if this isn't run in the finalization, it will 
            probably break the next time

        :raises AssertionError: When tags don't match the mainsolvers first tag
        """
        # Late import because preprocess is loaded before optimize
        solver = sys.modules["seisflows_solver"]

        # Relative pathing from here on out boys
        unix.cd(self.path_figures)
        sources = []
        for source_name in solver.source_names:
            sources += glob(os.path.join(source_name, "*.pdf"))

        # Incase this is run out of turn and pdfs were already deleted
        if not sources:
            return

        iter_steps = set([os.path.basename(_).split("_")[0] for _ in sources])
        for iter_step in iter_steps:
            # Merge pdfs that correspond to the same iteration and step count 
            fids = [_ for _ in sources if iter_step in _]
            merge_pdfs(fids=sorted(fids), fid_out=f"{iter_step}.pdf")
            unix.rm(fids)


