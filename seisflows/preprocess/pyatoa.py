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
from seisflows.config import Dict
from pyatoa.utils.images import merge_pdfs
from seisflows.tools.err import ParameterError

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]


class Pyatoa:
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
        self.data = None
        self.figures = None

    def check(self):
        """ 
        Checks Parameter and Path files, will be run at the start of a Seisflows
        workflow to ensure that things are set appropriately.
        """
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
                               "ADJ_SRC_TYPE", "PYFLEX_PRESET",
                               "FIX_WINDOWS", "PLOT", "FORMAT"
                               ]
        for req in required_parameters:
            if req not in PAR:
                raise ParameterError(PAR, req)

        # Check specific parameter requirements
        if PAR.FORMAT != "ascii":
            raise ValueError("Pyatoa preprocess currently only works with "
                             "the 'ascii' format")

        # Set default values parameters for any non-set parameters
        if "PLOT" not in PAR:
            setattr(PAR, "PLOT", True)

        if "LOGGING" not in PAR:
            setattr(PAR, "LOGGING", "DEBUG")

        if "MAP_CORNERS" not in PAR:
            setattr(PAR, "MAP_CORNERS", None)

        if "CLIENT" not in PAR:
            setattr(PAR, "CLIENT", None)

        if "SNAPSHOT" not in PAR:
            setattr(PAR, "SNAPSHOT", True)

        # Used to define the start time of fetched observation waveforms
        if "START_PAD" not in PAR:
            setattr(PAR, "START_PAD", 20)

        # Used to define the end time of fetched observation waveforms
        if "END_PAD" not in PAR:
            setattr(PAR, "END_PAD", PAR.DT * PAR.NT + PAR.START_PAD + 5)
        else:
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

    def prepare_eval_grad(self, path, source_name):
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
        pyaflowa = pyatoa.Pyaflowa(structure="seisflows", sfpaths=self.paths, 
                                   sfpar=PAR, iteration=optimize.iter,
                                   step_count=optimize.line_search.step_count)

        # Process all the stations for a given event using Pyaflowa
        misfit = pyaflowa.process(source_name, fix_windows=PAR.FIX_WINDOWS)

        # Generate the necessary files to continue the inversion
        if misfit:
            # Event misfit defined by Tape et al. (2010)
            self.write_residuals(path=path, scaled_misfit=misfit,
                                 source_name=source_name)
        
        self.snapshot()

    def finalize(self):
        """
        Run some serial finalization tasks specific to Pyatoa, which will help
        aggregate the collection of output information:
            Aggregate misfit windows using the Inspector class
            Generate PDFS of waveform figures for easy access
        """
        insp = pyatoa.Inspector(PAR.TITLE, verbose=False)
        insp.discover(path=os.path.join(self.paths.data, PAR.TITLE))
        insp.save(path=self.paths.data) 

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
            srcs = glob(os.path.join(self.paths.data, "*.h5"))
            for src in srcs:
                dst = os.path.join(self.paths.data, "snapshot")
                unix.cp(src, dst)

    def make_final_pdfs(self):
        """
        Utility function to combine all pdfs for a given event, iteration, and
        step count into a single pdf. To reduce on file count and provide easier
        visualization. Removes the original event-based pdfs
        """
        unix.cd(self.paths.figures)

        # This glob list contains all pdfs tagged e.g.
        # '{iter}{step}_{event_id}.pdf', we need to break apart by iter and step
        sources = glob("*.pdf")
        iterstep_tags = set([_.split("_")[0] for _ in sources])
        for is_tag in iterstep_tags:
            event_pdfs = glob(f"{is_tag}_*.pdf")
            iter_ = is_tag[:3]  # e.g. i01
            step_ = is_tag[3:]  # e.g. s00

            if not os.path.exists(iter_):
                os.makedirs(iter_)

            merge_pdfs(fids=event_pdfs,
                       fid_out=os.path.join(iter_, f"{iter_}{step_}.pdf")
                       )
            unix.rm(event_pdfs)



