#!/usr/bin/env python3
"""
The Pyatoa preprocessing module abstracts all preprocessing functionality
onto Pyatoa (https://github.com/bch0w/pyatoa/). The module defined below is
meant to set up and execute Pyatoa within a running SeisFlows3 workflow.

Pyatoa itself aggregates all of its connection with SeisFlows3 in the Pyaflowa
class, a purpose built object used to simplify calling Pyatoa from within
a SeisFlows3 workflow.
"""
import os
import sys
import logging
import numpy as np
from glob import glob

from pyatoa import Pyaflowa, Inspector

from seisflows3.tools import unix, msg
from seisflows3.config import custom_import, SeisFlowsPathsParameters, CFGPATHS

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]


class Pyatoa(custom_import("preprocess", "base")):
    """
    Data preprocessing class using the Pyaflowa class within the Pyatoa package.
    In charge of data discovery, preprocessing, filtering, misfiti
    quantification and data storage. The User does not need to implement Pyatoa,
    but rather interacts with it via the parameters and paths of SeisFlows3.
    """
    logger = logging.getLogger(__name__).getChild(__qualname__)

    def __init__(self):
        """
        These parameters should not be set by __init__!
        Attributes are just initialized as NoneTypes for clarity and docstrings

        :param logger: Class-specific logging module, log statements pushed
            from this logger will be tagged by its specific module/classname
        """
        pass

    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters()

        # Define the Parameters required by this module
        sf.par("UNIT_OUTPUT", required=True, par_type=str,
               docstr="Data units. Must match the synthetic output of external "
                      "solver. Available: ['DISP': displacement, "
                      "'VEL': velocity, 'ACC': acceleration]")

        # TODO Check this against T0 in check()
        sf.par("START_PAD", required=False, default=0, par_type=float,
               docstr="For data gathering; time before origin time to gather. "
                      "START_PAD >= T_0 in SPECFEM constants.h.in. "
                      "Positive values only")

        # TODO set this automatically by setting equal NT * DT
        sf.par("END_PAD", required=True, par_type=float,
               docstr="For data gathering; time after origin time to gather. "
                      "END_PAD >= NT * DT (of Par_file). Positive values only")

        sf.par("MIN_PERIOD", required=False, default="", par_type=float,
               docstr="Minimum filter corner in unit seconds. Bandpass filter "
                      "if set with `MAX_PERIOD`, highpass filter if set "
                      "without `MAX_PERIOD`, no filtering if not set and "
                      "`MAX_PERIOD also not set")

        sf.par("MAX_PERIOD", required=False, default="", par_type=float,
               docstr="Maximum filter corner in unit seconds. Bandpass filter "
                      "if set with `MIN_PERIOD`, lowpass filter if set "
                      "without `MIN_PERIOD`, no filtering if not set and "
                      "`MIN_PERIOD also not set")

        sf.par("CORNERS", required=False, default=4, par_type=int,
               docstr="Number of filter corners applied to filtering")

        sf.par("CLIENT", required=False, par_type=str,
               docstr="Client name for ObsPy FDSN data gathering. Pyatoa will "
                      "attempt to collect waveform and metadata based on "
                      "network and station codes provided in the SPECFEM "
                      "STATIONS file. If set None, no FDSN gathering will be "
                      "attempted")

        sf.par("ROTATE", required=False, default=False, par_type=bool,
               docstr="Attempt to rotate waveform components from NEZ -> RTZ")

        sf.par("PYFLEX_PRESET", required=False, default="default", 
               par_type=str,
               docstr="Parameter map for misfit window configuration defined "
                      "by Pyflex. IF None, misfit and adjoint sources will be "
                      "calculated on whole traces. For available choices, "
                      "see Pyatoa docs page (pyatoa.rtfd.io)")

        sf.par("FIX_WINDOWS", required=False, default=False,
               par_type="bool or str",
               docstr="How to address misfit window evaluation at each "
                      "evaluation. Options to re-use misfit windows collected "
                      "during an inversion, available options: "
                      "[True, False, 'ITER', 'ONCE'] "
                      "True: Re-use windows after first evaluation (i01s00); "
                      "False: Calculate new windows each evaluation; "
                      "'ITER': Calculate new windows at first evaluation of "
                      "each iteration (e.g., i01s00... i02s00..."
                      "'ONCE': Calculate new windows at first evaluation of "
                      "the workflow, i.e., at PAR.BEGIN")

        sf.par("ADJ_SRC_TYPE", required=False, default="cc",  par_type=str,
               docstr="Adjoint source type to evaluate misfit, defined by "
                      "Pyadjoint. Currently available options: "
                      "['cc': cross-correlation, 'mt': multitaper, "
                      "wav: waveform']")

        sf.par("PLOT", required=False, default=True, par_type=bool,
               docstr="Attempt to plot waveforms and maps as PDF files at each "
                      "function evaluation")

        sf.par("PYATOA_LOG_LEVEL", required=False, default="DEBUG", 
               par_type=str,
               docstr="Log level to set Pyatoa, Pyflex, Pyadjoint. Available: "
                      "['null': no logging, 'warning': warnings only, "
                      "'info': task tracking, "
                      "'debug': log all small details (recommended)]")

        # Parameters to control saving scratch/preprocess files to work dir.
        sf.par("SAVE_DATASETS", required=False, default=True, par_type=bool,
               docstr="Save PyASDF HDF5 datasets to disk. These datasets store "
                      "waveform data, metadata, misfit windows, adjoint "
                      "sources and configuration parameters")

        sf.par("SAVE_FIGURES", required=False, default=True, par_type=bool,
               docstr="Save output waveform figures to disk as PDFs")

        sf.par("SAVE_LOGS", required=False, default=True, par_type=bool,
               docstr="Save event-specific Pyatoa logs to disk as .txt files")

        # Define the Paths required by this module
        sf.path("PREPROCESS", required=False,
                default=os.path.join(PATH.SCRATCH, "preprocess"),
                docstr="scratch/ path to store waveform data and figures. "
                       "Pyatoa will generate an internal directory structure "
                       "here")

        sf.path("DATA", required=False,
                docstr="Directory to locally stored data. Pyatoa looks for "
                       "waveform and metadata in the 'PATH.DATA/mseed' and "
                       "'PATH.DATA/seed', directories respectively.")

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

        assert(PAR.FORMAT.upper() == "ASCII"), \
            "Pyatoa preprocess requires PAR.FORMAT=='ASCII'"

        assert((PAR.DT * PAR.NT) <= (PAR.START_PAD + PAR.END_PAD)), \
            ("Pyatoa preprocess must have (PAR.START_PAD + PAR.END_PAD) >= "
             "(PAR.DT * PAR.NT), current values will not provide sufficiently "
             f"long data traces (DT*NT={PAR.DT * PAR.NT}; "
             f"START+END={PAR.START_PAD + PAR.END_PAD}")

    def setup(self):
        """
        Sets up data preprocessing machinery by establishing an internally
        defined directory structure that will be used to store the outputs 
        of the preprocessing workflow

        Akin to an __init__ class, but to be called externally by the workflow.
        """
        unix.mkdir(PATH.PREPROCESS)

    def prepare_eval_grad(self, cwd, source_name, taskid, **kwargs):
        """
        Prepare the gradient evaluation by gathering, preprocessing waveforms, 
        and measuring misfit between observations and synthetics using Pyatoa.

        Reads in observed and synthetic waveforms, applies optional
        preprocessing, assesses misfit, and writes out adjoint sources and
        STATIONS_ADJOINT file.

        .. note::
            Meant to be called by solver.eval_func(), may have unused arguments
            to keep functions general across preprocessing subclasses.

        :type cwd: str
        :param cwd: current specfem working directory containing observed and
            synthetic seismic data to be read and processed. Should be defined
            by solver.cwd
        :type source_name: str
        :param source_name: the event id to be used for tagging and data lookup.
            Should be defined by solver.source_name
        :type taskid: int
        :param taskid: identifier of the currently running solver instance.
            Should be defined by solver.taskid
        :type filenames: list of str
        :param filenames: [not used] list of filenames defining the files in
            traces
        """
        if taskid == 0:
            self.logger.debug("preparing files for gradient evaluation with "
                              "Pyaflowa")

        # Process all the stations for a given event using Pyaflowa
        pyaflowa = self.setup_event_pyaflowa(source_name)
        scaled_misfit = pyaflowa.process(nproc=PAR.NPROC)

        if scaled_misfit is None:
            print(msg.cli(f"Event {source_name} returned no misfit, you may "
                          f"want to check logs and waveform figures, "
                          f"or consider discarding this event from your "
                          f"workflow", 
                          items=[pyaflowa.paths.logs, pyaflowa.paths.figures],
                          header="pyatoa preprocessing error", border="="))
            sys.exit(-1)

        # Event misfit defined by Tape et al. (2010) written to solver dir.
        self.write_residuals(path=cwd, scaled_misfit=scaled_misfit)

    def setup_event_pyaflowa(self, source_name=None):
        """
        A convenience function to set up a Pyaflowa processing instance for
        a specific event. 

        .. note::
            This is meant to be called by preprocess.prepare_eval_grad() but its
            also useful for debugging and manual processing where you can simply
            return a formatted Pyaflowa object and debug it directly.

        :type source_name: str
        :param source_name: solver source name to evaluate setup for. Must 
            match from list defined by: solver.source_names
        """
        # Late import because preprocess is loaded before optimize,
        # Optimize required to know which iteration/step_count we are at
        solver = sys.modules["seisflows_solver"]
        optimize = sys.modules["seisflows_optimize"]

        iteration = optimize.iter
        if source_name is None:
            source_name = solver.source_names[0]

        # Deal with the migration case where no step count given
        try:
            step_count = optimize.line_search.step_count
        except AttributeError:
            step_count = ""

        # Outsource data processing to an event-specfic Pyaflowa instance
        pyaflowa = Pyaflowa(sfpar=PAR, sfpath=PATH)
        pyaflowa.setup(source_name=source_name, iteration=iteration, 
                       step_count=step_count, loc="*", cha="*")
        
        return pyaflowa

    def finalize(self):
        """
        Run some serial finalization tasks specific to Pyatoa, which will help
        aggregate the collection of output information.

        .. note::
            This finalize function performs the following tasks:
            * Generate .csv files using the Inspector
            * Aggregate event-specific PDFs into a single evaluation PDF
            * Save scratch/ data into output/ if requested
        """
        # Initiate Pyaflowa to get access to path structure
        pyaflowa = Pyaflowa(sfpar=PAR, sfpath=PATH)
        unix.cd(pyaflowa.paths.datasets)

        # Generate the Inspector from existing datasets and save to disk
        # Allow this is fail, which might happen if we don't have enough data
        # or the Dataset is not formatted as expected
        insp = Inspector(PAR.TITLE, verbose=False)
        try:
            insp.discover()
            insp.save()
        except Exception as e:
            self.logger.warning(f"Uncontrolled exception in Inspector creation "
                                f"will not create inspector:\n{e}")

        # Make the final PDF for easier User ingestion of waveform/map figures
        pyaflowa.make_evaluation_composite_pdf()

        # Move scratch/ directory results into more permanent storage
        if PAR.SAVE_DATASETS:
            datasets = glob(os.path.join(pyaflowa.paths.datasets, "*.h5"))
            self._save_quantity(datasets, tag="datasets")
        
        if PAR.SAVE_FIGURES:
            figures = glob(os.path.join(pyaflowa.paths.figures, "*.pdf"))
            self._save_quantity(figures, tag="figures")

        if PAR.SAVE_LOGS:
            logs = glob(os.path.join(pyaflowa.paths.logs, "*.txt"))
            path_out = os.path.join(PATH.WORKDIR, CFGPATHS.LOGDIR)
            self._save_quantity(logs, path_out=path_out)
    
    def _save_quantity(self, filepaths, tag="", path_out=""):
        """
        Repeatable convenience function to save quantities from the scratch/
        directory to the output/ directory

        :type filepaths: list
        :param filepaths: full path to files that should be saved to output/
        :type tag: str  
        :param tag: tag for saving the files in PATH.OUTPUT. If not given, will
            save directly into the output/ directory
        :type path_out: str
        :param path_out: overwrite the default output path file naming
        """       
        if not path_out:
            path_out = os.path.join(PATH.OUTPUT, tag)

        if not os.path.exists(path_out):
            unix.mkdir(path_out)

        for src in filepaths:
            dst = os.path.join(path_out, os.path.basename(src))
            unix.cp(src, dst) 

    def write_residuals(self, path, scaled_misfit):
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
        residuals_file = os.path.join(path, "residuals")        
        np.savetxt(residuals_file, [scaled_misfit], fmt="%11.6e")

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
        if len(files) != PAR.NTASK:
            print(msg.cli(f"Pyatoa preprocessing module did not recover the "
                          f"correct number of residual files "
                          f"({len(files)}/{PAR.NTASK}). Please check that "
                          f"the preprocessing logs", header="error")
                  )
            sys.exit(-1)

        total_misfit = 0
        for filename in files:
            total_misfit += np.sum(np.loadtxt(filename))

        total_misfit /= PAR.NTASK

        return total_misfit

