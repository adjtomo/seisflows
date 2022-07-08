#!/usr/bin/env python3
"""
The Pyatoa preprocessing module abstracts all preprocessing functionality
onto Pyatoa (https://github.com/adjtomo/pyatoa/). The module defined below is
meant to set up and execute Pyatoa within a running SeisFlows workflow.

Pyatoa itself aggregates all of its connection with SeisFlows in the Pyaflowa
class, a purpose built object used to simplify calling Pyatoa from within
a SeisFlows workflow.
"""
import os
import sys
import numpy as np
from glob import glob
from pyatoa import Pyaflowa, Inspector

from seisflows import logger
from seisflows.tools import unix, msg
from seisflows.config import CFGPATHS


class Pyatoa:
    """
    Data preprocessing class using the Pyaflowa class within the Pyatoa package.
    In charge of data discovery, preprocessing, filtering, misfiti
    quantification and data storage. The User does not need to implement Pyatoa,
    but rather interacts with it via the parameters and paths of SeisFlows.
    """
    def __init__(self, data_format="ascii", components=None, ntask=1, nproc=1,
                 min_period=None, max_period=None, filter_corners=4,
                 client=None, rotate=False, pyflex_preset="default",
                 fix_windows=False, adj_src_type="cc", plot=True,
                 pyatoa_log_level="DEBUG", unit_output="VEL",
                 start_pad_s=0., end_pad_s=None, path_preprocess=None,
                 path_data=None, path_output=None, save_datasets=True,
                 save_figures=True, save_log_files=True, **kwargs):
        """
        Pyatoa preprocessing parameters that will be passed to Pyaflowa

        :type data_format: str
        :param data_format: data format for reading traces into memory. Pyatoa
            only works with 'ASCII' currently.
        :type components: str
        :param components: components to consider and tag data with. Should be
            string of letters such as 'RTZ'
        :type min_period: float
        :param min_period: Minimum filter corner in unit seconds. Bandpass
        filter if set with `max_period`, highpass filter if set without
        `max_period`, no filtering if not set and `max_period also not set
        :type max_period: float
        :param max_period: Maximum filter corner in unit seconds. Bandpass
            filter if set with `min_period`, lowpass filter if set without
            `min_period`, no filtering if not set and `min_period also not set
        :type filter_corners: int
        :param filter_corners: number of filter corners applied to filtering
        :type client: str
        :param client: Client name for ObsPy FDSN data gathering. Pyatoa will
            attempt to collect waveform and metadata based on network and
            station codes provided in the SPECFEM STATIONS file. If set None,
            no FDSN gathering will be attempted
        :type rotate: bool
        :param rotate: Attempt to rotate waveform components from NEZ -> RTZ
        :type pyflex_preset: str
        :param pyflex_preset: Parameter map for misfit window configuration
            defined by Pyflex. IF None, misfit and adjoint sources will be
            calculated on whole traces. For available choices, see Pyatoa docs
            page (pyatoa.rtfd.io)
        :type fix_windows: bool or str
        :param fix_windows: How to address misfit window evaluation at each
            evaluation. Options to re-use misfit windows collected during an
            inversion, available options:
            [True, False, 'ITER', 'ONCE']
            True: Re-use windows after first evaluation (i01s00);
            False: Calculate new windows each evaluation;
            'ITER': Calculate new windows at first evaluation of
            each iteration (e.g., i01s00... i02s00...
            'ONCE': Calculate new windows at first evaluation of
            the workflow, i.e., at self.par.BEGIN
        :type adj_src_type: str
        :param adj_src_type: Adjoint source type to evaluate misfit, defined by
            Pyadjoint. Currently available options: ['cc': cross-correlation,
            'mt': multitaper, 'wav': waveform']
        :type plot: bool
        :param plot: plot waveform figures and source receiver maps during
            the preprocessing stage
        :type pyatoa_log_level: str
        :param pyatoa_log_level: Log level to set Pyatoa, Pyflex, Pyadjoint.
            Available: ['null': no logging, 'warning': warnings only,
            'info': task tracking, 'debug': log all small details (recommended)]
        :type start_pad_s: int
        :param start_pad_s: seconds BEFORE origin time to gather data. Must be
            >= T_0 specificed in SPECFEM constants.h. Positive values only
        :type end_pad_s: int
        :param end_pad_s: seconds AFTER origin time to gather data. Must be
            >= NT * DT (from SPECFEM Par_file) postive values only.
        :type unit_output: str
        :param unit_output: Data units. Must match the synthetic output of
            external solver. Available: ['DISP': displacement, 'VEL': velocity,
            'ACC': acceleration]
        :type save_datasets: bool
        :param save_datasets: periodically save the output ASDFDataSets which
            contain data, metadata and results collected during the
            preprocessing procedure
        :type save_figures: bool
        :param save_figures: periodically save the output basemaps and
            data-synthetic waveform comparison figures
        :type save_log_files: bool
        :param save_log_files: periodically save log files created by Pyatoa
        :type path_preprocess: str
        :param path_preprocess: scratch path for preprocessing related steps
        :type path_data: str
        :param path_data: optional path for preprocessing module to discover
            waveform and meta-data.
        """
        # Shared SeisFlows parameters
        self.data_format = data_format
        self.components = components
        self.ntask = ntask
        self.nproc = nproc

        # Pyatoa specific parameters
        self.min_period = min_period
        self.max_period = max_period
        self.filter_corners = filter_corners
        self.client = client
        self.rotate = rotate
        self.pyflex_preset = pyflex_preset
        self.fix_windows = fix_windows
        self.adj_src_type = adj_src_type
        self.plot = plot
        self.pyatoa_log_level = pyatoa_log_level
        self.unit_output = unit_output
        self.start_pad_s = start_pad_s
        self.end_pad_s = end_pad_s

        self.path = path_preprocess or \
                    os.path.join(os.getcwd(), "scratch", "preprocess")
        self.path_output = path_output or os.path.join(os.getcwd(), "output")
        self.path_data = path_data

        self.save_datasets = save_datasets
        self.save_figures = save_figures
        self.save_log_files = save_log_files

    def check(self):
        """ 
        Checks Parameter and Path files, will be run at the start of a Seisflows
        workflow to ensure that things are set appropriately.
        """
        assert(self.data_format.upper() == "ASCII"), \
            "Pyatoa preprocess requires `data_format`=='ASCII'"

    def setup(self):
        """
        Sets up data preprocessing machinery by establishing an internally
        defined directory structure that will be used to store the outputs 
        of the preprocessing workflow
        """
        unix.mkdir(self.path)

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
        pyaflowa = Pyaflowa(sfpar=self.par, sfpath=self.path)
        unix.cd(pyaflowa.paths.datasets)

        # Generate the Inspector from existing datasets and save to disk
        # Allow this is fail, which might happen if we don't have enough data
        # or the Dataset is not formatted as expected
        insp = Inspector(self.par.TITLE, verbose=False)  # !!! TODO
        try:
            insp.discover()
            insp.save()
        except Exception as e:
            logger.warning(f"Uncontrolled exception in Inspector creation "
                                f"will not create inspector:\n{e}")

        # Make the final PDF for easier User ingestion of waveform/map figures
        pyaflowa.make_evaluation_composite_pdf()

        # Move scratch/ directory results into more permanent storage
        if self.save_datasets:
            datasets = glob(os.path.join(pyaflowa.paths.datasets, "*.h5"))
            self._save_quantity(datasets, tag="datasets")

        if self.save_figures:
            figures = glob(os.path.join(pyaflowa.paths.figures, "*.pdf"))
            self._save_quantity(figures, tag="figures")

        if self.save_log_files:
            logs = glob(os.path.join(pyaflowa.paths.logs, "*.txt"))
            path_out = os.path.join(self.path_output, CFGPATHS.LOGDIR)
            self._save_quantity(logs, path_out=path_out)

    def prepare_eval_grad(self, cwd, taskid, source_name, **kwargs):
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
            logger.debug("preparing files for gradient evaluation with "
                              "Pyaflowa")

        # Process all the stations for a given event using Pyaflowa
        pyaflowa = self._setup_event_pyaflowa(source_name)
        scaled_misfit = pyaflowa.process(nproc=self.nproc)

        if scaled_misfit is None:
            print(msg.cli(f"Event {source_name} returned no misfit, you may "
                          f"want to check logs and waveform figures, "
                          f"or consider discarding this event from your "
                          f"workflow", 
                          items=[pyaflowa.paths.logs, pyaflowa.paths.figures],
                          header="pyatoa preprocessing error", border="="))
            sys.exit(-1)

        # Event misfit defined by Tape et al. (2010) written to solver dir.
        self._write_residuals(path=cwd, scaled_misfit=scaled_misfit)

    def _setup_event_pyaflowa(self, source_name, iteration, step_count=""):
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
        # Outsource data processing to an event-specfic Pyaflowa instance
        pyaflowa = Pyaflowa(sfpar=self.par, sfpath=self.path)
        pyaflowa.setup(source_name=source_name, iteration=iteration, 
                       step_count=step_count, loc="*", cha="*")
        
        return pyaflowa
    
    def _save_quantity(self, filepaths, tag="", path_out=""):
        """
        Repeatable convenience function to save quantities from the scratch/
        directory to the output/ directory

        :type filepaths: list
        :param filepaths: full path to files that should be saved to output/
        :type tag: str  
        :param tag: tag for saving the files in self.path.OUTPUT. If not given, will
            save directly into the output/ directory
        :type path_out: str
        :param path_out: overwrite the default output path file naming
        """       
        if not path_out:
            path_out = os.path.join(self.path_output, tag)

        if not os.path.exists(path_out):
            unix.mkdir(path_out)

        for src in filepaths:
            dst = os.path.join(path_out, os.path.basename(src))
            unix.cp(src, dst) 

    @staticmethod
    def _write_residuals(path, scaled_misfit):
        """
        Computes residuals and saves them to a text file in the appropriate path

        :type path: str        
        :param path: scratch directory path, e.g. self.path.GRAD or self.path.FUNC
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
        if len(files) != self.ntask:
            print(msg.cli(f"Pyatoa preprocessing module did not recover the "
                          f"correct number of residual files "
                          f"({len(files)}/{self.ntask}). Please check that "
                          f"the preprocessing logs", header="error")
                  )
            sys.exit(-1)

        total_misfit = 0
        for filename in files:
            total_misfit += np.sum(np.loadtxt(filename))

        total_misfit /= self.ntask

        return total_misfit

