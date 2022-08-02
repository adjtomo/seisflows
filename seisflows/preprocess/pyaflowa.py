#!/usr/bin/env python3
"""
The Pyaflowa preprocessing module for waveform gathering, preprocessing and
misfit quantification. We use the name 'Pyaflowa' to avoid any potential
name overlaps with the actual pyatoa package.
"""
import os
import logging
import time
import random
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from glob import glob
from pyasdf import ASDFDataSet
from pyatoa import Config, Manager, Inspector, ManagerError
from pyatoa.utils.read import read_station_codes
from pyatoa.utils.images import imgs_to_pdf, merge_pdfs

from seisflows import logger
from seisflows.tools import unix
from seisflows.tools.config import Null, Dict, get_task_id
from seisflows.tools.specfem import check_source_names


class Pyaflowa:
    """
    [preprocess.pyaflowa] preprocessing and misfit quantification using Pyatoa

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
    :type start_pad_s: float
    :param start_pad_s: seconds BEFORE origin time to gather data. Must be
        >= T_0 specificed in SPECFEM constants.h. Positive values only
    :type end_pad_s: int
    :param end_pad_s: seconds AFTER origin time to gather data. Must be
        >= NT * DT (from SPECFEM Par_file) postive values only.
    :type unit_output: str
    :param unit_output: Data units. Must match the synthetic output of
        external solver. Available: ['DISP': displacement, 'VEL': velocity,
        'ACC': acceleration]
    :type export_datasets: bool
    :param export_datasets: periodically save the output ASDFDataSets which
        contain data, metadata and results collected during the
        preprocessing procedure
    :type export_figures: bool
    :param export_figures: periodically save the output basemaps and
        data-synthetic waveform comparison figures
    :type export_log_files: bool
    :param export_log_files: periodically save log files created by Pyatoa
    :type path_preprocess: str
    :param path_preprocess: scratch path for preprocessing related steps
    :type path_data: str
    :param path_data: optional path for preprocessing module to discover
        waveform and meta-data.
    """
    def __init__(self, min_period=1., max_period=10., filter_corners=4,
                 client=None, rotate=False, pyflex_preset="default",
                 fix_windows=False, adj_src_type="cc", plot=True,
                 pyatoa_log_level="DEBUG", unit_output="VEL", start_pad_s=0.,
                 end_pad_s=None, workdir=os.getcwd(), path_preprocess=None,
                 path_solver=None, path_specfem_data=None, path_data=None,
                 path_output=None, export_datasets=True, export_figures=True,
                 export_log_files=True, data_format="ascii",
                 data_case="data", components=None,
                 start=None, ntask=1, nproc=1, source_prefix=None,
                 **kwargs):
        """Pyatoa preprocessing parameters"""
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

        self.path = Dict(
            scratch=path_preprocess or os.path.join(workdir, "scratch",
                                                    "preprocess"),
            solver=path_solver or os.path.join(workdir, "scratch", "solver"),
            output=path_output or os.path.join(workdir, "output"),
            specfem_data=path_specfem_data,
            data=path_data,
        )

        # How to handle saving output data to disk
        self.export_datasets = export_datasets
        self.export_figures = export_figures
        self.export_log_files = export_log_files

        # Pyatoa-specific internal path structure for storing data etc.
        self.path["_logs"] = os.path.join(self.path.scratch, "logs")
        self.path["_tmplogs"] = os.path.join(self.path._logs, "tmp")
        self.path["_datasets"] = os.path.join(self.path.scratch, "datasets")
        self.path["_figures"] = os.path.join(self.path.scratch, "figures")

        # Where to look for externally stored waveform data and response files
        if self.path.data:
            self.path["_waveforms"] = os.path.join(self.path.data, "mseed")
            self.path["_responses"] = os.path.join(self.path.data, "seed")
        else:
            self.path["_waveforms"] = None
            self.path["_responses"] = None

        # SeisFlows parameters that should be set by other modules. Keep hidden
        # so `seisflows configure` doesn't attribute these to preprocess.
        self._data_format = data_format.upper()
        self._data_case = data_case.lower()
        self._components = components
        self._start = start
        self._ntask = ntask
        self._nproc = nproc
        self._source_prefix = source_prefix

        # Internal parameters to check against user-set parameters
        self._acceptable_data_formats = ["ASCII"]
        self._acceptable_source_prefixes = ["SOURCE", "FORCESOLUTION",
                                            "CMTSOLUTION"]

        # Internal attributes to be filled in by setup()
        self._config = None
        self._fix_windows = False
        self._station_codes = []
        self._source_names = []

    def check(self):
        """ 
        Checks Parameter and Path files, will be run at the start of a Seisflows
        workflow to ensure that things are set appropriately.
        """
        assert(self._data_format.upper() == "ASCII"), \
            "Pyatoa preprocess requires `data_format`=='ASCII'"

        assert(self.path.specfem_data is not None and
               os.path.exists(self.path.specfem_data)), (
            f"Pyatoa requires `path_specfem_data` to exist"
        )

        assert(os.path.exists(os.path.join(self.path.specfem_data,
                                           "STATIONS"))), \
            f"Pyatoa preprocessing requires that the `STATIONS` file exists " \
            f"within `path_specfem_data`"

        assert(self._source_prefix in self._acceptable_source_prefixes), (
            f"Pyatoa can only accept `source_prefix` in " 
            f"{self._acceptable_source_prefixes}, not '{self._source_prefix}'"
        )

    def setup(self):
        """
        Sets up data preprocessing machinery by establishing an internally
        defined directory structure that will be used to store the outputs 
        of the preprocessing workflow
        """
        for pathname in ["scratch", "_logs", "_tmplogs", "_datasets",
                         "_figures"]:
            unix.mkdir(self.path[pathname])

        # Generalized Config object that can be shared among all child processes
        # Contains paths to look for data and metadata
        self._config = Config(
            min_period=self.min_period, max_period=self.max_period,
            filter_corners=self.filter_corners, client=self.client,
            rotate=self.rotate, pyflex_preset=self.pyflex_preset,
            fix_windows=self.fix_windows, adj_src_type=self.adj_src_type,
            log_level=self.pyatoa_log_level, unit_output=self.unit_output,
            start_pad_s=self.start_pad_s, end_pad_s=self.end_pad_s,
            component_list=list(self._components),
            synthetics_only=bool(self._data_case == "synthetic"),
            paths={"waveforms": self.path["_waveforms"] or [],
                   "responses": self.path["_responses"] or [],
                   "events": [self.path.specfem_data]
                   }
        )
        # Generate a list of station codes that will be used to search for data
        self._station_codes = read_station_codes(
            path_to_stations=os.path.join(self.path.specfem_data, "STATIONS"),
            loc="*", cha="*"
        )
        # Get an internal list of source names. Will be the same as solver
        self._source_names = check_source_names(
            path_specfem_data=self.path.specfem_data,
            source_prefix=self._source_prefix, ntask=self._ntask
        )

    def quantify_misfit(self, source_name=None, save_residuals=None,
                        save_adjsrcs=None, iteration=1, step_count=0,
                        **kwargs):
        """
        Prepares solver for gradient evaluation by evaluating data-synthetic
        misfit and writing residuals and adjoint traces. Meant to be called by
        `workflow.evaluate_objective_function`.

        .. note::
            meant to be run on system using system.run() with access to solver

        :type source_name: str
        :param source_name: name of the event to quantify misfit for. If not
            given, will attempt to gather event id from the given task id which
            is assigned by system.run()
        :type save_residuals: str
        :param save_residuals: if not None, path to write misfit/residuls to
        :type save_adjsrcs: str
        :param save_adjsrcs: if not None, path to write adjoint sources to
        :type iteration: int
        :param iteration: current iteration of the workflow, information should
            be provided by `workflow` module if we are running an inversion.
            Defaults to 1 if not given (1st iteration)
        :type step_count: int
        :param step_count: current step count of the line search. Information
            should be provided by the `optimize` module if we are running an
            inversion. Defaults to 0 if not given (1st evaluation)
        """
        # Set the individual Config class for our given event and evaluation
        config = self._config.copy()
        config.event_id = source_name or self._source_names[get_task_id()]
        config.iteration = iteration
        config.step_count = step_count

        # Force the Manager to look in the solver directory for data
        # note: we are assuming the SeisFlows `solver` directory structure here.
        #   If we change how the default `solver` directory is named (defined by
        #   `solver.initialize_solver_directories()`), then this will break
        config.paths["waveforms"].append(
            os.path.join(self.path.solver, source_name, "traces", "obs")
        )
        config.paths["synthetics"].append(
            os.path.join(self.path.solver, source_name, "traces", "syn")
        )

        # Run misfit quantification for each station concurrently
        misfit, nwin = 0, 0
        with ProcessPoolExecutor(max_workers=unix.nproc() - 1) as executor:
            futures = [
                executor.submit(
                    self._quantify_misfit_station, config, code, save_adjsrcs)
                for code in self._station_codes
            ]
            # We only need to return misfit information. All data/results are
            # saved to the ASDFDataSet and status is logged to separate log file
            for future in futures:
                _misfit, _nwin = future.result()
                if _misfit is not None:
                    misfit += _misfit
                    nwin += _nwin

        # Calculate misfit based on the raw misfit and total number of windows
        if save_residuals:
            # Calculate the misfit based on the number of windows. Equation from
            # Tape et al. (2010). If no windows, misfit is simply raw misfit
            try:
                residuals = 0.5 * misfit / nwin
            except ZeroDivisionError:
                # Dealing with the case where nwin==0 (signifying either no
                # windows found, or calc'ing misfit on whole trace)
                residuals = misfit
            with open(save_residuals, "a") as f:
                f.write(f"{residuals:.2E}\n")

        # Combine all the individual .png files created into a single PDF
        if self.plot:
            output_fid = os.path.join(self.path._figures, self._ftag(config))
            self._make_event_figure_pdf(source_name, output_fid)

        # Finally, collect all the temporary log files and write a main log file
        pyatoa_logger = self._config_pyatoa_logger(
            fid=os.path.join(self.path._logs, f"{self._ftag(config)}.log")
        )
        pyatoa_logger.info(
            f"\n{'=' * 80}\n{'SUMMARY':^80}\n{'=' * 80}\n"
            f"SOURCE NAME: {config.event_id}\n"
            f"WINDOWS: {nwin}\n"
            f"RAW MISFIT: {misfit:.4f}\n"
            f"\n{'=' * 80}\n{'RAW LOGS':^80}\n{'=' * 80}"
            )
        self._collect_tmp_log_files(pyatoa_logger, config.event_id)

    @staticmethod
    def _ftag(config):
        """
        Create a re-usable file tag from the Config object as multiple functions
        will use this tag for file naming and file discovery.

        :type config: pyatoa.core.config.Config
        :param config: Configuration object that must contain the 'event_id',
            iteration and step count
        """
        return f"{config.event_id}_{config.iter_tag}_{config.step_tag}"

    def _quantify_misfit_station(self, config, station_code,
                                 save_adjsrcs=False):
        """
        Run misfit quantification for a single event-station pair. Gathers,
        preprocesses, windows and measures data, saves adjoint source if
        requested, and then returns the total misfit and the collected
        windows for the station.

        :type config: pyatoa.core.config.Config
        :param config: Config object that defines all the processing parameters
            required by the Pyatoa workflow
        :type station_code: str
        :param station_code: chosen station to quantify misfit for. Should be
            in the format 'NN.SSS.LL.CCC'
        :type save_adjsrcs: str
        :param save_adjsrcs: path to directory where adjoint sources should be
            saved. Filenames will be generated automatically by Pyatoa to fit
            the naming schema required by SPECFEM. If False, no adjoint sources
            will be saved. They of course can be saved manually later using
            Pyatoa + PyASDF
        """
        # Unique identifier for the given source-receiver pair for file naming
        # Something like 001_i01_s00_XX_XYZ
        net, sta, loc, cha = station_code.split(".")
        tag = f"{self._ftag(config)}_{net}_{sta}"

        # Configure a single source-receiver pair logger which will be collected
        # later by the main function
        log_file = os.path.join(self.path._tmplogs, f"{tag}.log")
        station_logger = self._config_pyatoa_logger(fid=log_file)
        station_logger.info(f"\n{'/' * 80}\n{station_code:^80}\n{'/' * 80}")

        # Begin data gathering/processing and misfit quantification
        mgmt = Manager(config=config)
        mgmt.gather(choice=["event"], event_id=config.event_id,
                    prefix=f"{self._source_prefix}_")

        # Attempt to gather data. If fail, return because theres nothing else
        # we can do without data
        _processed = False
        try:
            mgmt.gather(choice=["inv", "st_obs", "st_syn"], code=station_code)
        except ManagerError as e:
            station_logger.critical(e)
            return None, None

        # If any part of the processing fails, move on to plotting because we
        # will have gathered waveform data so a figure is still useful.
        try:
            _fix_windows = self._check_fixed_windows(
                iteration=mgmt.config.iteration,
                step_count=mgmt.config.step_count,
                logger=station_logger,
            )
            mgmt.standardize()
            mgmt.preprocess()
            mgmt.window(fix_windows=_fix_windows)
            mgmt.measure()
            _processed = True
        except ManagerError as e:
            station_logger.warning(e)
            pass

        # Plot waveform + map figure. Map may fail if we don't have appropriate
        # metdata, in which case we fall back to plotting waveform only
        if self.plot:
            # e.g., 001_i01_s00_XX_ABC.png
            save = os.path.join(self.path["_figures"], f"{tag}.png")
            try:
                mgmt.plot(choice="both", show=False, save=save)
            except ManagerError as e:
                station_logger.warning(e)
                mgmt.plot(choice="wav", show=False, save=save)

        # Write out the .adj adjoint source files for solver to discover.
        # Write empty adjoint sources for components with no adjoint sources
        if _processed and save_adjsrcs:
            mgmt.write_adjsrcs(path=save_adjsrcs, write_blanks=True)

        # Wait until the very end to write to the HDF5 file, then do it
        # pseudo-serially to get around trying to parallel write to HDF5 file
        while True:
            try:
                with ASDFDataSet(os.path.join(self.path["_datasets"],
                                              f"{config.event_id}.h5")) as ds:
                    mgmt.write(ds=ds)
                break
            except (BlockingIOError, FileExistsError):
                # Random sleep time [0,1]s to decrease chances of two processes
                # attempting to access at exactly the same time
                time.sleep(random.random())

        return mgmt.stats.misfit, mgmt.stats.nwin

    def sum_residuals(self, residuals):
        """
        Return summed residuals devided by number of events following equation
        in Tape et al. 2010

        :type residuals: np.array
        :param residuals: list of residuals from each NTASK event
        :rtype: float
        :return: sum of squares of residuals
        """
        assert(len(residuals) == self._ntask), \
            f"recovered an incorrect number of residual values"
        summed_residuals = np.sum(residuals)
        return summed_residuals / self._ntask

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
        # Generate the Inspector from existing datasets and save to disk
        # Allow this is fail, which might happen if we don't have enough data
        # or the Dataset is not formatted as expected
        unix.cd(self.path._datasets)
        insp = Inspector("inspector", verbose=False)
        try:
            insp.discover()
            insp.save()
        except Exception as e:
            logger.warning(f"Uncontrolled exception in Pyatoa Inspector "
                           f"creation -- will not create inspector:\n{e}")

        # Make the final PDF for easier User ingestion of waveform/map figures
        self._make_evaluation_composite_pdf()

        # Move scratch/ directory results into more permanent storage
        if self.export_datasets:
            src = glob(os.path.join(self.path._datasets, "*.h5"))
            dst = os.path.join(self.path.output, "datasets", "")
            unix.mkdir(dst)
            unix.cp(src, dst)

        if self.export_figures:
            src = glob(os.path.join(self.path._figures, "*.pdf"))
            dst = os.path.join(self.path.output, "figures", "")
            unix.mkdir(dst)
            unix.cp(src, dst)

        if self.export_log_files:
            src = glob(os.path.join(self.path._logs, "*.txt"))
            dst = os.path.join(self.path.output, "logs", "")
            unix.mkdir(dst)
            unix.cp(src, dst)

    def _check_fixed_windows(self, iteration, step_count, logger=Null()):
        """
        Determine how to address re-using misfit windows during an inversion
        workflow. Throw some log messages out to let the User know whether or
        not misfit windows will be re used throughout an inversion.

            True: Always fix windows except for i01s00 because we don't have any
                  windows for the first function evaluation
            False: Don't fix windows, always choose a new set of windows
            Iter: Pick windows only on the initial step count (0th) for each
                  iteration. WARNING - does not work well with Thrifty Inversion
                  because the 0th step count is usually skipped
            Once: Pick new windows on the first function evaluation and then fix
                  windows. Useful for when parameters have changed, e.g. filter
                  bounds

        :type iteration: int
        :param iteration: The current iteration of the SeisFlows3 workflow,
            within SeisFlows3 this is defined by `optimize.iter`
        :type step_count: int
        :param step_count: Current line search step count within the SeisFlows3
            workflow. Within SeisFlows3 this is defined by
            `optimize.line_search.step_count`
        :type logger: logging.Logger
        :param logger: The main logger for a given event, should be
            defined by `pyaflowa.quantify_misfit()`. If not provided, logs will
            get sent to DevNull
        :rtype: bool
        :return: bool on whether to use windows from the previous step
        """
        fix_windows = False
        # First function evaluation never fixes windows
        if iteration == 1 and step_count == 0:
            fix_windows = False
            logger.info("new windows; first evaluation")
        elif isinstance(self.fix_windows, str):
            # By 'iter'ation only pick new windows on the first step count
            if self.fix_windows.upper() == "ITER":
                if step_count == 0:
                    fix_windows = False
                    logger.info("new windows; first step count")
                else:
                    fix_windows = True
                    logger.info("fix windows; mid line search")
            # 'Once' picks windows only for the first function evaluation of
            # the current set of iterations.
            elif self.fix_windows.upper() == "ONCE":
                if iteration == self._start and step_count == 0:
                    fix_windows = False
                    logger.info("new windows; first workflow evaluation")
                else:
                    fix_windows = True
                    logger.info("fix windows; mid workflow")
        # Bool fix windows simply sets the parameter
        elif isinstance(self.fix_windows, bool):
            fix_windows = self.fix_windows
            logger.info(f"fixed windows flag set: {self.fix_windows}")

        return fix_windows

    def _config_pyatoa_logger(self, fid):
        """
        Create a log file to track processing of a given source-receiver pair.
        Because each station is processed asynchronously, we don't want them to
        log to the main file at the same time, otherwise we get a random mixing
        of log messages. Instead we have them log to temporary files, which
        are combined at the end of the processing script in serial.

        :type fid: str
        :param fid: full path and filename for logger that will be configured
        :rtype: logging.Logger
        :return: a logger which does NOT log to stdout and only logs to
            the given file defined by `fid`
        """
        handler = logging.FileHandler(fid, mode="w")
        logfmt = "[%(asctime)s] - %(name)s - %(levelname)s: %(message)s"
        formatter = logging.Formatter(logfmt, datefmt="%Y-%m-%d %H:%M:%S")
        handler.setFormatter(formatter)
        for log in ["pyflex", "pyadjoint", "pyatoa"]:
            # Set the overall log level
            logger = logging.getLogger(log)
            # Turn off any existing handlers (stream and file)
            while logger.hasHandlers():
                logger.removeHandler(logger.handlers[0])
            # Log to new temporary file
            logger.setLevel(self.pyatoa_log_level)
            logger.addHandler(handler)

        return logger

    def _collect_tmp_log_files(self, pyatoa_logger, event_id):
        """
        Each source-receiver pair has made its own log file. This function
        collects these files and writes their content back into the main log.
        This is a lot of IO but should be okay since the files are small.

        .. note::
            This was the most foolproof method for having multiple parallel
            processes write to the same file. I played around with StringIO
            buffers and file locks, but they became overly complicated and
            ultimately did not work how I wanted them to. This function trades
            filecount and IO overhead for simplicity.

        .. warning::
            The assumption here is that the number of source-receiver pairs
            is manageable (in the thousands). If we start reaching file count
            limits on the cluster then this method for logging may have to be
            re-thought. See link for example:
            https://stackless.readthedocs.io/en/3.7-slp/howto/
              logging-cookbook.html#using-concurrent-futures-processpoolexecutor

        :type pyatoa_logger: logging.Logger
        :param pyatoa_logger: The main logger for a given event, should be
            defined by `pyaflowa.quantify_misfit()`
        :type event_id: str
        :param event_id: given event id that we are concerned with. Used to
            search for matching log files in the temporary log file directory
        """
        tmp_logs = sorted(glob(os.path.join(self.path._tmplogs,
                                            f"*{event_id}_*.log")))
        with open(pyatoa_logger.handlers[0].baseFilename, "a") as fw:
            for tmp_log in tmp_logs:
                with open(tmp_log, "r") as fr:
                    fw.writelines(fr.readlines())
                unix.rm(tmp_log)  # delete after writing

    def _make_event_figure_pdf(self, source_name, output_fid):
        """
        Combine a list of single source-receiver PNGs into a single PDF file
        for the given event. Mostly a convenience function to make it easier
        to ingest waveform figures during a workflow.

        """
        # Sorrted by network and station name
        input_fids = sorted(glob(os.path.join(self.path._figures,
                                              f"{source_name}*.png")))
        if not input_fids:
            logger.warning(f"Pyatoa found no event figures for {source_name} "
                           f"to combine")
            return

        # Assuming the file naming format defined by `_quantify_misfit_station`
        source_name, iteration, step_count, *_ = input_fids[0].split("_")

        # e.g. i01s00_2018p130600.pdf
        output_fid = (f"{source_name}_{iteration}_{step_count}.pdf")

        # Merge all output pdfs into a single pdf, delete originals
        save = os.path.join(self.path._figures, output_fid)
        imgs_to_pdf(fids=sorted(input_fids), fid_out=save)
        for fid in input_fids:
            os.remove(fid)

    def _make_evaluation_composite_pdf(self):
        """
        Utility function to combine all PDFs generated by
        make_event_figure_pdf() into a single PDF tagged by the current
        evaluation (iteration, step count). This is meant to make it easier
        for the User to scroll through figures. Option to delete the original
        event-specific PDFs which are now redundant
        .. note::
            This can be run without running Pyaflowa.format()
        :type delete_originals: bool
        :param delete_originals: delete original pdf files after mergin
        """
        event_figures = glob(os.path.join(self.path._figures, "*", "*.pdf"))
        if not event_figures:
            logger.warning("Pyatoa could not find event PDFs to merge")
            return
        # Collecting evaluation tags, e.g., ['i01s00', 'i01s01']
        tags = set([os.path.basename(_).split("_")[0] for _ in event_figures])
        for tag in tags:
            fids = [fid for fid in event_figures if tag in fid]
            fid_out = os.path.join(self.path._figures, f"{tag}.pdf")
            merge_pdfs(fids=sorted(fids), fid_out=fid_out)
            for fid in fids:
                os.remove(fid)
