#!/usr/bin/env python3
"""
The Pyaflowa preprocessing module for waveform gathering, preprocessing and
misfit quantification. We use the name 'Pyaflowa' to avoid any potential
name overlaps with the actual pyatoa package.
"""
import os
import numpy as np
from pyasdf import ASDFDataSet
from pyatoa import Config, Manager, ManagerError
from pyatoa.utils.read import read_station_codes

from seisflows import logger
from seisflows.tools import unix
from seisflows.tools.config import Dict, get_task_id
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
        for pathname in ["scratch", "_logs", "_datasets", "_figures"]:
            unix.mkdir(self.path[pathname])

        # Generalized Config object that can be shared among all child processes
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

    def quantify_misfit(self, source_name, save_residuals=None,
                        save_adjsrcs=None, iteration=1, step_count=0,
                        **kwargs):
        """
        Prepares solver for gradient evaluation by writing residuals and
        adjoint traces. Meant to be called by
        `workflow.evaluate_objective_function`

        Reads in observed and synthetic waveforms, applies optional
        preprocessing, assesses misfit, and writes out adjoint sources and
        STATIONS_ADJOINT file.

        .. note::
            Meant to be called by solver.eval_func(), may have unused arguments
            to keep functions general across subclasses.

        :type save_residuals: str
        :param save_residuals: if not None, path to write misfit/residuls to
        :type save_adjsrcs: str
        :param save_adjsrcs: if not None, path to write adjoint sources to
        """
        # Set the individual Config class for our given event and evaluation
        config = self._config.copy()
        config.event_id = source_name
        config.iteration = iteration
        config.step_count = step_count

        # Force the Manager to look in the solver directory for data
        # Note: we are assuming the SeisFlows solver directory structure here
        config.paths["waveforms"].append(
            os.path.join(self.path.solver, source_name, "traces", "obs")
        )
        config.paths["synthetics"].append(
            os.path.join(self.path.solver, source_name, "traces", "syn")
        )

        # Set up the Pyatoa workflow, attempt to gather event metadata
        ds = ASDFDataSet(
            os.path.join(self.path["_datasets"], f"{source_name}.h5")
        )
        mgmt = Manager(config=config, ds=ds)
        mgmt.gather(choice=["event"], event_id=source_name,
                    prefix=f"{self._source_prefix}_")

        # Run data/metadata gathering, processing and misfit quantification
        misfit, nwin = 0, 0
        for station_code in self._station_codes:
            net, sta, loc, cha = station_code.split(".")
            _processed = False
            # Will gather data and metadata based on the station codes and
            # input paths from the Configuration object
            try:
                mgmt.gather(choice=["inv", "st_obs", "st_syn"],
                            code=station_code)
            except ManagerError as e:
                continue
            # If any part of the processing fails, move on to plotting
            try:
                mgmt.standardize()
                mgmt.preprocess()
                mgmt.window(
                    fix_windows=self._check_fixed_windows(iteration, step_count)
                )
                mgmt.measure()
                _processed = True
            except ManagerError as e:
                pass

            if self.plot:
                plot_fid = (
                    f"{source_name}_{config.iter_tag}_{config.step_tag}_"
                    f"{net}_{sta}.png"
                )
                save = os.path.join(self.path["_figures"], plot_fid)
                try:
                    mgmt.plot(choice="both", show=False, save=save)
                except ManagerError:
                    mgmt.plot(choice="wav", show=False, save=save)

            # Write out the .adj adjoint source files
            if _processed and save_adjsrcs:
                mgmt.write_adjsrcs(path=save_adjsrcs, write_blanks=True)

            misfit += mgmt.stats.misfit
            nwin += mgmt.stats.nwin

        if save_residuals:
            try:
                residuals = 0.5 * misfit / nwin
            except ZeroDivisionError:
                # Dealing with the case where nwin==0 (signifying either no
                # windows found, or calc'ing misfit on whole trace)
                residuals = misfit
            np.savetxt(save_residuals, [residuals], fmt="%11.6e")

    def _check_fixed_windows(self, iteration, step_count):
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

    #
    # def finalize(self):
    #     """
    #     Run some serial finalization tasks specific to Pyatoa, which will help
    #     aggregate the collection of output information.
    #
    #     .. note::
    #         This finalize function performs the following tasks:
    #         * Generate .csv files using the Inspector
    #         * Aggregate event-specific PDFs into a single evaluation PDF
    #         * Save scratch/ data into output/ if requested
    #     """
    #     # Initiate Pyaflowa to get access to path structure
    #     pyaflowa = Pyaflowa(sfpar=self.par, sfpath=self.path)
    #     unix.cd(pyaflowa.paths.datasets)
    #
    #     # Generate the Inspector from existing datasets and save to disk
    #     # Allow this is fail, which might happen if we don't have enough data
    #     # or the Dataset is not formatted as expected
    #     insp = Inspector(self.par.TITLE, verbose=False)  # !!! TODO
    #     try:
    #         insp.discover()
    #         insp.save()
    #     except Exception as e:
    #         logger.warning(f"Uncontrolled exception in Inspector creation "
    #                             f"will not create inspector:\n{e}")
    #
    #     # Make the final PDF for easier User ingestion of waveform/map figures
    #     pyaflowa.make_evaluation_composite_pdf()
    #
    #     # Move scratch/ directory results into more permanent storage
    #     if self.export_datasets:
    #         datasets = glob(os.path.join(pyaflowa.paths.datasets, "*.h5"))
    #         self._save_quantity(datasets, tag="datasets")
    #
    #     if self.export_figures:
    #         figures = glob(os.path.join(pyaflowa.paths.figures, "*.pdf"))
    #         self._save_quantity(figures, tag="figures")
    #
    #     if self.export_log_files:
    #         logs = glob(os.path.join(pyaflowa.paths.logs, "*.txt"))
    #         path_out = os.path.join(self.path_output, CFGPATHS.LOGDIR)
    #         self._save_quantity(logs, path_out=path_out)
    #
    # def prepare_eval_grad(self, cwd, taskid, source_name, **kwargs):
    #     """
    #     Prepare the gradient evaluation by gathering, preprocessing waveforms,
    #     and measuring misfit between observations and synthetics using Pyatoa.
    #
    #     Reads in observed and synthetic waveforms, applies optional
    #     preprocessing, assesses misfit, and writes out adjoint sources and
    #     STATIONS_ADJOINT file.
    #
    #     .. note::
    #         Meant to be called by solver.eval_func(), may have unused arguments
    #         to keep functions general across preprocessing subclasses.
    #
    #     :type cwd: str
    #     :param cwd: current specfem working directory containing observed and
    #         synthetic seismic data to be read and processed. Should be defined
    #         by solver.cwd
    #     :type source_name: str
    #     :param source_name: the event id to be used for tagging and data lookup.
    #         Should be defined by solver.source_name
    #     :type taskid: int
    #     :param taskid: identifier of the currently running solver instance.
    #         Should be defined by solver.taskid
    #     :type filenames: list of str
    #     :param filenames: [not used] list of filenames defining the files in
    #         traces
    #     """
    #     if taskid == 0:
    #         logger.debug("preparing files for gradient evaluation with "
    #                           "Pyaflowa")
    #
    #     # Process all the stations for a given event using Pyaflowa
    #     pyaflowa = self._setup_event_pyaflowa(source_name)
    #     scaled_misfit = pyaflowa.process(nproc=self.nproc)
    #
    #     if scaled_misfit is None:
    #         print(msg.cli(f"Event {source_name} returned no misfit, you may "
    #                       f"want to check logs and waveform figures, "
    #                       f"or consider discarding this event from your "
    #                       f"workflow",
    #                       items=[pyaflowa.paths.logs, pyaflowa.paths.figures],
    #                       header="pyatoa preprocessing error", border="="))
    #         sys.exit(-1)
    #
    #     # Event misfit defined by Tape et al. (2010) written to solver dir.
    #     self._write_residuals(path=cwd, scaled_misfit=scaled_misfit)
    #
    # def _setup_event_pyaflowa(self, source_name, iteration, step_count=""):
    #     """
    #     A convenience function to set up a Pyaflowa processing instance for
    #     a specific event.
    #
    #     .. note::
    #         This is meant to be called by preprocess.prepare_eval_grad() but its
    #         also useful for debugging and manual processing where you can simply
    #         return a formatted Pyaflowa object and debug it directly.
    #
    #     :type source_name: str
    #     :param source_name: solver source name to evaluate setup for. Must
    #         match from list defined by: solver.source_names
    #     """
    #     # Outsource data processing to an event-specfic Pyaflowa instance
    #     pyaflowa = Pyaflowa(sfpar=self.par, sfpath=self.path)
    #     pyaflowa.setup(source_name=source_name, iteration=iteration,
    #                    step_count=step_count, loc="*", cha="*")
    #
    #     return pyaflowa
    #
    # def _save_quantity(self, filepaths, tag="", path_out=""):
    #     """
    #     Repeatable convenience function to save quantities from the scratch/
    #     directory to the output/ directory
    #
    #     :type filepaths: list
    #     :param filepaths: full path to files that should be saved to output/
    #     :type tag: str
    #     :param tag: tag for saving the files in self.path.OUTPUT. If not given, will
    #         save directly into the output/ directory
    #     :type path_out: str
    #     :param path_out: overwrite the default output path file naming
    #     """
    #     if not path_out:
    #         path_out = os.path.join(self.path_output, tag)
    #
    #     if not os.path.exists(path_out):
    #         unix.mkdir(path_out)
    #
    #     for src in filepaths:
    #         dst = os.path.join(path_out, os.path.basename(src))
    #         unix.cp(src, dst)
    #
    # @staticmethod
    # def _write_residuals(path, scaled_misfit):
    #     """
    #     Computes residuals and saves them to a text file in the appropriate path
    #
    #     :type path: str
    #     :param path: scratch directory path, e.g. self.path.GRAD or self.path.FUNC
    #     :type scaled_misfit: float
    #     :param scaled_misfit: the summation of misfit from each
    #         source-receiver pair calculated by prepare_eval_grad()
    #     :type source_name: str
    #     :param source_name: name of the source related to the misfit, used
    #         for file naming
    #     """
    #     residuals_file = os.path.join(path, "residuals")
    #     np.savetxt(residuals_file, [scaled_misfit], fmt="%11.6e")
    #
    # def sum_residuals(self, files):
    #     """
    #     Averages the event misfits and returns the total misfit.
    #     Total misfit defined by Tape et al. (2010)
    #
    #     :type files: str
    #     :param files: list of single-column text files containing residuals
    #         that will have been generated using prepare_eval_grad()
    #     :rtype: float
    #     :return: average misfit
    #     """
    #     if len(files) != self.ntask:
    #         print(msg.cli(f"Pyatoa preprocessing module did not recover the "
    #                       f"correct number of residual files "
    #                       f"({len(files)}/{self.ntask}). Please check that "
    #                       f"the preprocessing logs", header="error")
    #               )
    #         sys.exit(-1)
    #
    #     total_misfit = 0
    #     for filename in files:
    #         total_misfit += np.sum(np.loadtxt(filename))
    #
    #     total_misfit /= self.ntask
    #
    #     return total_misfit
    #
