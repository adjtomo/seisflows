#!/usr/bin/env python3
"""
The Pyaflowa preprocessing module for waveform gathering, preprocessing and
misfit quantification. We use the name 'Pyaflowa' to avoid any potential
name overlaps with the actual Pyatoa package.
"""
import os
import logging
import time
import random
import numpy as np
from concurrent.futures import ProcessPoolExecutor, wait
from glob import glob
from pyasdf import ASDFDataSet

from pyatoa import Config, Manager, Inspector, ManagerError
from pyatoa.utils.read import read_station_codes, read_events_plus

from seisflows import logger
from seisflows.tools import unix
from seisflows.tools.config import Dict
from seisflows.tools.graphics import imgs_to_pdf, merge_pdfs
from seisflows.tools.specfem import (check_source_names,
                                     return_matching_waveform_files,
                                     get_src_rcv_lookup_table)
from seisflows.preprocess.default import read, initialize_adjoint_traces


class Pyaflowa:
    """
    Pyaflowa Preprocess
    -------------------
    Preprocessing and misfit quantification using Python's Adjoint Tomography
    Operations Assistant (Pyatoa)

    Parameters
    ----------
    :type min_period: float
    :param min_period: Minimum filter corner in unit seconds. Bandpass
    filter if set with `max_period`, highpass filter if set without
    `max_period`, no filtering if not set and `max_period also not set
    :type pyflex_parameters: dict
    :param pyflex_parameters: overwrite for Pyflex parameters defined
        in the Pyflex.Config object. Incorrectly defined argument names
        will raise a TypeError. See Pyflex docs for detailed parameter defs:
        http://adjtomo.github.io/pyflex/#config-object
    :type pyadjoint_parameters: dict
    :param pyadjoint_parameters: overwrite for Pyadjoint parameters defined
        in the Pyadjoint.Config object for the given `adj_src_type`.
        Incorrectly defined argument names will raise a TypeError. See
        Pyadjoint docs for detailed parameter definitions:
        https://adjtomo.github.io/pyadjoint/
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
        Pyadjoint. See `pyadjoint.config.ADJSRC_TYPES` for detailed options list
        - 'waveform': waveform misfit function
        - 'convolution': convolution misfit function
        - 'exponentiated_phase': exponentiated phase from Yuan et al. 2020
        - 'cc_traveltime': cross-correlation traveltime misfit
        - 'multitaper': multitaper misfit function
    :type plot_waveforms: bool
    :param plot_waveforms: plot waveform figures and source receiver maps during
        the preprocessing stage. Maps require metadata, and if they are not
        provided then only waveforms + windows + adjoint sources will be plotted
    :type preprocess_log_level: str
    :param preprocess_log_level: Log level to set Pyatoa, Pyflex, Pyadjoint.
        Available: ['null': no logging, 'warning': warnings only,
        'info': task tracking, 'debug': log all small details (recommended)]
    :type unit_output: str
    :param unit_output: Data units. Must match the synthetic output of
        external solver. Available: ['DISP': displacement, 'VEL': velocity,
        'ACC': acceleration]. Requires metadata.
    :type export_datasets: bool
    :param export_datasets: periodically save the output ASDFDataSets which
        contain data, metadata and results collected during the
        preprocessing procedure
    :type export_figures: bool
    :param export_figures: periodically save the output basemaps and
        data-synthetic waveform comparison figures
    :type export_log_files: bool
    :param export_log_files: periodically save log files created by Pyatoa

    Paths
    -----
    :type path_preprocess: str
    :param path_preprocess: scratch path for preprocessing related steps
    ***
    """
    def __init__(self, min_period=1., max_period=10.,
                 pyflex_parameters=None, pyadjoint_parameters=None,
                 fix_windows=False, adj_src_type="cc_traveltime",
                 plot_waveforms=True, preprocess_log_level="DEBUG",
                 export_datasets=True, export_figures=True,
                 export_log_files=True, workdir=os.getcwd(),
                 path_preprocess=None, path_solver=None, path_specfem_data=None,
                 path_data=None, path_output=None, obs_data_format="SAC",
                 syn_data_format="ascii", data_case="data", components=None,
                 start=None, ntask=1, nproc=1, source_prefix=None, **kwargs):
        """
        Pyatoa preprocessing parameters

        .. note::
            Paths and parameters listed here are shared with other modules and
            so are not included in the main class docstring.

        :type syn_data_format: str
        :param syn_data_format: data format for reading synthetic traces into
            memory. Shared with solver module. Pyatoa only works with 'ASCII'
            currently.
        :type data_case: str
        :param data_case: How to address 'data' in the workflow, options:
            'data': real data will be provided by the user in
            `path_data/{source_name}` in the same format that the solver will
            produce synthetics (controlled by `solver.format`) OR
            synthetic': 'data' will be generated as synthetic seismograms using
            a target model provided in `path_model_true`. If None, workflow will
            not attempt to generate data.
        :type components: str
        :param components: components to search for synthetic data with. None by
            default which uses a wildcard when searching for synthetics. If
            provided, User only wants to use a subset of components generated by
            SPECFEM. In that case, `components` should be string of letters such
            as 'ZN' (for up and north components)
        :type workdir: str
        :param workdir: working directory in which to look for data and store
            results. Defaults to current working directory
        :type path_solver: str
        :param path_solver: scratch path for all solver related tasks
        :type path_data: str
        :param path_data: path to any externally stored data required by the 
            solver
        """
        # Pyatoa related parameters
        self.min_period = min_period
        self.max_period = max_period
        self.fix_windows = fix_windows
        self.adj_src_type = adj_src_type
        self.plot_waveforms = plot_waveforms
        self.preprocess_log_level = preprocess_log_level

        self.path = Dict(
            scratch=path_preprocess or os.path.join(workdir, "scratch",
                                                    "preprocess"),
            solver=path_solver or os.path.join(workdir, "scratch", "solver"),
            output=path_output or os.path.join(workdir, "output"),
            specfem_data=path_specfem_data,
            data=path_data,
        )

        # Set the Pyflex and Pyadjoint external parameters
        _cfg = Config(pyflex_parameters=pyflex_parameters,
                      pyadjoint_parameters=pyadjoint_parameters)
        self.pyflex_parameters = {
            key: val for key, val in _cfg.pfcfg.items() if key not in
            ["min_period", "max_period"]}
        # Ignore meta-config parameters that are hardcoded in lower levels
        self.pyadjoint_parameters = {
            key: val for key, val in _cfg.pacfg.items() if key not in
            ["min_period", "max_period", "adjsrc_type", "double_difference"]
        }

        # How to handle saving output data to disk
        self.export_datasets = export_datasets
        self.export_figures = export_figures
        self.export_log_files = export_log_files

        # Pyatoa-specific internal directory structure for storing data etc.
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

        # Pyaflowa parameters that should be set by other modules.
        self.syn_data_format = syn_data_format.upper()
        self.obs_data_format = obs_data_format.upper()
        self.source_prefix = source_prefix
        self._data_case = data_case.lower()
        if components is not None:
            self._components = list(components)  # e.g. 'RTZ' -> ['R', 'T', 'Z']
        else:
            self._components = components

        self._start = start
        self._ntask = ntask
        self._nproc = nproc
        self._source_prefix = source_prefix

        # Internal parameters to check against user-set parameters
        self._syn_acceptable_data_formats = ["ASCII"]
        self._acceptable_source_prefixes = ["SOURCE", "FORCESOLUTION",
                                            "CMTSOLUTION"]
        self._acceptable_fix_windows = ["ITER", "ONCE", True, False]
        self._acceptable_unit_outputs = ["VEL", "DISP", "ACC"]

        # Internal attributes to be filled in by setup()
        self._srcrcv_stats = None
        self._config = None
        self._fix_windows = False
        self._station_codes = []
        self._source_names = []

        # Turn off main logger, which will be toggled on by later jobs
        logging.getLogger("pyatoa").setLevel("CRITICAL")

    def check(self):
        """ 
        Checks Parameter and Path files, will be run at the start of a Seisflows
        workflow to ensure that things are set appropriately.
        """
        assert(self.syn_data_format.upper() == "ASCII"), \
            "Pyatoa preprocess requires `syn_data_format`=='ASCII'"

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
        assert(self._fix_windows in self._acceptable_fix_windows), \
            f"Pyatoa `fix_windows` must be in {self._acceptable_fix_windows}"

    def setup(self):
        """
        Sets up data preprocessing machinery by establishing an internally
        defined directory structure that will be used to store the outputs 
        of the preprocessing workflow
        """
        for pathname in ["scratch", "_logs", "_tmplogs", "_datasets",
                         "_figures"]:
            unix.mkdir(self.path[pathname])

        if self._data_case == "synthetic":
            st_obs_type = "syn"
        else:
            st_obs_type = "obs"

        # Convert SeisFlows user parameters into Pyatoa config parameters
        self._config = Config(
            min_period=self.min_period, max_period=self.max_period,
            adj_src_type=self.adj_src_type, component_list=self._components,
            st_obs_type=st_obs_type, st_syn_type="syn",
            pyflex_parameters=self.pyflex_parameters,
            pyadjoint_parameters=self.pyadjoint_parameters
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

        # Get a lookup table providing relationships between each source and sta
        self._srcrcv_stats = get_src_rcv_lookup_table(
            path_to_data=self.path.specfem_data,
            source_prefix=self.source_prefix
        )

    @staticmethod
    def ftag(config):
        """
        Create a re-usable file tag from the Config object as multiple functions
        will use this tag for file naming and file discovery.

        :type config: pyatoa.core.config.Config
        :param config: Configuration object that must contain the 'event_id',
            iteration and step count
        """
        return f"{config.event_id}_{config.iter_tag}_{config.step_tag}"

    def quantify_misfit(self, source_name=None, save_residuals=None,
                        export_residuals=None, save_adjsrcs=None,
                        components=None, iteration=1, step_count=0,
                        _serial=False, **kwargs):
        """
        Main processing function to be called by Workflow module. Generates
        total misfit and adjoint sources for a given event with name 
        `source_name`.

        .. note::

            Meant to be called by `workflow.evaluate_objective_function` and
            run on system using system.run() to get access to compute nodes.

        :type source_name: str
        :param source_name: name of the event to quantify misfit for. If not
            given, will attempt to gather event id from the given task id which
            is assigned by system.run()
        :type save_residuals: str
        :param save_residuals: if not None, path to write misfit/residuls to
        :type export_residuals: str
        :param export_residuals: export all residuals (data-synthetic misfit)
            that are generated by the external solver to `path_output`. If
            False, residuals stored in scratch may be discarded at any time in
            the workflow
        :type save_adjsrcs: str
        :param save_adjsrcs: if not None, path to write adjoint sources to
        :type components: list
        :param components: optional list of components to ignore preprocessing
            traces that do not have matching components. The adjoint sources for
            these components will be 0. E.g., ['Z', 'N']. If None, all available
            components will be considered.
        :type iteration: int
        :param iteration: current iteration of the workflow, information should
            be provided by `workflow` module if we are running an inversion.
            Defaults to 1 if not given (1st iteration)
        :type step_count: int
        :param step_count: current step count of the line search. Information
            should be provided by the `optimize` module if we are running an
            inversion. Defaults to 0 if not given (1st evaluation)
        :type _serial: bool
        :param _serial: debug function to turn preprocessing to a serial task
            whereas it is normally a multiprocessed parallel task
        """
        # Set a unique Config object to specify which Source we want to process
        config = self._config.copy()
        config.event_id = source_name
        config.iteration = iteration
        config.step_count = step_count
        if components is not None:
            config.component_list = components

        # Generate empty adjoint sources and return a matching list of files
        # that will be fed into the misfit quantification machinery
        obs, syn = self._setup_quantify_misfit(source_name, save_adjsrcs,
                                               components)

        # Process each pair in serial.
        if _serial:
            total_misfit, total_windows = 0, 0
            for o, s in zip(obs, syn):
                misfit, nwin = self._quantify_misfit_single(o, s, config,
                                                            save_adjsrcs)

                total_misfit += misfit or 0
                total_windows += nwin or 0
        # Process each pair in parallel. Max workers is total num. of cores
        else:
            with ProcessPoolExecutor(max_workers=unix.nproc()) as executor:
                futures = [
                    executor.submit(self._quantify_misfit_single, o, s, config,
                                    save_adjsrcs) for o, s in zip(obs, syn)
                ]
            wait(futures)

            # Initialize empty values to store statistics on entire misfit quant
            total_misfit, total_windows = 0, 0
            for future in futures:
                misfit, nwin = future.result()

                total_misfit += misfit or 0
                total_windows += nwin or 0

        # Save residuals to external file for Workflow to calculate misfit `f`
        # Slightly different than Default preprocessing because we need to
        # normalize by the total number of windows
        if save_residuals:
            with open(save_residuals, "a") as f:
                f.write(f"{total_misfit / total_windows:.2E}\n")
        if export_residuals:
            if not os.path.exists(export_residuals):
                unix.mkdir(export_residuals)
            unix.cp(src=save_residuals, dst=export_residuals)

        # Combine all the individual .png files created into a single PDF for
        # easier scrolling convenience
        if self.plot_waveforms:
            fids = sorted(glob(os.path.join(self.path._figures,
                                            f"{source_name}*.png")))
            fid_out = os.path.join(self.path._figures,
                                   f"{self.ftag(config)}.pdf")
            imgs_to_pdf(fids, fid_out, remove_fids=True)

        # Collect all temp log files into a single log file
        self._finalize_logging(config, total_windows, total_misfit)

    def _setup_quantify_misfit(self, source_name, save_adjsrcs=None,
                               components=None):
        """
        Gather a list of filenames of matching waveform IDs that can be
        run through the misfit quantification step, and generate empty adjoint
        sources so that Solver knows which components are zero'd out.

        :type source_name: str
        :param source_name: the name of the source to process
        :type components: list
        :param components: optional list of components to ignore preprocessing
            traces that do not have matching components. The adjoint sources for
            these components will be 0. E.g., ['Z', 'N']. If None, all available
            components will be considered.
        :rtype: list of tuples
        :return: [(observed filename, synthetic filename)]. tuples will contain
            filenames for matching stations + component for obs and syn
        """
        obs_path = os.path.join(self.path.solver, source_name, "traces", "obs")
        syn_path = os.path.join(self.path.solver, source_name, "traces", "syn")

        # Initialize empty adjoint sources for all synthetics that may or may
        # not be overwritten by the misfit quantification step
        if save_adjsrcs is not None:
            syn_filenames = glob(os.path.join(syn_path, "*"))
            initialize_adjoint_traces(data_filenames=syn_filenames,
                                      fmt=self.syn_data_format,
                                      path_out=save_adjsrcs)

        # Return a matching list of observed and synthetic waveform filenames
        observed, synthetic = return_matching_waveform_files(
            obs_path, syn_path, obs_fmt=self.obs_data_format,
            syn_fmt=self.syn_data_format, components=components
        )

        return observed, synthetic

    def _quantify_misfit_single(self, obs_fid, syn_fid, config,
                                save_adjsrcs=False):
        """
        Main Pyatoa processing function to quantify misfit + generation adjsrc.

        Run misfit quantification for a single event-station pair. Gather data,
        preprocess, window and measure data, save adjoint source if
        requested, and then returns the total misfit and the collected
        windows for the station.

        :type obs_fid: str
        :param obs_fid: filename for the observed waveform to be processed
        :type syn_fid: str
        :param syn_fid: filename for the synthetic waveform to be procsesed
        :type config: pyatoa.core.config.Config
        :param config: Config object that defines all the processing parameters
            required by the Pyatoa workflow
        :type save_adjsrcs: str
        :param save_adjsrcs: path to directory where adjoint sources should be
            saved. Filenames will be generated automatically by Pyatoa to fit
            the naming schema required by SPECFEM. If False, no adjoint sources
            will be saved. They of course can be saved manually later using
            Pyatoa + PyASDF
        """
        # Read in data required for processing
        cat = read_events_plus(
            fid=os.path.join(self.path.specfem_data,
                             f"{self._source_prefix}_{config.event_id}"),
            fmt=self._source_prefix
        )
        # parameter `origintime` will only be applied if format=='ASCII'
        obs = read(fid=obs_fid, data_format=self.obs_data_format,
                   origintime=cat[0].preferred_origin().time)
        syn = read(fid=syn_fid, data_format=self.syn_data_format,
                   origintime=cat[0].preferred_origin().time)

        # Attempt to gather station metadata from data directory, or SAC header
        if self.obs_data_format == "SAC":
            inv = get_inv_from_sac_header(obs)
        else:



        # Unique identifier for the given source-receiver pair for file naming
        # Something like: 001_i01_s00_XX_XYZ
        tag = f"{self.ftag(config)}_{syn[0].id.replace('.', '_')}"

        # Configure a single source-receiver pair temporary logger
        station_logger = self._config_auxiliary_logger(
            fid=os.path.join(self.path._tmplogs, f"{tag}.log")
        )
        # Temporary log header is just the name of the station we're after
        station_logger.info(f"\n{'/' * 80}\n{syn[0].id:^80}\n{'/' * 80}")

        # Check whether or not we want to use misfit windows from last eval.
        _fix_win, _msg = self._check_fixed_windows(iteration=config.iteration,
                                                   step_count=config.step_count)
        station_logger.info(_msg)

        # Setup ASDFDataSet in read only so we can pull data/windows in parallel
        ds_fid = os.path.join(self.path["_datasets"], f"{config.event_id}.h5")
        if os.path.exists(ds_fid):
            while True:
                try:
                    ds = ASDFDataSet(ds_fid, mode="r")  # NOTE: read only mode
                    break
                except BlockingIOError:
                    pass
        else:
            ds = None

        mgmt = Manager(st_obs=obs, st_syn=syn, event=cat[0],
                       config=config, ds=ds)

        # If any part of this processing fails, move on to plotting because we
        # will have gathered waveform data so a figure is still useful.
        try:
            mgmt.standardize()
            mgmt.preprocess(remove_response=False, normalize_to="syn")
            mgmt.window(fix_windows=_fix_win)
            mgmt.measure()
        except ManagerError as e:
            station_logger.warning(e)
            pass

        # Plot waveform + map figure. Map may fail if we don't have appropriate
        # metadata, in which case we fall back to plotting waveform only
        if self.plot_waveforms:
            # e.g., 001_i01_s00_XX_ABC.png
            save = os.path.join(self.path["_figures"], f"{tag}.png")
            try:
                mgmt.plot(choice="both", show=False, save=save)
            except ManagerError as e:
                station_logger.warning(e)
                mgmt.plot(choice="wav", show=False, save=save)

        # Write out the .adj adjoint source files for solver to discover.
        if mgmt.stats.misfit is not None and save_adjsrcs:
            mgmt.write_adjsrcs(path=save_adjsrcs, write_blanks=False)

        # Wait until the very end to write to the HDF5 file, then do it
        # pseudo-serially to get around trying to parallel write to HDF5 file
        # WARNING: This is a big potential bottleneck here
        if ds is not None:
            ds._close()  # close the read-only version so we can open in write
        while True:
            try:
                with ASDFDataSet(ds_fid, mode="a") as ds:
                    mgmt.write_to_dataset(ds=ds)
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
        summed_residuals = np.sum(residuals)
        return summed_residuals / self._ntask

    def finalize(self):
        """
        Run serial finalization tasks at the end of a given iteration. These 
        tasks are specific to Pyatoa, used to aggregate figures and data.

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
        if self.plot_waveforms:
            event_pdfs = sorted(glob(os.path.join(self.path._figures,
                                                  "*_*_*.pdf")))

            # Strip off event name to get evaluation tag, i.e.: i01_s00.pdf
            fid_out = "_".join(os.path.basename(event_pdfs[0]).split("_")[1:])
            path_out = os.path.join(self.path._figures, f"{fid_out}")

            # Merge PDFs into a single PDF, delete originals
            merge_pdfs(fids=event_pdfs, fid_out=path_out, remove_fids=True)

        # Move scratch/ directory results into more permanent storage
        if self.export_datasets:
            src = glob(os.path.join(self.path._datasets, "*.h5"))
            src += glob(os.path.join(self.path._datasets, "*.csv"))  # inspector
            dst = os.path.join(self.path.output, "pyaflowa", "datasets", "")
            unix.mkdir(dst)
            unix.cp(src, dst)

        if self.export_figures:
            src = glob(os.path.join(self.path._figures, "*.pdf"))
            dst = os.path.join(self.path.output, "pyaflowa", "figures", "")
            unix.mkdir(dst)
            unix.cp(src, dst)

        if self.export_log_files:
            src = glob(os.path.join(self.path._logs, "*.log"))
            dst = os.path.join(self.path.output, "pyaflowa", "logs", "")
            unix.mkdir(dst)
            unix.cp(src, dst)

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
        :rtype: tuple (bool, str)
        :return: (bool on whether to use windows from the previous step,
            and a message that can be sent to the logger)
        """
        fix_windows = False
        msg = ""

        # First function evaluation never fixes windows
        if iteration == 1 and step_count == 0:
            fix_windows = False
            msg = "first evaluation of workflow, selecting new windows"
        elif isinstance(self.fix_windows, str):
            # By 'iter'ation only pick new windows on the first step count
            if self.fix_windows.upper() == "ITER":
                if step_count == 0:
                    fix_windows = False
                    msg = "first step of line search, will select new windows"
                else:
                    fix_windows = True
                    msg = "mid line search, fix windows from last evaluation"
            # 'Once' picks windows only for the first function evaluation of
            # the current set of iterations.
            elif self.fix_windows.upper() == "ONCE":
                if iteration == self._start and step_count == 0:
                    fix_windows = False
                    msg = "first evaluation of workflow, selecting new windows"
                else:
                    fix_windows = True
                    msg = "mid workflow, fix windows from last evaluation"
        # Bool fix windows simply sets the parameter
        elif isinstance(self.fix_windows, bool):
            fix_windows = self.fix_windows
            msg = f"fixed windows flag set constant: {self.fix_windows}"

        return fix_windows, msg

    def _config_auxiliary_logger(self, fid):
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
        for log in ["pyflex", "pyadjoint", "pysep", "pyatoa"]:
            # Set the overall log level
            logger = logging.getLogger(log)
            # Turn off any existing handlers (stream and file)
            while logger.hasHandlers():
                logger.removeHandler(logger.handlers[0])
            # Log to new temporary file
            logger.setLevel(self.preprocess_log_level)
            logger.addHandler(handler)

        return logger

    def _finalize_logging(self, config, total_windows, total_misfit):
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

        :type config: pyatoa.core.config.Config
        :param config: Config object that will be queried for iteration, step
            count and event ID information
        :type total_windows: int
        :param total_windows: total number of windows collected for a given
            source. this will be written to the final log message
        :type total_misfit: float
        :param total_misfit: total misfit for a given source. this will be
            written to the final log message
        """
        pyatoa_logger = self._config_auxiliary_logger(
            fid=os.path.join(self.path._logs, f"{self.ftag(config)}.log")
        )
        # Summary log message so that User can quickly assess each source
        pyatoa_logger.info(
            f"\n{'=' * 80}\n{'SUMMARY':^80}\n{'=' * 80}\n"
            f"SOURCE NAME: {config.event_id}\n"
            f"WINDOWS: {total_windows}\n"
            f"RAW MISFIT: {total_misfit:.4f}\n"
            f"\n{'=' * 80}\n{'RAW LOGS':^80}\n{'=' * 80}"
            )

        # Collect each station's log file and write them into the main file
        tmp_logs = sorted(glob(os.path.join(self.path._tmplogs,
                                            f"*{config.event_id}_*.log")))

        with open(pyatoa_logger.handlers[0].baseFilename, "a") as fw:
            for tmp_log in tmp_logs:
                with open(tmp_log, "r") as fr:
                    fw.writelines(fr.readlines())
                unix.rm(tmp_log)  # delete after writing



