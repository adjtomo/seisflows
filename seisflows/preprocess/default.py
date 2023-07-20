#!/usr/bin/env python3
"""
The SeisFlows Preprocessing module is in charge of interacting with seismic
data (observed and synthetic). It should contain functionality to read and write
seismic data, apply preprocessing such as filtering, quantify misfit,
and write adjoint sources that are expected by the solver.
"""
import os
import numpy as np
from concurrent.futures import ProcessPoolExecutor, wait
from glob import glob

from obspy import read as obspy_read
from obspy import Stream, Trace, UTCDateTime

from seisflows import logger
from seisflows.tools import signal, unix
from seisflows.tools.config import Dict, get_task_id

from seisflows.plugins.preprocess import misfit as misfit_functions
from seisflows.plugins.preprocess import adjoint as adjoint_sources


class Default:
    """
    Default Preprocess
    ------------------
    Data processing for seismic traces, with options for data misfit,
    filtering, normalization and muting.

    Parameters
    ----------
    :type obs_data_format: str
    :param obs_data_format: data format for reading observed traces into
        memory. Available formats: 'su', 'ascii', 'sac'
    :type unit_output: str
    :param unit_output: Data units. Must match the synthetic output of
        external solver. Available: ['DISP': displacement, 'VEL': velocity,
        'ACC': acceleration, 'PRE': pressure]
    :type misfit: str
    :param misfit: misfit function for waveform comparisons. For available
        see seisflows.plugins.preprocess.misfit
    :type adjoint: str
    :param adjoint: adjoint source misfit function (backprojection function for 
        migration, or the objective function in FWI). For available see
        seisflows.plugins.preprocess.adjoint
    :type normalize: str
    :param normalize: Data normalization parameters used to normalize the
        amplitudes of waveforms. By default, set to NoneType, which means no 
        normalization is applied. User can choose from one of the following 
        options to normalize BOTH `obs` and `syn` data:

        - TNORML1: normalize per trace by the L1 norm of itself
        - TNORML2: normalize per trace by the L2 norm of itself
        - TNORM_MAX: normalize by the maximum positive amplitude in the trace
        - TNORM_ABSMAX: normalize by the absolute maximum amplitude in the trace
        - TNORM_MEAN: normalize by the mean of the absolute trace

        Note: options ENORML? are not currently available. If this is a 
        feature you would like to see, please open a GitHub Issue.
        - ENORML1: normalize per event by L1 of traces; OR
        - ENORML2: normalize per event by L2 of traces;
    :type filter: str
    :param filter: Data filtering type, by default no filtering is applied.
        Available options for user to choose are:

        - BANDPASS (requires: MIN_FREQ + MAX_FREQ OR MIN_PERIOD + MAX PERIOD);
        - LOWPASS (requires: MAX_FREQ OR MIN_PERIOD);
        - HIGHPASS (requires: MIN_FREQ OR MAX_PERIOD);
    :type min_period: float
    :param min_period: Minimum filter period applied to time series.
        See also MIN_FREQ, MAX_FREQ, if User defines FREQ parameters, they
        will overwrite PERIOD parameters.
    :type max_period: float
    :param max_period: Maximum filter period applied to time series. See
        also MIN_FREQ, MAX_FREQ, if User defines FREQ parameters, they will
        overwrite PERIOD parameters.
    :type min_freq: float
    :param min_freq: Maximum filter frequency applied to time series,
        See also MIN_PERIOD, MAX_PERIOD, if User defines FREQ parameters,
        they will overwrite PERIOD parameters.
    :type max_freq: float
    :param max_freq: Maximum filter frequency applied to time series,
        See also MIN_PERIOD, MAX_PERIOD, if User defines FREQ parameters,
        they will overwrite PERIOD parameters.
    :type mute: list
    :param mute: Data mute parameters used to zero out early / late
        arrivals or offsets. Choose any number of:
        
        - EARLY: mute early arrivals;
        - LATE: mute late arrivals;
        - SHORT: mute short source-receiver distances;
        - LONG: mute long source-receiver distances

    Paths
    -----
    :type path_preprocess: str
    :param path_preprocess: scratch path for all preprocessing processes,
        including saving files
    ***
    """
    def __init__(self, syn_data_format="ascii", obs_data_format="ascii",
                 unit_output="VEL", misfit="waveform",
                 adjoint="waveform", normalize=None, filter=None,
                 min_period=None, max_period=None, min_freq=None, max_freq=None,
                 mute=None, early_slope=None, early_const=None, late_slope=None,
                 late_const=None, short_dist=None, long_dist=None,
                 workdir=os.getcwd(), path_preprocess=None, path_solver=None,
                 **kwargs):
        """
        Preprocessing module parameters

        .. note::
            Paths and parameters listed here are shared with other modules and 
            so are not included in the class docstring.

        :type syn_data_format: str
        :param syn_data_format: data format for reading synthetic traces into
            memory. Shared with solver module. Available formats: 'su', 'ascii'
        :type workdir: str
        :param workdir: working directory in which to look for data and store
        results. Defaults to current working directory
        :type path_preprocess: str
        :param path_preprocess: scratch path for all preprocessing processes,
            including saving files
        """
        self.syn_data_format = syn_data_format.upper()
        self.obs_data_format = obs_data_format.upper()
        self.unit_output = unit_output.upper()
        self.misfit = misfit
        self.adjoint = adjoint
        self.normalize = normalize

        self.filter = filter
        self.min_period = min_period
        self.max_period = max_period
        self.min_freq = min_freq
        self.max_freq = max_freq
        self.mute = mute or []
        self.normalize = normalize or []

        # Set the min/max frequencies and periods, frequency takes priority
        if self.filter:
            if self.min_freq is not None:
                self.max_period = 1 / self.min_freq
            elif self.max_period is not None:
                self.min_freq = 1 / self.max_period

            if self.max_freq is not None:
                self.min_period = 1 / self.max_freq
            elif self.min_period is not None:
                self.max_freq = 1 / self.min_period

        # Mute arrivals sub-parameters
        self.early_slope = early_slope
        self.early_const = early_const
        self.late_slope = late_slope
        self.late_const = late_const
        self.short_dist = short_dist
        self.long_dist = long_dist

        self.path = Dict(
            scratch=path_preprocess or os.path.join(workdir, "scratch",
                                                    "preprocess"),
            solver=path_solver or os.path.join(workdir, "scratch", "solver")
        )

        # The list <_obs_acceptable_data_formats> always includes
        # <_syn_acceptable_data_formats> in addition to more formats
        self._syn_acceptable_data_formats = ["SU", "ASCII"]
        self._obs_acceptable_data_formats = ["SU", "ASCII", "SAC"]

        self._acceptable_unit_output = ["DISP", "VEL", "ACC", "PRE"]

        # Misfits and adjoint sources are defined by the available functions
        # in each of these plugin files. Drop hidden variables from dir()
        self._acceptable_misfits = [_ for _ in dir(misfit_functions)
                                    if not _.startswith("_")]
        self._acceptable_adjsrcs = [_ for _ in dir(adjoint_sources)
                                    if not _.startswith("_")]

        # Acceptable preprocessing parameter options
        self._acceptable_norms = {"TNORML1", "TNORML2", "TNORM_MAX",
                                  "TNORM_ABSMAX", "TNORM_MEAN"}  
                                  #, "ENORML1", "ENORML2"}
        self._acceptable_mutes = {"EARLY", "LATE", "LONG", "SHORT"}
        self._acceptable_filters = {"BANDPASS", "LOWPASS", "HIGHPASS"}

        # Internal attributes used to keep track of inversion workflows
        self._iteration = None
        self._step_count = None
        self._source_names = None

    def check(self):
        """ 
        Checks parameters and paths
        """
        if self.misfit:
            assert(self.misfit in self._acceptable_misfits), \
                f"preprocess.misfit must be in {self._acceptable_misfits}"
        if self.adjoint:
            assert(self.adjoint in self._acceptable_adjsrcs), \
                f"preprocess.misfit must be in {self._acceptable_adjsrcs}"

        # Data normalization option
        if self.normalize:
            assert(self.normalize.upper() in self._acceptable_norms)

        # Data muting options
        if self.mute:
            chosen_mutes = [_.upper() for _ in self.mute]
            assert(set(chosen_mutes).issubset(self._acceptable_mutes))
            if "EARLY" in chosen_mutes:
                assert(self.early_slope is not None)
                assert(self.early_slope >= 0.)
                assert(self.early_const is not None)
            if "LATE" in chosen_mutes:
                assert(self.late_slope is not None)
                assert(self.late_slope >= 0.)
                assert(self.late_const is not None)
            if "SHORT" in chosen_mutes:
                assert(self.short_dist is not None)
                assert (self.short_dist > 0)
            if "LONG" in chosen_mutes:
                assert(self.long_dist is not None)
                assert (self.long_dist > 0)

        # Data filtering options that will be passed to ObsPy filters
        if self.filter:
            assert self.filter.upper() in self._acceptable_filters, \
                f"self.filter must be in {self._acceptable_filters}"

            # Check that the correct filter bounds have been set
            if self.filter.upper() == "BANDPASS":
                assert(self.min_freq is not None and
                       self.max_freq is not None), \
                    ("BANDPASS filter PAR.MIN_PERIOD and PAR.MAX_PERIOD or " 
                     "PAR.MIN_FREQ and PAR.MAX_FREQ")
            elif self.filter.upper() == "LOWPASS":
                assert(self.max_freq is not None or
                       self.min_period is not None),\
                    "LOWPASS requires PAR.MAX_FREQ or PAR.MIN_PERIOD"
            elif self.filter.upper() == "HIGHPASS":
                assert(self.min_freq is not None or
                       self.max_period is not None),\
                    "HIGHPASS requires PAR.MIN_FREQ or PAR.MAX_PERIOD"

            # Check that filter bounds make sense, by this point, MIN and MAX
            # FREQ and PERIOD should be set, so we just check the FREQ
            assert(0 < self.min_freq < np.inf), "0 < PAR.MIN_FREQ < inf"
            assert(0 < self.max_freq < np.inf), "0 < PAR.MAX_FREQ < inf"
            assert(self.min_freq < self.max_freq), (
                "PAR.MIN_FREQ < PAR.MAX_FREQ"
            )
    
        # Check that User-chosen data formats are acceptable
        assert(self.syn_data_format.upper() in
                self._syn_acceptable_data_formats), (
            f"synthetic data format must be in "
            f"{self._syn_acceptable_data_formats}"
            )

        assert(self.obs_data_format.upper() in
                self._obs_acceptable_data_formats), (
            f"observed data format must be in "
            f"{self._obs_acceptable_data_formats}"
            )

        assert(self.unit_output.upper() in self._acceptable_unit_output), \
            f"unit output must be in {self._acceptable_unit_output}"

    def setup(self):
        """
        Sets up data preprocessing machinery
        """
        unix.mkdir(self.path.scratch)

    def finalize(self):
        """
        Teardown procedures for the default preprocessing class. Required
        to keep things general because Pyaflowa preprocessing module has
        some finalize procedures.
        """
        pass

    def read(self, fid, data_format):
        """
        Waveform reading functionality. Imports waveforms as Obspy streams

        :type fid: str
        :param fid: path to file to read data from
        :type data_format: str
        :param data_format: format of the file to read data from
        :rtype: obspy.core.stream.Stream
        :return: ObsPy stream containing data stored in `fid`
        """
        st = None
        if data_format.upper() == "SU":
            st = obspy_read(fid, format="SU", byteorder="<")
        elif data_format.upper() == "SAC":
            st = obspy_read(fid, format="SAC")
        elif data_format.upper() == "ASCII":
            st = read_ascii(fid)
        return st

    def write(self, st, fid):
        """
        Waveform writing functionality. Writes waveforms back to format that
        SPECFEM recognizes

        :type st: obspy.core.stream.Stream
        :param st: stream to write
        :type fid: str
        :param fid: path to file to write stream to
        """
        if self.syn_data_format.upper() == "SU":
            for tr in st:
                # Work around for ObsPy data type conversion
                tr.data = tr.data.astype(np.float32)
            max_delta = 0.065535
            dummy_delta = max_delta

            if st[0].stats.delta > max_delta:
                for tr in st:
                    tr.stats.delta = dummy_delta

            # Write data to file
            st.write(fid, format="SU")

        elif self.syn_data_format.upper() == "ASCII":
            for tr in st:
                # Float provides time difference between starttime and default
                time_offset = float(tr.stats.starttime)
                data_out = np.vstack((tr.times() + time_offset, tr.data)).T
                np.savetxt(fid, data_out, ["%13.7f", "%17.7f"])

    def _calculate_misfit(self, **kwargs):
        """Wrapper for plugins.preprocess.misfit misfit/objective function"""
        if self.misfit is not None:
            return getattr(misfit_functions, self.misfit)(**kwargs)
        else:
            return None

    def _generate_adjsrc(self, **kwargs):
        """Wrapper for plugins.preprocess.adjoint source/backproject function"""
        if self.adjoint is not None:
            return getattr(adjoint_sources, self.adjoint)(**kwargs)
        else:
            return None

    def initialize_adjoint_traces(self, data_filenames, output):
        """
        SPECFEM requires that adjoint traces be present for every matching
        synthetic seismogram. If an adjoint source does not exist, it is
        simply set as zeros. This function creates all adjoint traces as
        zeros, to be filled out later. Does this in parallel for speedup

        :type data_filenames: list of str
        :param data_filenames: existing solver waveforms (synthetic) to read.
            These will be copied, zerod out, and saved to path `save`. Should
            come from solver.data_filenames
        :type output: str
        :param output: path to save the new adjoint traces to. Ideally this is
            set to 'solver/traces/adj'
        """
        # Read in a dummy synthetic file and zero out all data to write
        st = self.read(fid=data_filenames[0],
                       data_format=self.syn_data_format).copy()
        for tr in st:
            tr.data *= 0
        # Write adjoint sources in parallel using an empty Stream object
        with ProcessPoolExecutor(max_workers=unix.nproc()) as executor:
            futures = [
                executor.submit(self._write_adjsrc_single, st, fid, output) 
                for fid in data_filenames
                ]
        # Simply wait until this task is completed
        wait(futures)

    def _write_adjsrc_single(self, st, fid, output):
        """Parallelizable function to write out empty adjoint source"""
        adj_fid = self.rename_as_adjoint_source(os.path.basename(fid))
        self.write(st=st, fid=os.path.join(output, adj_fid))

    def quantify_misfit(self, source_name=None, save_residuals=None,
                        export_residuals=None, save_adjsrcs=None,
                        components=None, iteration=1, step_count=0, **kwargs):
        """
        Prepares solver for gradient evaluation by writing residuals and
        adjoint traces. Meant to be called by solver.eval_func().

        Reads in observed and synthetic waveforms, applies optional
        preprocessing, assesses misfit, and writes out adjoint sources and
        STATIONS_ADJOINT file. Processing for each station is done in parallel
        using concurrent.futures.

        .. warning::

            The concurrent processing in this function may fail in the case
            that a User is running N>1 events using the 'Cluster' system but
            on a local workstation, because each event is also run with
            multiprocessing, so their compute may quickly run out of RAM or
            cores. Might need to introduce `max_workers_preproc` and
            `max_workers_solver` to ensure that there is a good balance
            between the two values.

        :type source_name: str
        :param source_name: name of the event to quantify misfit for. If not
            given, will attempt to gather event id from the given task id which
            is assigned by system.run()
        :type save_residuals: str
        :param save_residuals: if not None, path to write misfit/residuls to
        :type export_residuals: bool
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
        """
        # Retrieve matching obs and syn trace filenames to run through misfit
        # and initialize empty adjoint sources
        obs, syn = self._setup_quantify_misfit(source_name, save_adjsrcs,
                                               components)
        
        # Process each pair in parallel. Max workers is the total num. of cores
        # !!! see note in docstring !!!
        with ProcessPoolExecutor(max_workers=unix.nproc()) as executor:
            futures = [
                executor.submit(self._quantify_misfit_single, o, s, 
                                save_residuals, save_adjsrcs)
                for (o, s) in zip(obs, syn)
            ]
        if save_residuals:
            # Results are returned in the order that they were submitted and
            # writte to text file for other modules to find
            with open(save_residuals, "a") as f:
                for future in futures:
                    residual = future.result()
                    f.write(f"{residual:.2E}\n")
        else:
            # Simply wait for all processes to finish before proceeding
            wait(futures)

        # Exporting residuals to disk (output/) for more permanent storage
        if export_residuals:
            if not os.path.exists(export_residuals):
                unix.mkdir(export_residuals)
            unix.cp(src=save_residuals, dst=export_residuals)

    def _setup_quantify_misfit(self, source_name, save_adjsrcs=None,
                               components=None):
        """
        Gather a list of filenames of matching waveform IDs that can be
        run through the misfit quantification step. Perform some checks to
        ensure that the data is provided in usable formats and that
        preprocessing is only run for synthetics which have a matching
        observation.

        .. note::

            Obs and syn waveform files are expected to be in the format
            NN.SSS.CCc.* (N=network, S=station, C=channel, c=component;
            following SPECFEM ASCII formatting). They will be matched on
            `NN.SSS.c` (dropping channel naming because SEED convention may have
            different channel naming). For example, synthetic name
            'AA.S001.BXZ.semd' will be converted to 'AA.S001.Z', and matching
            observation 'AA.S001.HHZ.SAC' will be converted to 'AA.S001.Z'.
            These two will be matched.

        .. note::

            Assumes the directory structure that has been set up by the solver.
            That is, that waveforms are stored in directories:
            `scratch/solver/{source_name}/traces/obs`  for observations
            `scratch/solver/{source_name}/traces/syn`  for synthetics

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
        # Get organized by looking for available data
        source_name = source_name or self._source_names[get_task_id()]

        obs_path = os.path.join(self.path.solver, source_name, "traces", "obs")
        syn_path = os.path.join(self.path.solver, source_name, "traces", "syn")

        observed = sorted(os.listdir(obs_path))
        synthetic = sorted(os.listdir(syn_path))

        logger.debug(f"found {len(observed)} obs and {len(synthetic)} syn "
                     f"waveforms for event {source_name}")

        assert (len(observed) != 0 and len(synthetic) != 0), \
            f"cannot quantify misfit, missing observed or synthetic traces"

        # Initialize empty adjoint sources for all synthetics that may or may
        # not be overwritten by the misfit quantification step
        if save_adjsrcs is not None:
            syn_filenames = glob(os.path.join(syn_path, "*"))
            self.initialize_adjoint_traces(data_filenames=syn_filenames,
                                           output=save_adjsrcs)

        # Verify observed traces format is acceptable within this module
        obs_ext = list(set([os.path.splitext(x)[-1] for x in observed]))
        assert(len(obs_ext) == 1), (
            f"'{source_name}/traces/obs' has > 1 file formats available, but "
            f"only 1 ({self.obs_data_format}) was expected"
        )
        # Check if the expected file format matches the provided one
        if self.obs_data_format.upper() == "ASCII":
            obs_ext_ok = obs_ext[0].upper() in [".ASCII",
                                                f".SEM{self.unit_output[0]}"]
        else:
            obs_ext_ok = obs_ext[0].upper() == f".{self.obs_data_format}"
        assert obs_ext_ok, (
            f"{source_name}/traces/obs unexpected file format "
            f"{obs_ext[0].upper()} != {self.obs_data_format}"
        )

        # Do the same checks but for the synthetic waveforms
        syn_ext = list(set([os.path.splitext(x)[-1] for x in synthetic]))
        assert(len(syn_ext) == 1), (
            f"'{source_name}/traces/syn' has > 1 file formats available, but "
            f"only 1 ({self.syn_data_format}) was expected"
        )
        # Check if the expected file format matches the provided one
        if self.syn_data_format.upper() == "ASCII":
            syn_ext_ok = syn_ext[0].upper() in [".ASCII",
                                                f".SEM{self.unit_output[0]}"]
        else:
            syn_ext_ok = syn_ext[0].upper() == f".{self.syn_data_format}"
        assert syn_ext_ok, (
            f"{source_name}/traces/syn unexpected file format "
            f"{syn_ext[0].upper()} != {self.syn_data_format}"
        )

        # fmt path/to/NN.SSS.CCc* -> NN.SSS.c (see docstring note for details)
        match_obs = self._curtail_fids_for_file_matching(observed)
        match_syn = self._curtail_fids_for_file_matching(synthetic)

        # only return traces that have both observed and synthetic file match
        match_traces = \
                sorted(list(set(match_syn).intersection(set(match_obs))))
        logger.info(f"{len(match_traces)} traces matching between obs and syn")

        # Curtail based on user-chosen components if required
        if components:
            match_comps = []
            for comp in components:
                match_comps += [_ for _ in match_traces if _.endswith(comp)]
            match_traces = match_comps
            logger.info(f"{len(match_traces)} traces with comps {components}")

        # Final check to makes sure that we still have data that can be compared
        assert(len(match_traces) != 0), (
            f"there are no traces with both observed and synthetic files for "
            f"source: {source_name}; verify that 'traces/obs' and 'traces/syn' "
            f"have the format 'NN.SSS.CCc*', and match on variables 'N', 'S', "
            f"and 'c'"
        )
        logger.info(f"{source_name} has {len(match_traces)} matching traces "
                    f"for preprocessing")

        # Generate the list of full path waveform fids for matching obs + syn
        obs_paths, syn_paths = [], []
        for short_fid in match_traces:
            # Find the corresponding full path name based on the short match vrs
            obs_fid = observed[match_obs.index(short_fid)]
            obs_paths.append(os.path.join(obs_path, obs_fid))

            syn_fid = synthetic[match_syn.index(short_fid)]
            syn_paths.append(os.path.join(syn_path, syn_fid))

        return obs_paths, syn_paths

    def _quantify_misfit_single(self, obs_fid, syn_fid, save_residuals=None,
                                save_adjsrcs=None):
        """
        Run misfit quantification for one pair of data-synthetic waveforms.
        This is kept in a separate function so that it can be parallelized for
        more efficient processing.

        Order of processing steps is:
        resample -> filter (optional) -> mute (optional) -> normalize (optional)
            -> calculate misfit (optional) -> create adjoint source (optional)

        :type obs_fid: str
        :param obs_fid: filename for the observed waveform to be processed
        :type syn_fid: str
        :param syn_fid: filename for the synthetic waveform to be procsesed
        :type save_residuals: str
        :param save_residuals: if not None, path to write misfit/residuls to
        :type save_adjsrcs: str
        :param save_adjsrcs: if not None, path to write adjoint sources to
        :rtype: float
        :return: residual value, calculated by the chosen misfit function
            comparing `obs` and `syn`
        """
        # Read in waveforms based on the User-defined format(s)
        obs = self.read(fid=obs_fid, data_format=self.obs_data_format)
        syn = self.read(fid=syn_fid, data_format=self.syn_data_format)

        # Process observations and synthetics identically
        obs, syn = self._apply_resample(obs, syn)
        if self.filter:
            obs = self._apply_filter(obs)
            syn = self._apply_filter(syn)
        if self.mute:
            obs = self._apply_mute(obs)
            syn = self._apply_mute(syn)
        if self.normalize:
            obs = self._apply_normalize(obs)
            syn = self._apply_normalize(syn)

        # Write the residuals/misfit and adjoint sources for each component
        # The assumption here is that `obs` and `syn` are length=1
        residual = 0
        for tr_obs, tr_syn in zip(obs, syn):
            # Simple check to make sure zip retains ordering. Only works if
            # both syn and data have component stat. This may not be the case
            # for poorly labelled data
            if tr_obs.stats.component and tr_syn.stats.component:
                assert (tr_obs.stats.component == tr_syn.stats.component), (
                    f"{obs_fid}, {syn_fid}"
                    f"Mismatched components for '{os.path.basename(obs_fid)}' "
                    f"obs: `{tr_obs.stats.component}` != " 
                    f"syn: `{tr_syn.stats.component}`. Please check `obs` data"
                    )

            # Calculate the misfit value and write to file
            if save_residuals and self._calculate_misfit:
                residual = self._calculate_misfit(
                    obs=tr_obs.data, syn=tr_syn.data,
                    nt=tr_syn.stats.npts, dt=tr_syn.stats.delta
                )
                logger.debug(f"{tr_syn.get_id()} residual={residual:.3E}")
            else:
                residual = 0

            # Generate an adjoint source trace, write to file in scratch dir.
            if save_adjsrcs and self._generate_adjsrc:
                adjsrc = tr_syn.copy()
                adjsrc.data = self._generate_adjsrc(
                    obs=tr_obs.data, syn=tr_syn.data,
                    nt=tr_syn.stats.npts, dt=tr_syn.stats.delta
                )
                adjsrc = Stream(adjsrc)
                fid = os.path.basename(syn_fid)
                fid = self.rename_as_adjoint_source(fid)
                self.write(st=adjsrc, fid=os.path.join(save_adjsrcs, fid))
                logger.debug(f"writing adjoint source: {fid}")

        return residual

    def _curtail_fids_for_file_matching(self, fid_list, components=None):
        """
        Convenience function to convert NN.SSS.CCc.* -> NN.SSS.c and also
        curtail list of file IDs by component

        :type fid_list: list of str
        :param fid_list: list of file path/IDs that need to be shortened
        :type components: list
        :param components: optional list of components to ignore preprocessing
            traces that do not have matching components. The adjoint sources for
            these components will be 0. E.g., ['Z', 'N']. If None, all available
            components will be considered.
        :rtype: list of str
        :return: list of shortened file IDs
        """
        # Drop full path incase these are given as absolute paths
        full_fids = [os.path.basename(_) for _ in fid_list]
        # Split into expected format NN.SSS.CCc, drop extension
        fids = [_.split(".")[:3] for _ in full_fids]
        short_fids = []
        for fid in fids:
            net, sta, cha = fid
            comp = cha[-1]
            # NN.SSS.c
            short_fids.append(f"{net}.{sta}.{comp}")

        return short_fids

    def rename_as_adjoint_source(self, fid):
        """
        Rename synthetic waveforms into filenames consistent with how SPECFEM
        expects adjoint sources to be named. Usually this just means adding
        a '.adj' to the end of the filename

        :type fid: str
        :param fid: file path to rename for adjoint source
        :rtype: str
        :return: renamed file that matches expected SPECFEM filename format
            for adjoint sources
        """
        if not fid.endswith(".adj"):
            if self.syn_data_format.upper() == "SU":
                fid = f"{fid}.adj"
            elif self.syn_data_format.upper() == "ASCII":
                # Differentiate between SPECFEM3D and 3D_GLOBE file naming
                # SPECFEM3D: NN.SSSS.CCC.sem?
                # SPECFEM3D_GLOBE: NN.SSSS.CCC.sem.ascii
                ext = os.path.splitext(fid)[-1]
                # SPECFEM3D
                if ".sem" in ext:
                    fid = fid.replace(ext, ".adj")
                # GLOBE (!!! Hardcoded to only work with ASCII format)
                elif ext == ".ascii":
                    root, ext1 = os.path.splitext(fid)  # .ascii
                    root, ext2 = os.path.splitext(root)  # .sem
                    fid = fid.replace(f"{ext2}{ext1}", ".adj")

        return fid

    def finalize(self):
        """Teardown procedures for the default preprocessing class"""
        pass

    @staticmethod
    def sum_residuals(residuals):
        """
        Returns the summed square of residuals for each event. Following
        Tape et al. 2007

        :type residuals: np.array
        :param residuals: list of residuals from each NTASK event
        :rtype: float
        :return: sum of squares of residuals
        """
        return np.sum(residuals ** 2.)

    def _apply_resample(self, st_a, st_b):
        """
        Resample all traces in `st_a` to the sampling rate of `st_b`. Resamples
        one to one, that is each trace in `obs` is resampled to the
        corresponding indexed trace in `syn`

        :type st_a: obspy.core.stream.Stream
        :param st_a: stream to be resampled using sampling rates from `st_b`
        :type st_b: obspy.core.stream.Stream
        :param st_b:  stream whose sampling rates will be used to resample
            `st_a`. Usually this is the synthetic data
        :rtype: (obspy.core.stream.Stream, obspy.core.stream.Stream)
        :return: `st_a` (resampled), `st_b`
        """
        for tr_a, tr_b in zip(st_a, st_b):
            sr_a = tr_a.stats.sampling_rate
            sr_b = tr_b.stats.sampling_rate
            # Is this the correct resampling method to use?
            if sr_a != sr_b:
                logger.debug(f"resampling '{tr_a.get_id()}' {sr_a}->{sr_b} Hz")
                # Resample in place
                tr_a.resample(sampling_rate=tr_b.stats.sampling_rate)

        return st_a, st_b

    def _apply_filter(self, st):
        """
        Apply a filter to waveform data using ObsPy, throw on a standard 
        demean, detrened and taper prior to filtering. Options for different
        filtering types. Uses default filter options from ObsPy.
        
        Zerophase enforced to be True to avoid phase shifting data.

        :type st: obspy.core.stream.Stream
        :param st: stream to be filtered
        :rtype: obspy.core.stream.Stream
        :return: filtered traces
        """
        # Pre-processing before filtering
        st.detrend("demean")
        st.detrend("linear")
        st.taper(0.05, type="hann")

        if self.filter.upper() == "BANDPASS":
            st.filter("bandpass", zerophase=True, freqmin=self.min_freq,
                      freqmax=self.max_freq)
        elif self.filter.upper() == "LOWPASS":
            st.filter("lowpass", zerophase=True, freq=self.max_freq)
        elif self.filter.upper() == "HIGHPASS":
            st.filter("highpass", zerophase=True, freq=self.min_freq)

        return st

    def _apply_mute(self, st):
        """
        Apply mute on data based on early or late arrivals, and short or long
        source receiver distances

        .. note::

            The underlying mute functions have been refactored but not tested
            as I was not aware of the intended functionality. Not gauranteed
            to work, use at your own risk.

        :type st: obspy.core.stream.Stream
        :param st: stream to mute
        :rtype: obspy.core.stream.Stream
        :return: muted stream object
        """
        mute_choices = [_.upper() for _ in self.mute]
        if "EARLY" in mute_choices:
            st = signal.mute_arrivals(st, slope=self.early_slope,
                                      const=self.early_const, choice="EARLY")
        if "LATE" in mute_choices:
            st = signal.mute_arrivals(st, slope=self.late_slope,
                                      const=self.late_const, choice="LATE")
        if "SHORT" in mute_choices:
            st = signal.mute_offsets(st, dist=self.short_dist, choice="SHORT")
        if "LONG" in mute_choices:
            st = signal.mute_offsets(st, dist=self.long_dist, choice="LONG")

        return st

    def _apply_normalize(self, st):
        """
        Normalize the amplitudes of waveforms based on user choice

        .. warning::

            Event normalization does not currently work as this requires 
            access to all waveform simultaneously whereas we do processing
            trace by trace. Will need to devise a method for calculating this
            in the future. For now, option has been remoevd

        :type st: obspy.core.stream.Stream
        :param st: All of the data streams to be normalized
        :rtype: obspy.core.stream.Stream
        :return: stream with normalized traces
        """
        # Normalize each trace by its L1 norm
        if self.normalize.upper() == "TNORML1":
            for tr in st:
                w = np.linalg.norm(tr.data, ord=1)
                if w < 0:
                    logger.warning(f"CAUTION: L1 Norm for {tr.get_id()} is "
                                   f"negative, this will result in "
                                   f"unintentional sign flip")
                tr.data /= w
        # Normalize each trace by its L2 norm
        elif self.normalize.upper() == "TNORML2":
            for tr in st:
                w = np.linalg.norm(tr.data, ord=2)
                if w < 0:
                    logger.warning(f"CAUTION: L2 Norm for {tr.get_id()} is "
                                   f"negative, this will result in "
                                   f"unintentional sign flip")
                tr.data /= w
        # Normalize each trace by its maximum positive amplitude
        elif self.normalize.upper() == "TNORM_MAX":
            for tr in st:
                w = np.max(tr.data)
                tr.data /= w
        # Normalize each trace by the maximum amplitude (neg or pos) 
        elif self.normalize.upper() == "TNORM_ABSMAX":
            for tr in st:
                w = np.abs(tr.max())
                tr.data /= w
        # Normalize by the mean of absolute trace amplitudes
        elif self.normalize.upper() == "TNORM_MEAN":
            for tr in st:
                w = np.mean(np.abs(tr.data))
                tr.data /= w

        # !!! These are not currently accessible. Open a GitHub issue if you
        # !!! would like to see event-based normalization
        # Normalize an event by the L1 norm of all traces
        # if 'ENORML1' in norm_choices:
        #     w = 0.
        #     for tr in st_out:
        #         w += np.linalg.norm(tr.data, ord=1)
        #     for tr in st_out:
        #         tr.data /= w
        # # Normalize an event by the L2 norm of all traces
        # elif "ENORML2" in norm_choices:
        #     w = 0.
        #     for tr in st_out:
        #         w += np.linalg.norm(tr.data, ord=2)
        #     for tr in st_out:
        #         tr.data /= w

        return st


def read_ascii(fid, origintime=None):
    """
    Read waveforms in two-column ASCII format. This is copied directly from
    pyatoa.utils.read.read_sem()
    """
    try:
        times = np.loadtxt(fname=fid, usecols=0)
        data = np.loadtxt(fname=fid, usecols=1)

    # At some point in 2018, the Specfem developers changed how the ascii files
    # were formatted from two columns to comma separated values, and repeat
    # values represented as 2*value_float where value_float represents the data
    # value as a float
    except ValueError:
        times, data = [], []
        with open(fid, 'r') as f:
            lines = f.readlines()
        for line in lines:
            try:
                time_, data_ = line.strip().split(',')
            except ValueError:
                if "*" in line:
                    time_ = data_ = line.split('*')[-1]
                else:
                    raise ValueError
            times.append(float(time_))
            data.append(float(data_))

        times = np.array(times)
        data = np.array(data)

    if origintime is None:
        origintime = UTCDateTime("1970-01-01T00:00:00")

    # We assume that dt is constant after 'precision' decimal points
    delta = round(times[1] - times[0], 4)

    # Honor that Specfem doesn't start exactly on 0
    origintime += times[0]

    # Write out the header information. Deal with the fact that SPECFEM2D/3D and
    # 3D_GLOBE have slightly different formats for their filenames
    net, sta, cha, *fmt = os.path.basename(fid).split('.')
    stats = {"network": net, "station": sta, "location": "",
             "channel": cha, "starttime": origintime, "npts": len(data),
             "delta": delta, "mseed": {"dataquality": 'D'},
             "time_offset": times[0], "format": fmt[0]
             }
    st = Stream([Trace(data=data, header=stats)])

    return st
