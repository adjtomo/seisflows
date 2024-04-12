#!/usr/bin/env python3
"""
The SeisFlows Preprocessing module is in charge of interacting with seismic
data (observed and synthetic). It should contain functionality to read and write
seismic data, apply preprocessing such as filtering, quantify misfit,
and write adjoint sources that are expected by the solver.
"""
import os
import sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor, wait
from glob import glob

from obspy import read as obspy_read
from obspy import Stream, Trace, UTCDateTime

from seisflows import logger
from seisflows.tools import unix
from seisflows.tools.config import Dict, get_task_id
from seisflows.tools.graphics import plot_waveforms
from seisflows.tools.signal import normalize, resample, filter, mute, trim
from seisflows.tools.specfem import (rename_as_adjoint_source,
                                     return_matching_waveform_files)

from seisflows.plugins.preprocess import misfit as misfit_functions
from seisflows.plugins.preprocess import adjoint as adjoint_sources


class Default:
    """
    Default Preprocess [Preprocess Base]
    ------------------------------------
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
        - RNORM_OBS_MAX: normalize synthetic traces by the maximum amplitude of
            the corresponding observed trace
        - RNORM_OBS_ABSMAX: normalize synthetic traces by the absolute maximum
            amplitude of the corresponding observed trace
        - RNORM_SYN_MAX: normalize observed traces by the maximum amplitude of
            the corresponding synthetic trace
        - RNORM_SYN_ABSMAX: normalize observed traces by the absolute maximum
            amplitude of the corresponding synthetic trace

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
        arrivals or offsets. Choose the following:
        - EARLY: mute early arrivals
        - LATE: mute late arrivals;
        - SHORT: mute short source-receiver distances;
        - LONG: mute long source-receiver distances
        input comma separated list of strings, (e.g., mute: early,late,short)
    :type plot_waveforms: bool
    :param plot_waveforms: plot waveforms from each evaluation of the misfit.
        By default turned off as this can produce many files for an interative
        inversion.

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
                 plot_waveforms=False, source_prefix=None,
                 workdir=os.getcwd(), path_preprocess=None, path_solver=None,
                 path_specfem_data=None,
                 **kwargs):
        """
        Preprocessing module parameters

        .. note::

            Paths and parameters listed here are shared with other modules and
            so are not included in the class docstring.

        :type syn_data_format: str
        :param syn_data_format: data format for reading synthetic traces into
            memory. Shared with solver module. Available formats: 'su', 'ascii'
        :type source_prefix: str
        :param source_prefix: prefix of source/event/earthquake files. Used for
            reading and determining source locations for preprocessing
        :type workdir: str
        :param workdir: working directory in which to look for data and store
        results. Defaults to current working directory
        :type path_solver: str
        :param path_solver: scratch path for solver used to access trace files
        :type path_specfem_data: str
        :param path_specfem_data: path to SPECFEM DATA/ directory which must
            contain the CMTSOLUTION, STATIONS and Par_file files used for
            running SPECFEM
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

        # Miscellaneous paramters
        self.source_prefix = source_prefix
        self.plot_waveforms = plot_waveforms

        self.path = Dict(
            scratch=path_preprocess or os.path.join(workdir, "scratch",
                                                    "preprocess"),
            solver=path_solver or os.path.join(workdir, "scratch", "solver"),
            specfem_data=path_specfem_data or None
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
                                  "TNORM_ABSMAX", "TNORM_MEAN", "RNORM_SYN_MAX",
                                  "RNORM_SYN_ABSMAX", "RNORM_OBS_MAX",
                                  "RNORM_OBS_ABSMAX"}
                                  #, "ENORML1", "ENORML2"}
        self._acceptable_mutes = {"EARLY", "LATE", "LONG", "SHORT"}
        self._acceptable_filters = {"BANDPASS", "LOWPASS", "HIGHPASS"}

        # Internal attributes used to keep track of inversion workflows
        self._iteration = None
        self._step_count = None

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

        # This is a redundant check on the DATA/STATIONS file (solver also
        # runs this check). This is required by noise workflows to determine
        # station rotation
        assert (self.path.specfem_data is not None and
                os.path.exists(self.path.specfem_data)), (
            f"`path_specfem_data` must exist and must point to directory "
            f"containing SPECFEM input files"
        )
        # Ensure STATIONS files exist as the locations are used for preproc,
        assert(os.path.exists(
            os.path.join(self.path.specfem_data, "STATIONS"))), (
            f"DATA/STATIONS does not exist but is required by preprocessing"
        )
        # Ensure source files exist as their locations are used for preproc.
        assert(glob(os.path.join(self.path.specfem_data,
                                 f"{self.source_prefix}_*"))), (
            f"DATA/{self.source_prefix}_* does not exist but is required"
        )

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

    def quantify_misfit(self, source_name=None, save_residuals=None,
                        export_residuals=None, save_adjsrcs=None,
                        components=None, iteration=1, step_count=0,
                        _serial=False, **kwargs):
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
        :param iteration: optional, current iteration of the workflow used for
            tagging output waveform figures if internal parameter 
            `plot_waveforms` is set
        :type step_count: int
        :param step_count: optional, current step count of the line search, 
            used for tagging output waveform figures if internal parameter
            `plot_waveforms` is set
        :type _serial: bool
        :param _serial: debug function to turn preprocessing to a serial task
            whereas it is normally a multiprocessed parallel task
        """
        self._iteration = iteration
        self._step_count = step_count

        # Retrieve matching obs and syn trace filenames to run through misfit
        # and initialize empty adjoint sources
        obs, syn = self._setup_quantify_misfit(source_name, save_adjsrcs,
                                               components)

        # Process each pair in parallel. Max workers is the total num. of cores
        # !!! see note in docstring !!!
        residuals = []
        if _serial:
            for o, s in zip(obs, syn):
                residual = self._quantify_misfit_single(o, s, source_name, 
                                                        save_residuals, 
                                                        save_adjsrcs)
                residuals.append(residual)
        # Process each pair in parallel. Max workers is total num. of cores        
        else:
            with ProcessPoolExecutor(max_workers=unix.nproc()) as executor:
                futures = [
                    executor.submit(self._quantify_misfit_single, o, s,
                                    source_name, save_residuals, save_adjsrcs)
                    for (o, s) in zip(obs, syn)
                ]
            wait(futures)
            # Results are returned in the order that they were submitted and
            for future in futures:
                try:
                    residual = future.result()
                    residuals.append(residual)
                except Exception as e:
                    logger.critical(f"PREPROC FAILED: {e}")
                    sys.exit(-1)

        # Write residuals to text file for other modules to find
        if save_residuals:
            with open(save_residuals, "a") as f:
                for residual in residuals:
                    f.write(f"{residual:.2E}\n")

            # Exporting residuals to disk (output/) for more permanent storage
            if export_residuals:
                unix.mkdir(export_residuals)
                unix.cp(src=save_residuals, dst=export_residuals)

        logger.info(f"FINISH QUANTIFY MISFIT: {source_name}")

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

    def _quantify_misfit_single(self, obs_fid, syn_fid, source_name=None,
                                save_residuals=None, save_adjsrcs=None):
        """
        Run misfit quantification for one pair of data-synthetic waveforms.
        This is kept in a separate function so that it can be parallelized for
        more efficient processing.

        :type obs_fid: str
        :param obs_fid: filename for the observed waveform to be processed
        :type syn_fid: str
        :param syn_fid: filename for the synthetic waveform to be procsesed
        :type source_name: str
        :param source_name: name of the source used for tagging output waveform
            figures if internal paramter `plot` is set True
        :type save_residuals: str
        :param save_residuals: if not None, path to write misfit/residuls to
        :type save_adjsrcs: str
        :param save_adjsrcs: if not None, path to write adjoint sources to
        :rtype: float
        :return: residual value, calculated by the chosen misfit function
            comparing `obs` and `syn`
        """
        # Read in waveforms based on the User-defined format(s)
        obs = read(fid=obs_fid, data_format=self.obs_data_format)
        syn = read(fid=syn_fid, data_format=self.syn_data_format)

        # Wrap all the preprocessing functions into single function so that
        # it can be more easily overwritten by overwriting classes
        # Default order of processing steps is:
        # resample, filter (optional), (optional), normalize (optional)
        obs, syn = self.preprocess(obs, syn)

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
                # Calculate misfit from measurement (e.g., Tromp et al. 2005 
                # Eq. 1), ensuring misfit is a positive value so that it can be
                # correctly compared during line search
                residual = 1/2 * residual ** 2
                logger.debug(f"{tr_syn.get_id()} misfit={residual:.3E}")
            else:
                residual = 0

            # Generate an adjoint source trace, write to file in scratch dir.
            # SPECFEM expects non time-reversed adjoint sources
            if save_adjsrcs and self._generate_adjsrc:
                adjsrc = tr_syn.copy()
                adjsrc.data = self._generate_adjsrc(
                    obs=tr_obs.data, syn=tr_syn.data,
                    nt=tr_syn.stats.npts, dt=tr_syn.stats.delta
                )
                fid = os.path.basename(syn_fid)
                fid = rename_as_adjoint_source(fid, fmt=self.syn_data_format)
                write(st=Stream(adjsrc), fid=os.path.join(save_adjsrcs, fid),
                      data_format=self.syn_data_format)
            else:
                adjsrc = None

            # Plot waveforms showing observed, synthetic and adjoint source
            if self.plot_waveforms:
                # If source name is not provided by calling function, assign to
                # the given task ID which can later be used to track source name
                if source_name is None:
                    source_name = get_task_id()

                # Determine how to tag the output files
                fig_path = os.path.join(self.path.scratch, source_name)
                if not os.path.exists(fig_path):
                    unix.mkdir(fig_path)
                # e.g., path/to/figure/NN_SSS_LL_CCC_IS.png 
                if self._iteration is not None:
                    _itr = f"_i{self._iteration:0>2}"
                else:
                    _itr = ""
                if self._step_count is not None:
                    _stp = f"s{self._step_count:0>2}"
                else:
                    _stp = ""
                fid_out = os.path.join(
                    fig_path, f"{tr_syn.id.replace('.', '_')}{_itr}{_stp}.png"
                )
                title = f"{tr_syn.id}; misfit={residual:.3E}"
                plot_waveforms(tr_obs=tr_obs, tr_syn=tr_syn, tr_adj=adjsrc,
                               fid_out=fid_out, title=title)

        return residual

    def preprocess(self, obs, syn):
        """
        Convenience function that wraps all preprocessing steps so that
        they can be more easily overwritten if the User wants preprocessing
        that does not follow this convention

        :type obs: obspy.core.stream.Stream
        :param obs: Stream containing observed waveforms
        :type syn: obspy.core.stream.Stream
        :param syn: Stream containing synthetic waveforms
        """
        # Resample observed waveforms to the same sampling rate as the
        # synthetics because the output adjoint sources will need this samp rate
        obs, syn = resample(st_a=obs, st_b=syn)

        # Trim the observed seismograms to the length of the synthetics. 
        # Returned in the same order as input so syn first since we are trimming
        # to syn
        syn, obs = trim(st=syn, st_trim=obs)

        # Apply some basic detrends to clean up data
        for st in [obs, syn]:
            st.detrend("demean")
            st.detrend("linear")
            st.taper(0.05, type="hann")

        if self.filter:
            obs = filter(obs, choice=self.filter, min_freq=self.min_freq,
                         max_freq=self.max_freq)
            syn = filter(syn, choice=self.filter, min_freq=self.min_freq,
                         max_freq=self.max_freq)
        if self.mute:
            obs = mute(obs)
            syn = mute(syn)
        if self.normalize:
            if "RNORM_SYN" in self.normalize:
                # normalize `obs` relative to the synthetic trace only
                obs = normalize(st=obs, choice=self.normalize, st_rel=syn)
            elif "RNORM_OBS" in self.normalize:
                # normalize `syn` relative to the observed trace only
                syn = normalize(st=syn, choice=self.normalize, st_rel=obs)
            else:
                # else, normalize traces relative to themselves
                obs = normalize(st=obs, choice=self.normalize)
                syn = normalize(st=syn, choice=self.normalize)

        return obs, syn

def read(fid, data_format, **kwargs):
    """
    Waveform reading functionality. Imports waveforms as Obspy streams

    :type fid: str
    :param fid: path to file to read data from
    :type data_format: str
    :param data_format: format of the file to read data from
    :rtype: obspy.core.stream.Stream
    :return: ObsPy stream containing data stored in `fid`
    :raises TypeError: if the provided data format cannot be read
    """
    if data_format.upper() == "SU":
        st = obspy_read(fid, format="SU", byteorder="<", **kwargs)
    elif data_format.upper() == "ASCII":
        st = read_ascii(fid, **kwargs)
    else:
        st = obspy_read(fid, format=data_format.upper(), **kwargs)
    return st

def write(st, fid, data_format):
    """
    Waveform writing functionality. Writes waveforms back to format that
    SPECFEM recognizes

    :type st: obspy.core.stream.Stream
    :param st: stream to write
    :type fid: str
    :param fid: path to file to write stream to
    :type data_format: str
    :param data_format: format of the file to write to
    """
    if data_format.upper() == "SU":
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

    elif data_format.upper() == "ASCII":
        for tr in st:
            # Time offset should have been set by `read_ascii` when reading in
            # the original ASCII waveform file
            try:
                time_offset = tr.stats.time_offset
            except AttributeError:
                time_offset = 0
            data_out = np.vstack((tr.times() + time_offset, tr.data)).T
            np.savetxt(fid, data_out, ["%13.7f", "%17.7f"])

def read_ascii(fid, origintime=None, **kwargs):
    """
    Converts SPECFEM synthetics into ObsPy Stream objects with the correct
    header information.

    .. note::

        This is a trimmed down version of pysep.utils.io.read_sem() which is
        copied here for better visibility within SeisFlows, and ignores things
        like SAC headers appending

    :type fid: str
    :param fid: path of the given ascii file
    :type origintime: obspy.UTCDateTime
    :param origintime: UTCDatetime object for the origintime of the event
    :rtype st: obspy.Stream.stream
    :return st: stream containing header and data info taken from ascii file
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
        origintime = UTCDateTime("2000-01-01T00:00:00")

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

def initialize_adjoint_traces(data_filenames, fmt, path_out="./"):
    """
    SPECFEM requires that adjoint traces be present for every matching
    synthetic seismogram. If an adjoint source does not exist, it is
    simply set as zeros. This function creates all adjoint traces as
    zeros, to be filled out later. Does this in parallel for speedup

    :type data_filenames: list of str
    :param data_filenames: existing solver waveforms (synthetic) to read.
        These will be copied, zerod out, and saved to path `save`. Should
        come from solver.data_filenames
    :type fmt: str
    :param fmt: format of the input waveforms which will be fed to the `read`
        function to get a working Stream object used to write empty adjsrcs
    :type path_out: str
    :param path_out: path to save the new adjoint traces to. Ideally this is
        set to 'solver/traces/adj'
    """
    # Read in a dummy synthetic file and zero out all data to write
    st = read(fid=data_filenames[0], data_format=fmt).copy()
    for tr in st:
        tr.data *= 0
    # Write adjoint sources in parallel using an empty Stream object
    with ProcessPoolExecutor(max_workers=unix.nproc()) as executor:
        futures = [
            executor.submit(_write_adjsrc_single, st, fid, path_out, fmt)
            for fid in data_filenames
            ]
    # Simply wait until this task is completed
    wait(futures)

def _write_adjsrc_single(st, fid, output, fmt):
    """Parallelizable function to write out empty adjoint source"""
    adj_fid = rename_as_adjoint_source(os.path.basename(fid), fmt=fmt)
    write(st=st, fid=os.path.join(output, adj_fid), data_format=fmt)

