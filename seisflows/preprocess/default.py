#!/usr/bin/env python3
"""
The SeisFlows Preprocessing module is in charge of interacting with seismic
data (observed and synthetic). It should contain functionality to read and write
seismic data, apply preprocessing such as filtering, quantify misfit,
and write adjoint sources that are expected by the solver.
"""
import os
import obspy
import numpy as np

from seisflows import logger
from seisflows.tools import signal, unix
from seisflows.plugins.preprocess import adjoint, misfit, readers, writers


class Default:
    """
    Default SeisFlows preprocessing class

    Provides data processing functions for seismic traces, with options for
    data misfit, filtering, normalization and muting
    """
    def __init__(self, data_format="ascii", misfit="waveform", backproject=None,
                 normalize=None, filter=None, min_period=None, max_period=None,
                 min_freq=None, max_freq=None, mute=None, path_preprocess=None,
                 **kwargs):
        """
        Preprocessing module parameters

        :type data_format: str
        :param data_format: data format for reading traces into memory. For
            available see: seisflows.plugins.preprocess.readers
        :type misfit: str
        :param misfit: misfit function for waveform comparisons. For available
            see seisflows.plugins.preprocess.misfit
        :type backproject: str
        :param backproject: backprojection function for migration, or the
            objective function in FWI. For available see
            seisflows.plugins.preprocess.adjoint
        :type normalize: str
        :param normalize: Data normalization parameters used to normalize the
            amplitudes of waveforms. Choose from two sets:
            ENORML1: normalize per event by L1 of traces; OR
            ENORML2: normalize per event by L2 of traces;
            &
            TNORML1: normalize per trace by L1 of itself; OR
            TNORML2: normalize per trace by L2 of itself
        :type filter: str
        :param filter: Data filtering type, available options are:
            BANDPASS (req. MIN/MAX PERIOD/FREQ);
            LOWPASS (req. MAX_FREQ or MIN_PERIOD);
            HIGHPASS (req. MIN_FREQ or MAX_PERIOD)
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
            EARLY: mute early arrivals;
            LATE: mute late arrivals;
            SHORT: mute short source-receiver distances;
            LONG: mute long source-receiver distances
        :type path_preprocess: str
        :param path_preprocess: scratch path for all preprocessing processes,
            including saving files
        """
        self.data_format = data_format.title()
        self.misfit = misfit
        self.backproject = backproject,
        self.normalize = normalize
        self.filter = filter
        self.min_period = min_period
        self.max_period = max_period
        self.min_freq = min_freq
        self.max_freq = max_freq
        self.mute = mute or []
        self.normalize = normalize or []
        self.path = path_preprocess or \
                    os.path.join(os.getcwd(), "scratch", "preprocess")

        # TODO: Add the mute parameters here, const, slope and dist

        self._misfit = None
        self._adjoint = None
        self._reader = None
        self._writer = None

    def check(self, validate=True):
        """ 
        Checks parameters and paths
        """
        # Data normalization option
        if self.normalize:
            acceptable_norms = {"TNORML1", "TNORML2", "ENORML1", "ENORML2"}
            chosen_norms = [_.upper() for _ in self.normalize]
            assert(set(chosen_norms).issubset(acceptable_norms))

        # Data muting options
        if self.mute:
            acceptable_mutes = {"EARLY", "LATE", "LONG", "SHORT"}
            chosen_mutes = [_.upper() for _ in self.mute]
            assert(set(chosen_mutes).issubset(acceptable_mutes))
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
            acceptable_filters = ["BANDPASS", "LOWPASS", "HIGHPASS"]
            assert self.filter.upper() in acceptable_filters, \
                f"self.filter must be in {acceptable_filters}"

            # Set the min/max frequencies and periods, frequency takes priority
            if self.min_freq is not None:
                self.max_period = 1 / self.min_freq
            elif self.max_period is not None:
                self.min_freq = 1 / self.max_period

            if self.max_freq is not None:
                self.min_period = 1 / self.max_freq
            elif self.min_period is not None:
                self.max_freq =  1 / self.min_period

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

        # Assert that readers and writers available
        # TODO | This is a bit vague as dir contains imported modules and hidden
        # TODO | variables (e.g., np, __name__)
        assert(self.data_format.lower() in dir(readers)), (
            f"Reader {self.data_format} not found")
        assert(self.data_format.lower() in dir(writers)), (
            f"Writer {self.data_format} not found")

    def setup(self):
        """
        Sets up data preprocessing machinery by dynamicalyl loading the
        misfit, adjoint source type, and specifying the expected file type
        for input and output seismic data.
        """
        unix.mkdir(self.path)

        # Define misfit function and adjoint trace generator
        if self.misfit:
            logger.debug(f"misfit function is: '{self.misfit}'")
            self._misfit = getattr(misfit, self.misfit.lower())
            self._adjoint = getattr(adjoint, self.misfit.lower())
        elif self.backproject:
            logger.debug(
                f"backproject function is: '{self.backproject}'"
            )
            self._adjoint = getattr(adjoint, self.backproject.lower())

        # Define seismic data reader and writer
        self._reader = getattr(readers, self.data_format.lower())
        self._writer = getattr(writers, self.data_format.lower())

    def finalize(self):
        """
        Any finalization processes that need to take place at the end of an
        iteration
        """
        pass

    def initialize_adjoint_traces(self, filenames=None):
        """
        TO DO
        """
        for filename in self.data_filenames:
            st = self.preprocess.reader(path=os.path.join(self.cwd, "traces", "obs"),
                        filename=filename
                        )
            # Zero out data just so we have empty adjoint traces as SPECFEM
            # will expect all adjoint sources to have all components
            st *= 0

            # Write traces back to the adjoint trace directory
            preprocess.writer(st=st, filename=filename,
                              path=os.path.join(self.cwd, "traces", "adj")
                              )

    def quantify_misfit(self, observed, synthetic,
                        write_residuals=False, write_adjsrcs=False, **kwargs):
        """
        Prepares solver for gradient evaluation by writing residuals and
        adjoint traces. Meant to be called by solver.eval_func().

        Reads in observed and synthetic waveforms, applies optional
        preprocessing, assesses misfit, and writes out adjoint sources and
        STATIONS_ADJOINT file.

        .. note::
            Meant to be called by solver.eval_func(), may have unused arguments
            to keep functions general across subclasses.

        :type cwd: str
        :param cwd: current specfem working directory containing observed and
            synthetic seismic data to be read and processed. Should be defined
            by solver.cwd
        :type filenames: list of str
        :param filenames: list of filenames defining the files in traces
        """
        for obs_fid, syn_fid in zip(observed, synthetic):
            obs = self._reader(filename=obs_fid)
            syn = self._reader(filename=syn_fid)

            # Process observations and synthetics identically
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
            for tr_obs, tr_syn in zip(obs, syn):
                if write_residuals:
                    residual = self._misfit(obs=tr_obs.data, syn=tr_syn.data,
                                            nt=tr_syn.stats.npts,
                                            dt=tr_syn.stats.delta)
                    with open(write_residuals, "a") as f:
                        f.write(f"{residual:.2E}\n")
                if write_adjsrcs:
                    adjsrc = syn.copy()
                    adjsrc.data = self._adjoint(obs=tr_obs.data, syn=tr_syn.data,
                                                nt=tr_syn.stats.npts,
                                                dt=tr_syn.stats.delta)
                    if self.data_format.upper() == "ASCII":
                        # Change the extension to '.adj' from whatever it is
                        ext = os.path.splitext(os.path.basename(obs))[-1]
                        filename = os.path.basename(obs).replace(ext, ".adj")
                    elif self.data_format.upper() == "SU":
                        # TODO implement this
                        raise NotImplementedError
                    self._writer(st=adjsrc,
                                 filename=os.path.join(write_adjsrcs, filename)
                                 )

    def prepare_eval_grad(self, cwd, taskid, filenames, **kwargs):
        """
        Prepares solver for gradient evaluation by writing residuals and
        adjoint traces. Meant to be called by solver.eval_func().

        Reads in observed and synthetic waveforms, applies optional
        preprocessing, assesses misfit, and writes out adjoint sources and
        STATIONS_ADJOINT file.

        .. note::
            Meant to be called by solver.eval_func(), may have unused arguments
            to keep functions general across subclasses.

        :type cwd: str
        :param cwd: current specfem working directory containing observed and
            synthetic seismic data to be read and processed. Should be defined
            by solver.cwd
        :type filenames: list of str
        :param filenames: list of filenames defining the files in traces
        """
        if taskid == 0:
            logger.debug("preparing files for gradient evaluation")

        for filename in filenames:
            obs = self._reader(path=os.path.join(cwd, "traces", "obs"),
                               filename=filename)
            syn = self._reader(path=os.path.join(cwd, "traces", "syn"),
                               filename=filename)

            # Process observations and synthetics identically
            if self.filter:
                if taskid == 0:
                    logger.debug(f"applying {self.filter} filter to data")
                obs = self._apply_filter(obs)
                syn = self._apply_filter(syn)
            if self.mute:
                if taskid == 0:
                    logger.debug(f"applying {self.mute} mutes to data")
                obs = self._apply_mute(obs)
                syn = self._apply_mute(syn)
            if self.normalize:
                if taskid == 0:
                    logger.debug(
                        f"normalizing data with: {self.normalize}"
                    )
                obs = self._apply_normalize(obs)
                syn = self._apply_normalize(syn)

            if self.misfit is not None:
                self._write_residuals(cwd, syn, obs)

            # Write the adjoint traces. Rename file extension for Specfem
            if self.data_format.upper() == "ASCII":
                # Change the extension to '.adj' from whatever it is
                ext = os.path.splitext(filename)[-1]
                filename_out = filename.replace(ext, ".adj")
            elif self.data_format.upper() == "SU":
                # TODO implement this
                raise NotImplementedError

            self._write_adjoint_traces(path=os.path.join(cwd, "traces", "adj"),
                                       syn=syn, obs=obs, filename=filename_out)

        # Copy over the STATIONS file to STATIONS_ADJOINT required by Specfem
        # ASSUMING that all stations are used in adjoint simulation
        # TODO !!! This is SPECFEM dependent? Belongs in solver.specfem?
        src = os.path.join(cwd, "DATA", "STATIONS")
        dst = os.path.join(cwd, "DATA", "STATIONS_ADJOINT")
        unix.cp(src, dst)

    def sum_residuals(self, files):
        """
        Sums squares of residuals

        :type files: str
        :param files: list of single-column text files containing residuals
        :rtype: float
        :return: sum of squares of residuals
        """
        total_misfit = 0.
        for filename in files:
            total_misfit += np.sum(np.loadtxt(filename) ** 2.)

        return total_misfit

    def _write_residuals(self, obs, syn, output):
        """
        Computes residuals between observed and synthetic seismogram based on
        the misfit function self.misfit. Saves the residuals for each
        data-synthetic pair into a text file located at:

        ./scratch/solver/*/residuals

        The resulting file will be a single-column ASCII file that needs to be
        summed before use by the solver

        :type path: str
        :param path: location "adjoint traces" will be written
        :type syn: obspy.core.stream.Stream
        :param syn: synthetic data
        :type obs: obspy.core.stream.Stream
        :param syn: observed data
        """
        residuals = []
        for tr_obs, tr_syn in zip(obs, syn):
            residual = self._misfit(syn=tr_syn.data, obs=tr_obs.data,
                                    nt=tr_syn.stats.npts,
                                    dt=tr_syn.stats.delta
                                    )
            residuals.append(residual)

        filename = os.path.join(output, "residuals")
        if os.path.exists(filename):
            residuals = np.append(residuals, np.loadtxt(filename))

        np.savetxt(filename, residuals)

    def _write_adjoint_traces(self, path, syn, obs, filename):
        """
        Writes "adjoint traces" required for gradient computation

        :type path: str
        :param path: location "adjoint traces" will be written
        :type syn: obspy.core.stream.Stream
        :param syn: synthetic data
        :type obs: obspy.core.stream.Stream
        :param syn: observed data
        :type filename: str
        :param filename: filename to write adjoint traces to
        """
        # Use the synthetics as a template for the adjoint sources
        adj = syn.copy()
        for tr_adj, tr_obs, tr_syn in zip(adj, obs, syn):
            tr_adj.data = self._adjoint(syn=tr_syn.data, obs=tr_obs.data,
                                        nt=tr_syn.stats.npts,
                                        dt=tr_syn.stats.delta
                                        )

        self._writer(adj, path, filename)

    def _apply_filter(self, st):
        """
        Apply a filter to waveform data using ObsPy

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
                                      const=self.early_const,
                                      choice="EARLY")
        if "LATE" in mute_choices:
            st = signal.mute_arrivals(st, slope=self.late_slope,
                                      const=self.late_const,
                                      choice="LATE")
        if "SHORT" in mute_choices:
            st = signal.mute_offsets(st, dist=self.short_dist,
                                     choice="SHORT")
        if "LONG" in mute_choices:
            st = signal.mute_arrivals(st, dist=self.long_dist,
                                      choice="LONG")

        return st

    def _apply_normalize(self, st):
        """
        Normalize the amplitudes of waveforms based on user choice

        .. note::
            The normalization function has been refactored but not tested
            as I was not aware of the intended functionality. Not gauranteed
            to work, use at your own risk.

        :type st: obspy.core.stream.Stream
        :param st: All of the data streams to be normalized
        :rtype: obspy.core.stream.Stream
        :return: stream with normalized traces
        """
        st_out = st.copy()
        norm_choices = [_.upper() for _ in self.normalize]

        # Normalize an event by the L1 norm of all traces
        if 'ENORML1' in norm_choices:
            w = 0.
            for tr in st_out:
                w += np.linalg.norm(tr.data, ord=1)
            for tr in st_out:
                tr.data /= w
        # Normalize an event by the L2 norm of all traces
        elif "ENORML2" in norm_choices:
            w = 0.
            for tr in st_out:
                w += np.linalg.norm(tr.data, ord=2)
            for tr in st_out:
                tr.data /= w
        # Normalize each trace by its L1 norm
        if "TNORML1" in norm_choices:
            for tr in st_out:
                w = np.linalg.norm(tr.data, ord=1)
                if w > 0:
                    tr.data /= w
        elif "TNORML2" in norm_choices:
            # normalize each trace by its L2 norm
            for tr in st_out:
                w = np.linalg.norm(tr.data, ord=2)
                if w > 0:
                    tr.data /= w

        return st_out
