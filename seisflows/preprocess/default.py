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

from seisflows.core import Base
from seisflows.tools import signal, unix
from seisflows.plugins.preprocess import adjoint, misfit, readers, writers


class Default(Base):
    """
    Default SeisFlows preprocessing class

    Provides data processing functions for seismic traces, with options for
    data misfit, filtering, normalization and muting
    """
    def __init__(self):
        """
        These parameters should not be set by __init__!
        Attributes are just initialized as NoneTypes for clarity and docstrings
        """
        super().__init__()

        # Define the Parameters required by this module
        self.required.par(
            "MISFIT", required=False, default="waveform", par_type=str,
            docstr="Misfit function for waveform comparisons, for available "
                   "see seisflows.plugins.misfit"
        )
        self.required.par(
            "BACKPROJECT", required=False, default="null", par_type=str,
            docstr="Backprojection function for migration, for available see "
                   "seisflows.plugins.adjoint"
        )
        self.required.par(
            "NORMALIZE", required=False, default="null", par_type=str,
            docstr="Data normalization option"
        )
        self.required.par(
            "FILTER", required=False, default="null", par_type=str,
            docstr="Data filtering type, available options are:"
                   "BANDPASS (req. MIN/MAX PERIOD/FREQ);"
                   "LOWPASS (req. MAX_FREQ or MIN_PERIOD); "
                   "HIGHPASS (req. MIN_FREQ or MAX_PERIOD) "
        )
        self.required.par(
            "MIN_PERIOD", required=False, par_type=float,
            docstr="Minimum filter period applied to time series."
                   "See also MIN_FREQ, MAX_FREQ, if User defines FREQ "
                   "parameters, they will overwrite PERIOD parameters."
        )
        self.required.par(
            "MAX_PERIOD", required=False, par_type=float,
            docstr="Maximum filter period applied to time series."
                   "See also MIN_FREQ, MAX_FREQ, if User defines FREQ "
                   "parameters, they will overwrite PERIOD parameters."
        )
        self.required.par(
            "MIN_FREQ", required=False, par_type=float,
            docstr="Maximum filter frequency applied to time series."
                   "See also MIN_PERIOD, MAX_PERIOD, if User defines FREQ "
                   "parameters, they will overwrite PERIOD parameters."
        )
        self.required.par(
            "MAX_FREQ", required=False, par_type=float,
            docstr="Maximum filter frequency applied to time series,"
                   "See also MIN_PERIOD, MAX_PERIOD, if User defines FREQ "
                   "parameters, they will overwrite PERIOD parameters."
        )
        self.required.par(
            "MUTE", required=False, par_type=list, default=[],
            docstr="Data mute parameters used to zero out early / late "
                   "arrivals or offsets. Choose any number of: "
                   "EARLY: mute early arrivals; "
                   "LATE: mute late arrivals; "
                   "SHORT: mute short source-receiver distances; "
                   "LONG: mute long source-receiver distances"
        )
        self.required.par(
            "NORMALIZE", required=False, par_type=list, default=[],
            docstr="Data normalization parameters used to normalize the "
                   "amplitudes of waveforms. Choose from two sets: "
                   "ENORML1: normalize per event by L1 of traces; OR "
                   "ENORML2: normalize per event by L2 of traces; AND "
                   "TNORML1: normalize per trace by L1 of itself; OR "
                   "TNORML2: normalize per trace by L2 of itself"
        )
        # TODO: Add the mute parameters here, const, slope and dist

        # Define the Paths required by this module
        self.required.path(
            "PREPROCESS", required=False,
            default=os.path.join(self.path.WORKDIR, "scratch", "preprocess"),
            docstr="scratch path to store any preprocessing outputs"
        )

        self.misfit = None
        self.adjoint = None
        self.reader = None
        self.writer = None

    def check(self, validate=True):
        """ 
        Checks parameters and paths
        """
        super().check(validate=validate)

        # Data normalization option
        if self.par.NORMALIZE:
            acceptable_norms = {"TNORML1", "TNORML2", "ENORML1", "ENORML2"}
            chosen_norms = [_.upper() for _ in self.par.NORMALIZE]
            assert(set(chosen_norms).issubset(acceptable_norms))

        # Data muting options
        if self.par.MUTE:
            acceptable_mutes = {"EARLY", "LATE", "LONG", "SHORT"}
            chosen_mutes = [_.upper() for _ in self.par.MUTE]
            assert(set(chosen_mutes).issubset(acceptable_mutes))
            if "EARLY" in chosen_mutes:
                assert(self.par.EARLY_SLOPE is not None)
                assert(self.par.EARLY_SLOPE >= 0.)
                assert(self.par.EARLY_CONST is not None)
            if "LATE" in chosen_mutes:
                assert(self.par.LATE_SLOPE is not None)
                assert(self.par.LATE_SLOPE >= 0.)
                assert(self.par.LATE_CONST is not None)
            if "SHORT" in chosen_mutes:
                assert(self.par.SHORT_DIST is not None)
                assert (self.par.SHORT_DIST > 0)
            if "LONG" in chosen_mutes:
                assert(self.par.LONG_DIST is not None)
                assert (self.par.LONG_DIST > 0)

        # Data filtering options that will be passed to ObsPy filters
        if self.par.FILTER:
            acceptable_filters = ["BANDPASS", "LOWPASS", "HIGHPASS"]
            assert self.par.FILTER.upper() in acceptable_filters, \
                f"self.par.FILTER must be in {acceptable_filters}"

            # Set the min/max frequencies and periods, frequency takes priority
            if self.par.MIN_FREQ is not None:
                self.par.MAX_PERIOD = 1 / self.par.MIN_FREQ
            elif self.par.MAX_PERIOD is not None:
                self.par.MIN_FREQ = 1 / self.par.MAX_PERIOD

            if self.par.MAX_FREQ is not None:
                self.par.MIN_PERIOD = 1 / self.par.MAX_FREQ
            elif self.par.MIN_PERIOD is not None:
                self.par.MAX_FREQ =  1 / self.par.MIN_PERIOD

            # Check that the correct filter bounds have been set
            if self.par.FILTER.upper() == "BANDPASS":
                assert(self.par.MIN_FREQ is not None and
                       self.par.MAX_FREQ is not None), \
                    ("BANDPASS filter PAR.MIN_PERIOD and PAR.MAX_PERIOD or " 
                     "PAR.MIN_FREQ and PAR.MAX_FREQ")
            elif self.par.FILTER.upper() == "LOWPASS":
                assert(self.par.MAX_FREQ is not None or
                       self.par.MIN_PERIOD is not None),\
                    "LOWPASS requires PAR.MAX_FREQ or PAR.MIN_PERIOD"
            elif self.par.FILTER.upper() == "HIGHPASS":
                assert(self.par.MIN_FREQ is not None or
                       self.par.MAX_PERIOD is not None),\
                    "HIGHPASS requires PAR.MIN_FREQ or PAR.MAX_PERIOD"

            # Check that filter bounds make sense, by this point, MIN and MAX
            # FREQ and PERIOD should be set, so we just check the FREQ
            assert(0 < self.par.MIN_FREQ < np.inf), "0 < PAR.MIN_FREQ < inf"
            assert(0 < self.par.MAX_FREQ < np.inf), "0 < PAR.MAX_FREQ < inf"
            assert(self.par.MIN_FREQ < self.par.MAX_FREQ), (
                "PAR.MIN_FREQ < PAR.MAX_FREQ"
            )

        # Assert that readers and writers available
        # TODO | This is a bit vague as dir contains imported modules and hidden
        # TODO | variables (e.g., np, __name__)
        assert(self.par.FORMAT.lower() in dir(readers)), (
            f"Reader {self.par.FORMAT} not found")
        assert(self.par.FORMAT.lower() in dir(writers)), (
            f"Writer {self.par.FORMAT} not found")

        # Assert that either misfit or backproject exists
        if self.par.WORKFLOW.upper() == "INVERSION":
            assert(self.par.MISFIT is not None)

    def setup(self):
        """
        Sets up data preprocessing machinery by dynamicalyl loading the
        misfit, adjoint source type, and specifying the expected file type
        for input and output seismic data.
        """
        unix.mkdir(self.path.PREPROCESS)

        # Define misfit function and adjoint trace generator
        if self.par.MISFIT:
            self.logger.debug(f"misfit function is: '{self.par.MISFIT}'")
            self.misfit = getattr(misfit, self.par.MISFIT.lower())
            self.adjoint = getattr(adjoint, self.par.MISFIT.lower())
        elif self.par.BACKPROJECT:
            self.logger.debug(
                f"backproject function is: '{self.par.BACKPROJECT}'"
            )
            self.adjoint = getattr(adjoint, self.par.BACKPROJECT.lower())

        # Define seismic data reader and writer
        self.reader = getattr(readers, self.par.FORMAT.lower())
        self.writer = getattr(writers, self.par.FORMAT.lower())

    def finalize(self):
        """
        Any finalization processes that need to take place at the end of an
        iteration
        """
        super().finalize()

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
            self.logger.debug("preparing files for gradient evaluation")

        for filename in filenames:
            obs = self.reader(path=os.path.join(cwd, "traces", "obs"),
                              filename=filename)
            syn = self.reader(path=os.path.join(cwd, "traces", "syn"),
                              filename=filename)

            # Process observations and synthetics identically
            if self.par.FILTER:
                if taskid == 0:
                    self.logger.debug(f"applying {self.par.FILTER} filter to data")
                obs = self._apply_filter(obs)
                syn = self._apply_filter(syn)
            if self.par.MUTE:
                if taskid == 0:
                    self.logger.debug(f"applying {self.par.MUTE} mutes to data")
                obs = self._apply_mute(obs)
                syn = self._apply_mute(syn)
            if self.par.NORMALIZE:
                if taskid == 0:
                    self.logger.debug(
                        f"normalizing data with: {self.par.NORMALIZE}"
                    )
                obs = self._apply_normalize(obs)
                syn = self._apply_normalize(syn)

            if self.par.MISFIT is not None:
                self._write_residuals(cwd, syn, obs)

            # Write the adjoint traces. Rename file extension for Specfem
            if self.par.FORMAT.upper() == "ASCII":
                # Change the extension to '.adj' from whatever it is
                ext = os.path.splitext(filename)[-1]
                filename_out = filename.replace(ext, ".adj")
            elif self.par.FORMAT.upper() == "SU":
                # TODO implement this
                raise NotImplementedError

            self._write_adjoint_traces(path=os.path.join(cwd, "traces", "adj"),
                                       syn=syn, obs=obs, filename=filename_out)

        # Copy over the STATIONS file to STATIONS_ADJOINT required by Specfem
        # ASSUMING that all stations are used in adjoint simulation
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

    def _write_residuals(self, path, syn, obs):
        """
        Computes residuals between observed and synthetic seismogram based on
        the misfit function self.par.MISFIT. Saves the residuals for each
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
            residual = self.misfit(syn=tr_syn.data, obs=tr_obs.data,
                                   nt=tr_syn.stats.npts,
                                   dt=tr_syn.stats.delta
                                   )
            residuals.append(residual)

        filename = os.path.join(path, "residuals")
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
            tr_adj.data = self.adjoint(syn=tr_syn.data, obs=tr_obs.data,
                                       nt=tr_syn.stats.npts,
                                       dt=tr_syn.stats.delta
                                       )

        self.writer(adj, path, filename)

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

        if self.par.FILTER.upper() == "BANDPASS":
            st.filter("bandpass", zerophase=True, freqmin=self.par.MIN_FREQ,
                      freqmax=self.par.FREQMAX)
        elif self.par.FILTER.upper() == "LOWPASS":
            st.filter("lowpass", zerophase=True, freq=self.par.MAX_FREQ)
        elif self.par.FILTER.upper() == "HIGHPASS":
            st.filter("highpass", zerophase=True, freq=self.par.MIN_FREQ)

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
        mute_choices = [_.upper() for _ in self.par.MUTE]
        if "EARLY" in mute_choices:
            st = signal.mute_arrivals(st, slope=self.par.EARLY_SLOPE,
                                      const=self.par.EARLY_CONST,
                                      choice="EARLY")
        if "LATE" in mute_choices:
            st = signal.mute_arrivals(st, slope=self.par.LATE_SLOPE,
                                      const=self.par.LATE_CONST,
                                      choice="LATE")
        if "SHORT" in mute_choices:
            st = signal.mute_offsets(st, dist=self.par.SHORT_DIST,
                                     choice="SHORT")
        if "LONG" in mute_choices:
            st = signal.mute_arrivals(st, dist=self.par.LONG_DIST,
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
        norm_choices = [_.upper() for _ in self.par.NORMALIZE]

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
