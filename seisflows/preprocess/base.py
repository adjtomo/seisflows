#!/usr/bin/env python3
"""
The SeisFlows Preprocessing module is in charge of interacting with seismic
data (observed and synthetic). It should contain functionality to read and write
seismic data, apply preprocessing such as filtering, quantify misfit,
and write adjoint sources that are expected by the solver.
"""
import os
import sys
import obspy
import logging
import numpy as np

from seisflows.tools import msg
from seisflows.tools import signal, unix
from seisflows.plugins.preprocess import adjoint, misfit, readers, writers
from seisflows.config import SeisFlowsPathsParameters

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]


class Default:
    """
    Default SeisFlows preprocessing class

    Provides data processing functions for seismic traces, with options for
    data misfit, filtering, normalization and muting
    """
    # Class-specific logger accessed using self.logger
    logger = logging.getLogger(__name__).getChild(__qualname__)

    def __init__(self):
        """
        These parameters should not be set by __init__!
        Attributes are just initialized as NoneTypes for clarity and docstrings
        """
        self.misfit = None
        self.adjoint = None
        self.reader = None
        self.writer = None

    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters()

        # Define the Parameters required by this module
        sf.par("MISFIT", required=False, default="waveform", par_type=str,
               docstr="Misfit function for waveform comparisons, for available "
                      "see seisflows.plugins.misfit")

        sf.par("BACKPROJECT", required=False, default="null", par_type=str,
               docstr="Backprojection function for migration, for available "
                      "see seisflows.plugins.adjoint")

        sf.par("NORMALIZE", required=False, default="null", par_type=str,
               docstr="Data normalization option")

        sf.par("FILTER", required=False, default="null", par_type=str,
               docstr="Data filtering type, available options are:"
                      "BANDPASS (req. MIN/MAX PERIOD/FREQ);"
                      "LOWPASS (req. MAX_FREQ or MIN_PERIOD); "
                      "HIGHPASS (req. MIN_FREQ or MAX_PERIOD) ")

        sf.par("MIN_PERIOD", required=False, par_type=float,
               docstr="Minimum filter period applied to time series."
                      "See also MIN_FREQ, MAX_FREQ, if User defines FREQ "
                      "parameters, they will overwrite PERIOD parameters.")

        sf.par("MAX_PERIOD", required=False, par_type=float,
               docstr="Maximum filter period applied to time series."
                      "See also MIN_FREQ, MAX_FREQ, if User defines FREQ "
                      "parameters, they will overwrite PERIOD parameters.")

        sf.par("MIN_FREQ", required=False, par_type=float,
               docstr="Maximum filter frequency applied to time series."
                      "See also MIN_PERIOD, MAX_PERIOD, if User defines FREQ "
                      "parameters, they will overwrite PERIOD parameters.")

        sf.par("MAX_FREQ", required=False, par_type=float,
               docstr="Maximum filter frequency applied to time series,"
                      "See also MIN_PERIOD, MAX_PERIOD, if User defines FREQ "
                      "parameters, they will overwrite PERIOD parameters.")

        sf.par("MUTE", required=False, par_type=list, default=[],
               docstr="Data mute parameters used to zero out early / late "
                      "arrivals or offsets. Choose any number of: "
                      "EARLY: mute early arrivals; "
                      "LATE: mute late arrivals; "
                      "SHORT: mute short source-receiver distances; "
                      "LONG: mute long source-receiver distances")
        sf.par("NORMALIZE", required=False, par_type=list, default=[],
               docstr="Data normalization parameters used to normalize the "
                      "amplitudes of waveforms. Choose from two sets: "
                      "ENORML1: normalize per event by L1 of traces; OR "
                      "ENORML2: normalize per event by L2 of traces; AND "
                      "TNORML1: normalize per trace by L1 of itself; OR "
                      "TNORML2: normalize per trace by L2 of itself")

        # TODO: Add the mute parameters here, const, slope and dist

        # Define the Paths required by this module
        sf.path("PREPROCESS", required=False,
                default=os.path.join(PATH.SCRATCH, "preprocess"),
                docstr="scratch path to store any preprocessing outputs")

        return sf

    def check(self, validate=True):
        """ 
        Checks parameters and paths
        """
        if validate:
            self.required.validate()

        # Data normalization option
        if PAR.NORMALIZE:
            acceptable_norms = {"TNORML1", "TNORML2", "ENORML1", "ENORML2"}
            chosen_norms = [_.upper() for _ in PAR.NORMALIZE]
            assert(set(chosen_norms).issubset(acceptable_norms))

        # Data muting options
        if PAR.MUTE:
            acceptable_mutes = {"EARLY", "LATE", "LONG", "SHORT"}
            chosen_mutes = [_.upper() for _ in PAR.MUTE]
            assert(set(chosen_mutes).issubset(acceptable_mutes))
            if "EARLY" in chosen_mutes:
                assert(PAR.EARLY_SLOPE is not None)
                assert(PAR.EARLY_SLOPE >= 0.)
                assert(PAR.EARLY_CONST is not None)
            if "LATE" in chosen_mutes:
                assert(PAR.LATE_SLOPE is not None)
                assert(PAR.LATE_SLOPE >= 0.)
                assert(PAR.LATE_CONST is not None)
            if "SHORT" in chosen_mutes:
                assert(PAR.SHORT_DIST is not None)
                assert (PAR.SHORT_DIST > 0)
            if "LONG" in chosen_mutes:
                assert(PAR.LONG_DIST is not None)
                assert (PAR.LONG_DIST > 0)

        # Data filtering options that will be passed to ObsPy filters
        if PAR.FILTER:
            acceptable_filters = ["BANDPASS", "LOWPASS", "HIGHPASS"]
            assert PAR.FILTER.upper() in acceptable_filters, \
                f"PAR.FILTER must be in {acceptable_filters}"

            # Set the min/max frequencies and periods, frequency takes priority
            if PAR.MIN_FREQ is not None:
                PAR.MAX_PERIOD = 1 / PAR.MIN_FREQ
            elif PAR.MAX_PERIOD is not None:
                PAR.MIN_FREQ = 1 / PAR.MAX_PERIOD

            if PAR.MAX_FREQ is not None:
                PAR.MIN_PERIOD = 1 / PAR.MAX_FREQ
            elif PAR.MIN_PERIOD is not None:
                PAR.MAX_FREQ =  1 / PAR.MIN_PERIOD

            # Check that the correct filter bounds have been set
            if PAR.FILTER.upper() == "BANDPASS":
                assert(PAR.MIN_FREQ is not None and PAR.MAX_FREQ is not None), \
                    ("BANDPASS filter PAR.MIN_PERIOD and PAR.MAX_PERIOD or " 
                     "PAR.MIN_FREQ and PAR.MAX_FREQ")
            elif PAR.FILTER.upper() == "LOWPASS":
                assert(PAR.MAX_FREQ is not None or PAR.MIN_PERIOD is not None),\
                    "LOWPASS requires PAR.MAX_FREQ or PAR.MIN_PERIOD"
            elif PAR.FILTER.upper() == "HIGHPASS":
                assert(PAR.MIN_FREQ is not None or PAR.MAX_PERIOD is not None),\
                    "HIGHPASS requires PAR.MIN_FREQ or PAR.MAX_PERIOD"

            # Check that filter bounds make sense, by this point, MIN and MAX
            # FREQ and PERIOD should be set, so we just check the FREQ
            assert(0 < PAR.MIN_FREQ < np.inf), "0 < PAR.MIN_FREQ < inf"
            assert(0 < PAR.MAX_FREQ < np.inf), "0 < PAR.MAX_FREQ < inf"
            assert(PAR.MIN_FREQ < PAR.MAX_FREQ), "PAR.MIN_FREQ < PAR.MAX_FREQ"

        # Assert that readers and writers available
        # TODO | This is a bit vague as dir contains imported modules and hidden
        # TODO | variables (e.g., np, __name__)
        assert(PAR.FORMAT.lower() in dir(readers)), (
            f"Reader {PAR.FORMAT} not found")
        assert(PAR.FORMAT.lower() in dir(writers)), (
            f"Writer {PAR.FORMAT} not found")

        # Assert that either misfit or backproject exists 
        if PAR.WORKFLOW.upper() == "INVERSION":
            assert(PAR.MISFIT is not None)

    def setup(self):
        """
        Sets up data preprocessing machinery by dynamicalyl loading the
        misfit, adjoint source type, and specifying the expected file type
        for input and output seismic data.
        """
        unix.mkdir(PATH.PREPROCESS)

        # Define misfit function and adjoint trace generator
        if PAR.MISFIT:
            self.logger.debug(f"misfit function is: '{PAR.MISFIT}'")
            self.misfit = getattr(misfit, PAR.MISFIT.lower())
            self.adjoint = getattr(adjoint, PAR.MISFIT.lower())
        elif PAR.BACKPROJECT:
            self.logger.debug(f"backproject function is: '{PAR.BACKPROJECT}'")
            self.adjoint = getattr(adjoint, PAR.BACKPROJECT.lower())

        # Define seismic data reader and writer
        self.reader = getattr(readers, PAR.FORMAT.lower())
        self.writer = getattr(writers, PAR.FORMAT.lower())

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
            if PAR.FILTER:
                if taskid == 0:
                    self.logger.debug(f"applying {PAR.FILTER} filter to data")
                obs = self._apply_filter(obs)
                syn = self._apply_filter(syn)
            if PAR.MUTE:
                if taskid == 0:
                    self.logger.debug(f"applying {PAR.MUTE} mutes to data")
                obs = self._apply_mute(obs)
                syn = self._apply_mute(syn)
            if PAR.NORMALIZE:
                if taskid == 0:
                    self.logger.debug(f"normalizing data with: {PAR.NORMALIZE}")
                obs = self._apply_normalize(obs)
                syn = self._apply_normalize(syn)

            if PAR.MISFIT is not None:
                self._write_residuals(cwd, syn, obs)

            # Write the adjoint traces. Rename file extension for Specfem
            if PAR.FORMAT.upper() == "ASCII":
                # Change the extension to '.adj' from whatever it is
                ext = os.path.splitext(filename)[-1]
                filename_out = filename.replace(ext, ".adj")
            elif PAR.FORMAT.upper() == "SU":
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

    def finalize(self):
        """
        Any finalization processes that need to take place at the end of an
        iteration
        """
        pass

    def _write_residuals(self, path, syn, obs):
        """
        Computes residuals between observed and synthetic seismogram based on
        the misfit function PAR.MISFIT. Saves the residuals for each
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

        if PAR.FILTER.upper() == "BANDPASS":
            st.filter("bandpass", zerophase=True, freqmin=PAR.MIN_FREQ,
                      freqmax=PAR.FREQMAX)
        elif PAR.FILTER.upper() == "LOWPASS":
            st.filter("lowpass", zerophase=True, freq=PAR.MAX_FREQ)
        elif PAR.FILTER.upper() == "HIGHPASS":
            st.filter("highpass", zerophase=True, freq=PAR.MIN_FREQ)

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
        mute_choices = [_.upper() for _ in PAR.MUTE]
        if "EARLY" in mute_choices:
            st = signal.mute_arrivals(st, slope=PAR.EARLY_SLOPE,
                                      const=PAR.EARLY_CONST,
                                      choice="EARLY")
        if "LATE" in mute_choices:
            st = signal.mute_arrivals(st, slope=PAR.LATE_SLOPE,
                                      const=PAR.LATE_CONST,
                                      choice="LATE")
        if "SHORT" in mute_choices:
            st = signal.mute_offsets(st, dist=PAR.SHORT_DIST,
                                     choice="SHORT")
        if "LONG" in mute_choices:
            st = signal.mute_arrivals(st, dist=PAR.LONG_DIST,
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
        norm_choices = [_.upper() for _ in PAR.NORMALIZE]

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
