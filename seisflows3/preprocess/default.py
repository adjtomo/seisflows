#!/usr/bin/env python
"""
This is the main seisflows.preprocess.base

This is a main Seisflows class, it controls the preprocessing.
"""
import os
import sys
import obspy
import numpy as np

from seisflows3.tools import msg
from seisflows3.tools import signal, unix
from seisflows3.config import custom_import
from seisflows3.tools.err import ParameterError
from seisflows3.tools.tools import exists, getset
from seisflows3.plugins import adjoint, misfit, readers, writers
from seisflows3.config import SeisFlowsPathsParameters

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]


class Default(custom_import("preprocess", "base")):
    """
    Default SeisFlows preprocessing class

    Provides data processing functions for seismic traces, with options for
    data misfit, filtering, normalization and muting
    """
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

        sf.par("MUTE", required=False, default="null", par_type=str,
               docstr="Data muting option")

        sf.par("FILTER", required=False, default="null", par_type=str,
               docstr="Data filtering option")

        return sf

    def check(self, validate=True):
        """ 
        Checks parameters and paths
        """
        if validate:
            self.required.validate()

        # Data normalization option
        if PAR.NORMALIZE:
            self.check_normalize_parameters()

        # Data muting option
        if PAR.MUTE:
            self.check_mute_parameters()

        # Data filtering option using Obspy
        if PAR.FILTER:
            self.check_filter_parameters()

        # Assert that readers and writers available
        if PAR.FORMAT not in dir(readers):
            print(msg.ReaderError)
            raise ParameterError()

        if PAR.FORMAT not in dir(writers):
            print(msg.WriterError)
            raise ParameterError()

        # Assert that either misfit or backproject exists 
        if PAR.WORKFLOW.upper() == "INVERSION" and not PAR.MISFIT:
            # !!! Need a better error here
            raise ParameterError("PAR.MISFIT must be set w/ default preprocess")

    def setup(self):
        """
        Sets up data preprocessing machinery
        """
        # Define misfit function and adjoint trace generator
        if PAR.MISFIT:
            self.misfit = getattr(misfit, PAR.MISFIT.lower())
            self.adjoint = getattr(adjoint, PAR.MISFIT.lower())
        elif PAR.BACKPROJECT:
            self.adjoint = getattr(adjoint, PAR.BACKPROJECT.lower())

        # Define seismic data reader and writer
        self.reader = getattr(readers, PAR.FORMAT)
        self.writer = getattr(writers, PAR.FORMAT)

    def prepare_eval_grad(self, cwd="./", **kwargs):
        """
        Prepares solver for gradient evaluation by writing residuals and
        adjoint traces

        :type cwd: str
        :param cwd: current specfem working directory containing observed and 
            synthetic seismic data to be read and processed
        """
        # Need to load solver mid-workflow as preprocess is loaded first
        solver = sys.modules["seisflows_solver"]

        for filename in solver.data_filenames:
            obs = self.reader(path=os.path.join(cwd, "traces", "obs"),
                              filename=filename)
            syn = self.reader(path=os.path.join(cwd, "traces", "syn"), 
                              filename=filename)

            # Process observations
            obs = self.apply_filter(obs)
            obs = self.apply_mute(obs)
            obs = self.apply_normalize(obs)

            # Process synthetics
            syn = self.apply_filter(syn)
            syn = self.apply_mute(syn)
            syn = self.apply_normalize(syn)

            if PAR.MISFIT:
                self.write_residuals(cwd, syn, obs)

            # Write the adjoint traces. Rename file extension for Specfem
            if PAR.FORMAT.upper() == "ASCII":
                # Change the extension to '.adj' from whatever it is
                ext = os.path.splitext(filename)[-1]
                filename_out = filename.replace(ext, ".adj")

            self.write_adjoint_traces(path=os.path.join(cwd, "traces", "adj"),
                                      syn=syn, obs=obs, filename=filename_out)

        # Copy over the STATIONS file to STATIONS_ADJOINT required by Specfem
        src = os.path.join(cwd, "DATA", "STATIONS")
        dst = os.path.join(cwd, "DATA", "STATIONS_ADJOINT")
        unix.cp(src, dst)
        

    def write_residuals(self, path, syn, obs):
        """
        Computes residuals

        :type path: str
        :param path: location "adjoint traces" will be written
        :type syn: obspy.core.stream.Stream
        :param syn: synthetic data
        :type obs: obspy.core.stream.Stream
        :param syn: observed data
        """
        nt, dt, _ = self.get_time_scheme(syn)
        nn, _ = self.get_network_size(syn)

        residuals = []
        for ii in range(nn):
            residuals.append(self.misfit(syn[ii].data, obs[ii].data, nt, dt))

        filename = os.path.join(path, "residuals")
        if exists(filename):
            residuals = np.append(residuals, np.loadtxt(filename))

        np.savetxt(filename, residuals)

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

    def write_adjoint_traces(self, path, syn, obs, filename):
        """
        Writes "adjoint traces" required for gradient computation

        :type path: str
        :param path: location "adjoint traces" will be written
        :type syn: obspy.core.stream.Stream
        :param syn: synthetic data
        :type obs: obspy.core.stream.Stream
        :param syn: observed data
        :type channel: str
        :param channel: channel or component code used by writer
        """
        nt, dt, _ = self.get_time_scheme(syn)
        nn, _ = self.get_network_size(syn)

        adj = syn
        for ii in range(nn):
            adj[ii].data = self.adjoint(syn[ii].data, obs[ii].data, nt, dt)

        self.writer(adj, path, filename)

    def apply_filter(self, st):
        """
        Apply a filter using Obspy

        :type st: obspy.core.stream.Stream
        :param st: stream to be filtered
        :rtype: obspy.core.stream.Stream
        :return: filtered traces
        """
        # If no filter given, don't do anything
        if PAR.FILTER is None:
            return st

        # Pre-processing before filtering
        st.detrend("demean")
        st.detrend("linear")
        st.taper(0.05, type="hann")

        if PAR.FILTER.upper() == "BANDPASS":
            st.filter("bandpass", zerophase=True,
                      freqmin=PAR.FREQMIN, freqmax=PAR.FREQMAX)
        elif PAR.FILTER.upper() == "LOWPASS":
            st.filter("lowpass", zerophase=True, freq=PAR.FREQ)
        elif PAR.FILTER.upper() == "HIGHPASS":
            st.filter("highpass", zerophase=True, freq=PAR.FREQ)

        return st

    def apply_mute(self, st):
        """
        Apply mute on data

        :type st: obspy.core.stream.Stream
        :param st: stream to mute
        :return:
        """
        if not PAR.MUTE:
            return st

        if 'MuteEarlyArrivals' in PAR.MUTE:
            traces = signal.mute_early_arrivals(
                st,
                PAR.MUTE_EARLY_ARRIVALS_SLOPE,  # (units: time/distance)
                PAR.MUTE_EARLY_ARRIVALS_CONST,  # (units: time)
                self.get_time_scheme(st),
                self.get_source_coords(st),
                self.get_receiver_coords(st))

        if 'MuteLateArrivals' in PAR.MUTE:
            traces = signal.mute_late_arrivals(
                st,
                PAR.MUTE_LATE_ARRIVALS_SLOPE,  # (units: time/distance)
                PAR.MUTE_LATE_ARRIVALS_CONST,  # (units: time)
                self.get_time_scheme(st),
                self.get_source_coords(st),
                self.get_receiver_coords(st))

        if 'MuteShortOffsets' in PAR.MUTE:
            traces = signal.mute_short_offsets(
                st,
                PAR.MUTE_SHORT_OFFSETS_DIST,
                self.get_source_coords(st),
                self.get_receiver_coords(st))

        if 'MuteLongOffsets' in PAR.MUTE:
            traces = signal.mute_long_offsets(
                st,
                PAR.MUTE_LONG_OFFSETS_DIST,
                self.get_source_coords(st),
                self.get_receiver_coords(st))

        return traces

    def apply_normalize(self, traces):
        """

        :param traces:
        :return:
        """
        if not PAR.NORMALIZE:
            return traces

        if 'NormalizeEventsL1' in PAR.NORMALIZE:
            # normalize event by L1 norm of all traces
            w = 0.
            for tr in traces:
                w += np.linalg.norm(tr.data, ord=1)
            for tr in traces:
                tr.data /= w

        elif 'NormalizeEventsL2' in PAR.NORMALIZE:
            # normalize event by L2 norm of all traces
            w = 0.
            for tr in traces:
                w += np.linalg.norm(tr.data, ord=2)
            for tr in traces:
                tr.data /= w

        if 'NormalizeTracesL1' in PAR.NORMALIZE:
            # normalize each trace by its L1 norm
            for tr in traces:
                w = np.linalg.norm(tr.data, ord=1)
                if w > 0:
                    tr.data /= w

        elif 'NormalizeTracesL2' in PAR.NORMALIZE:
            # normalize each trace by its L2 norm
            for tr in traces:
                w = np.linalg.norm(tr.data, ord=2)
                if w > 0:
                    tr.data /= w

        return traces

    def apply_filter_backwards(self, traces):
        """

        :param traces:
        :return:
        """
        for tr in traces:
            tr.data = np.flip(tr.data)

        traces = self.apply_filter()

        for tr in traces:
            tr.data = np.flip(tr.data)

        return traces

    def check_filter_parameters(self):
        """
        Checks filter settings based on user parameters and user provided
        filter corners
        """
        assert PAR.FILTER.upper() in ["BANDPASS", "LOWPASS", "HIGHPASS"]

        if PAR.FILTER.upper() == "BANDPASS":
            # Check that filter parameters are provided
            if "MIN_FREQ" not in PAR and "MAX_FREQ" not in PAR:
                raise ParameterError("MIN_FREQ / MAX_FREQ>")
            if "MIN_PERIOD" not in PAR and "MAX_PERIOD" not in PAR:
                raise ParameterError("MIN_PERIOD / MAX_PERIOD")

            # Assign the corresponding frequencies or periods
            if "MIN_FREQ" in PAR:
                PAR.MIN_PERIOD = 1 / PAR.MAX_FREQ
                PAR.MAX_PERIOD = 1 / PAR.MIN_FREQ
            elif "MIN_PERIOD" in PAR:
                PAR.MIN_FREQ = 1 / PAR.MAX_PERIOD
                PAR.MAX_FREQ = 1 / PAR.MIN_PERIOD

            # Check that the values provided make sense
            assert PAR.MIN_FREQ > 0., "Minimum frequency must be > 0"
            assert PAR.MIN_FREQ < PAR.MAX_FREQ, "Max freq < min freq"
            assert PAR.FREQMAX < np.inf, "Max freq > infity"

        elif PAR.FILTER.upper() in ["LOWPASS", "HIGHPASS"]:
            if "FREQ" not in PAR or "PERIOD" not in PAR:
                raise ParameterError("FREQ / PERIOD")

            if PAR.FREQ:
                PAR.PERIOD = 1 / PAR.FREQ
            elif PAR.PERIOD:
                PAR.FREQ = 1 / PAR.PERIOD

            assert PAR.FREQ > 0., "Freq must be > 0"
            assert PAR.FREQ < np.inf, "Freq > infinity"

    def check_mute_parameters(self):
        """
        Checks mute settings, which are used to zero out early or late arrivals
        or offsets
        """
        assert getset(PAR.MUTE) <= {'MuteEarlyArrivals',
                                    'MuteLateArrivals',
                                    'MuteShortOffsets',
                                    'MuteLongOffsets'}

        if 'MuteEarlyArrivals' in PAR.MUTE:
            assert 'MUTE_EARLY_ARRIVALS_SLOPE' in PAR
            assert 'MUTE_EARLY_ARRIVALS_CONST' in PAR
            assert PAR.MUTE_EARLY_ARRIVALS_SLOPE >= 0.

        if 'MuteLateArrivals' in PAR.MUTE:
            assert 'MUTE_LATE_ARRIVALS_SLOPE' in PAR
            assert 'MUTE_LATE_ARRIVALS_CONST' in PAR
            assert PAR.MUTE_LATE_ARRIVALS_SLOPE >= 0.

        if 'MuteShortOffsets' in PAR.MUTE:
            assert 'MUTE_SHORT_OFFSETS_DIST' in PAR
            assert 0 < PAR.MUTE_SHORT_OFFSETS_DIST

        if 'MuteLongOffsets' in PAR.MUTE:
            assert 'MUTE_LONG_OFFSETS_DIST' in PAR
            assert 0 < PAR.MUTE_LONG_OFFSETS_DIST

        if 'MuteShortOffsets' not in PAR.MUTE:
            setattr(PAR, 'MUTE_SHORT_OFFSETS_DIST', 0.)

        if 'MuteLongOffsets' not in PAR.MUTE:
            setattr(PAR, 'MUTE_LONG_OFFSETS_DIST', 0.)

    def check_normalize_parameters(self):
        """
        Check that the normalization parameters are properly set
        """
        assert getset(PAR.NORMALIZE) < {'NormalizeTracesL1',
                                        'NormalizeTracesL2',
                                        'NormalizeEventsL1',
                                        'NormalizeEventsL2'}

    def get_time_scheme(self, traces):
        """
        FIXME: extract time scheme from trace headers rather than parameters

        :param traces:
        :return:
        """
        nt = PAR.NT
        dt = PAR.DT
        t0 = 0.
        return nt, dt, t0

    def get_network_size(self, traces):
        """

        :param traces:
        :return:
        """
        nrec = len(traces)
        nsrc = 1
        return nrec, nsrc

    def get_receiver_coords(self, st):
        """
        Retrieve the coordinates from a Stream object

        :type st: obspy.core.stream.Stream
        :param st: a stream to query for coordinates
        :return:
        """
        if PAR.FORMAT.upper == "SU":
            rx, ry, rz = [], [], []

            for tr in st:
                rx += [tr.stats.su.trace_header.group_coordinate_x]
                ry += [tr.stats.su.trace_header.group_coordinate_y]
                rz += [0.]
            return rx, ry, rz
        else:
            raise NotImplementedError

    def get_source_coords(self, st):
        """
        Get the coordinates of the source object

        :type st: obspy.core.stream.Stream
        :param st: a stream to query for coordinates
        :return:
        """
        if PAR.FORMAT.upper() == "SU":
            sx, sy, sz = [], [], []
            for tr in st:
                sx += [tr.stats.su.trace_header.source_coordinate_x]
                sy += [tr.stats.su.trace_header.source_coordinate_y]
                sz += [0.]
            return sx, sy, sz
        else:
            raise NotImplementedError
    
    def finalize(self):
        """
        Any finalization processes that need to take place at the end of an iter
        """
        pass
