#!/usr/bin/env python
"""
This is the subclass seisflows.preprocess.pyatoa

This is a main Seisflows class, it controls the preprocessing.
This subclass uses the Python package Pyatoa to perform preprocessing, and
misfit measurement.
"""
import os
import sys
import time
import obspy
import numpy as np
from glob import glob

from seisflows.config import custom_import
from seisflows.tools import signal
from seisflows.tools.err import ParameterError
from seisflows.tools.tools import exists, getset

from pyatoa import Pyaflowa

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class Pyatoa(custom_import("preprocess", "base")):
    """
    Data preprocessing class

    Outsources data handling to Pyatoa via the class Pyaflowa. All calls are
    made external, this class is simply used as a Seisflows abstraction for
    calls to Pyatoa.
    """
    def __init__(self):
        """
        Initate empty variables
        """
        self.pyaflowa = None

    def check(self):
        """ 
        Checks parameters and paths
        """
        pass

    def setup(self):
        """
        Sets up data preprocessing machinery by initiating the Pyaflowa class,
        which will parse directory structure and generate the relevant
        directories that are to be filled during the workflow
        """
        self.pyaflowa = Pyaflowa(par=vars(PAR), paths=vars(PATH))

    def prepare_eval_grad(self, path='.'):
        """
        Prepares solver for gradient evaluation by writing residuals and
        adjoint traces.

        This functionality is already taken care of by Pyaflowa.process()

        :type path: str
        :param path: directory containing observed and synthetic seismic data
        """
        pass

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
            residuals.extend(list(np.loadtxt(filename)))

        np.savetxt(filename, residuals)

    def sum_residuals(self, files):
        """
        Sums the average misfit from each event, average by number of events

        :type files: str
        :param files: list of single-column text files containing residuals
        :rtype: float
        :return: average misfit
        """
        assert(len(files) == PAR.NTASK), \
            "Number of misfit files does not match the number of events"
        total_misfit = 0
        for filename in files:
            total_misfit += np.sum(np.loadtxt(filename))

        # Save the total misfit
        total_misfit = total_misfit / PAR.NTASK

        return total_misfit

    def write_adjoint_traces(self, path, syn, obs, channel):
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

        self.writer(adj, path, channel)

