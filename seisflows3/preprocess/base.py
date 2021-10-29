#!/usr/bin/env python
"""
This is the main seisflows.preprocess.base

This is a main Seisflows class, it controls the preprocessing of waveforms and
generation of adjoint sources
"""
from seisflows.config import SeisFlowsPathsParameters


class Base:
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

        return sf

    def check(self, validate=True):
        """ 
        Checks parameters and paths
        """
        pass

    def setup(self):
        """
        Sets up data preprocessing machinery
        """
        pass

    def prepare_eval_grad(self):
        """
        Prepares solver for gradient evaluation by quantifying misfit and 
        writing adjoint sources
        """
        pass

    def sum_residuals(self, files):
        """
        Sums residuals from each source to generate the total misfit
        """
        pass

    def finalize(self):
        """
        Any finalization processes that need to take place at the end of an iter
        """
        pass
