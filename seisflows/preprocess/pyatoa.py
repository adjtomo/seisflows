#!/usr/bin/env python
"""
This is the base class seisflows.preprocess.pyatoa

This is a main Seisflows class, it controls the preprocessing.
This class uses the Python package Pyatoa to perform preprocessing, and
misfit measurement.
"""
import os
import sys
import time
import copy
import pyasdf
import pyatoa
import logger
import numpy as np
from glob import glob

from pyatoa.utils.asdf.clean import clean_dataset
from pyatoa.utils.form import format_iter, format_step

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]
solver = sys.modules['seisflows_solver']
optimize = sys.modules['seisflows_optimize']


class Pyatoa(object):
    """
    Data preprocessing class

    Outsources data handling to Pyatoa via the class Pyaflowa. All calls are
    made external, this class is simply used as a Seisflows abstraction for
    calls to Pyatoa.
    """
    def check(self):
        """ 
        Checks Parameter and Path files, will be run at the start of a Seisflows
        workflow to ensure that things are set appropriately.
        """
        # Check the path requirements
        if "PYATOA" not in PATH:
            setattr(PATH, "PYATOA", os.path.join(PATH.SCRATCH, "pyatoa"))

        if "DATA" not in PATH:
            setattr(PATH, "DATA", [])

        # Check the parameter requirements
        if PAR.FORMAT != "ascii":
            raise ValueError("Pyatoa preprocess currently only works with "
                             "the 'ascii' format")

        if "UNIT_OUTPUT" not in PAR:
            raise ValueError("Pyatoa preprocess requires parameter "
                             "UNIT_OUTPUT to define the desired waveform "
                             "units")

        if "MIN_PERIOD" not in PAR:
            raise ValueError("Pyatoa preprocess requires parameter "
                             "MIN_PERIOD to define the lower filter corner")

        if "MAX_PERIOD" not in PAR:
            raise ValueError("Pyatoa preprocess requires parameter "
                             "MAX_PERIOD to define the upper filter corner")

        if "CORNERS" not in PAR:
            raise ValueError("Pyatoa preprocess requires parameter "
                             "CORNERS to define the filter corner order")

        if "ROTATE" not in PAR:
            raise ValueError("Pyatoa preprocess requires parameter "
                             "ROTATE to define the desired components")

        if "CLIENT" not in PAR:
            setattr(PAR, "CLIENT", None)

        # Used to define the start time of fetched observation waveforms
        if "START_PAD" not in PAR:
            setattr(PAR, "START_PAD", 20)

        # Used to define the end time of fetched observation waveforms
        if "END_PAD" not in PAR:
            setattr(PAR, "END_PAD", PAR.DT * PAR.NT + PAR.START_PAD + 5)
        else:
            if PAR.DT * PAR.NT >= PAR.START_PAD + PAR.END_PAD:
                raise ValueError("Pyatoa preprocess parameters START_PAD and "
                                 "END_PAD will not provide long enough obs."
                                 "traces to match the length of synthetics")

        if "ADJ_SRC_TYPE" not in PAR:
            raise ValueError("Pyatoa preprocess requires parameter "
                             "ADJ_SRC_TYPE to define the desired adjoint "
                             "source type")

        if "PYFLEX_PRESET" not in PAR:
            raise ValueError("Pyatoa preprocess requires parameter "
                             "PYFLEX_PRESET to define the desired map "
                             "for Pyflex parameters")

        if "FIX_WINDOWS" not in PAR:
            raise ValueError("Pyatoa preprocess requires parameter "
                             "FIX_WINDOWS to define windowing selection "
                             "behavior")

        if "PLOT" not in PAR:
            setattr(PAR, "PLOT", True)

        if "LOGGING" not in PAR:
            setattr(PAR, "LOGGING", "DEBUG")

        if "MAP_CORNERS" not in PAR:
            setattr(PAR, "MAP_CORNERS", None)

    def setup(self):
        """
        Sets up data preprocessing machinery by establishing an internally
        defined directory structure that will be used to store the outputs 
        of the preprocessing workflow

        Akin to an __init__ class, but to be called externally by the workflow.
        """
        self.data = os.path.join(PATH.PYATOA, "data")
        self.figures = os.path.join(PATH.PYATOA, "figures")

        # Set the logging level, to be outputted to stdout
        for log in ["pyflex", "pyflex", "pyadjoint"]:
            logging.getLogger(log).setLevel(PAR.LOGGING)

        # Establish the Pyatoa Configuration object using Seisflows parameters
        self.config = pyatoa.Config(
            event_id=None, iteration=None, step_count=None, 
            synthetics_only=bool(PAR.CASE.lower() == "synthetic"),
            min_period=PAR.MIN_PERIOD, max_period=PAR.MAX_PERIOD, 
            filter_corners=PAR.CORNERS, unit_output=PAR.UNIT_OUTPUT,
            rotate_to_rtz=PAR.ROTATE, client=PAR.CLIENT,
            start_pad=PAR.START_PAD, end_pad=PAR.END_PAD,
            adj_src_type=PAR.ADJ_SRC_TYPE, pyflex_preset=PAR.PYFLEX_PRESET,
            cfgpaths={"waveforms": [os.path.join(PATH.DATA, "mseeds")],
                      "responses": [os.path.join(PATH.DATA, "seed")]
                      }
            )

    @property
    def fix_windows(self):
        """
        A property to check if windows need to be fixed. Dependent on the 
        iteration/step count, as well as the User set parameter
        """
        if optimize.iter == 1 and optimize.line_search.step_count == 0:
            return False
        elif isinstance(PAR.FIX_WINDOWS, bool):
            return PAR.FIX_WINDOWS

    def prepare_eval_grad(self, path, source_name):
        """
        Prepare the gradient evaluation by gathering, preprocessing waveforms, 
        and measuring misfit between observations and synthetics using Pyatoa.

        :type path: str
        :param path: path to the current Specfem working directory
        :type source_name: str
        :param source_name: the event id to be used for tagging and data lookup
        """
        # Some internal path naming and parameter setting
        dataset = os.path.join(self.data, source_name)
        figures = os.path.join(self.figures, source_name)
        misfit = os.path.join(self.data, "misfit")

        # Establish event-specific configuration parameters
        config = copy.deepcopy(self.config)
        config.source_name = source_name
        config.iteration = optimize.iter
        config.step_count = optimize.line_search.step_count
        config.cfgpaths["waveforms"] += [os.path.join(path, "traces", "obs")]
        config.cfgpaths["synthetics"] += [os.path.join(path, "traces", "syn")]

        # Begin processing using Pyatoa
        misfit, nwin = 0, 0
        inv = pyatoa.read_stations(os.path.join(path, "DATA", "STATIONS"))

        with pyasdf.ASDFDataSet(dataset):
            clean_dataset(dataset, iteration=optimize.iter, 
                          step_count=optimize.line_search.step_count)

            config.write(write_to=dataset)
            mgmt = pytoa.Manager(ds=ds, config=config)
            for network in inventory:
                for station in network:
                    mgmt.reset()

                    # Gather data; if fail, move onto the next station
                    try:
                        mgmt.gather(station_code=f"{net.code}.{sta.code}.*.HH*")
                    except pyatoa.ManagerError as e:
                        pyatoa.logger.warning(e)
                        continue

                    # Process data; if fail, move onto waveform plotting
                    try:
                        mgmt.flow(preprocess_overwrite=overwrite,
                                  fix_windows=self.fix_windows)

                        self.write_adjoint_traces(
                                       path=os.path.join(path, "traces", "adj"),
                                       traces=mgmt.adjsrcs.values(), 
                                       offset=mgmt.stats.time_offset_sec
                                       )

                        misfit += mgmt.misfit
                        nwin += mgmt.nwin
                    except Exception as e:
                        pyatoa.logger.warning(e, exc_info=True)
                        pass

                    # Plot the data
                    if PAR.PLOT:
                        mgmt.plot(save=os.path.join(figures, f"{sta.code}.png"),
                                  map_corners=PAR.MAP_CORNERS, show=False, 
                                  return_figure=False)

        # Generate the necessary files to continue the inversion
        if misfit:
            # Write the event misfit
            np.savetxt(os.path.join(misfit, source_name), 
                       0.5 * misfit / num_win)

            # Create blank adjoint sources and STATIONS_ADJOINT
            self.write_additional_adjoint_files(
                        path=path, inv=inv, comp_list=mgmt.config.component_list
                        )

    def sum_residuals(self, files):
        """
        Averages the event misfits and returns the total misfit.

        Total misfit defined by Tape et al. (2007): 
            total_misfit = sum(event_misfits) / number_srcs

        :type files: str
        :param files: list of single-column text files containing residuals
            that will have been generated using prepare_eval_grad()
        :rtype: float
        :return: average misfit
        """
        assert(len(files) == PAR.NTASK), \
            "Number of misfit files does not match the number of events"

        total_misfit = 0
        for filename in files:
            total_misfit += np.sum(np.loadtxt(filename))

        total_misfit /= PAR.NTASK

        return total_misfit

    def write_adjoint_traces(self, path, adjsrcs, offset):
        """
        Writes "adjoint traces" required for gradient computation using the 
        functionality contained in the Pyadjint AdjointSource object

        :type path: str
        :param path: location "adjoint traces" will be written
        :type syn: obspy.core.stream.Stream
        :param syn: synthetic data
        :type obs: obspy.core.stream.Stream
        :param syn: observed data
        :type channel: str
        :param channel: channel or component code used by writer
        """
        for adj in adjsrcs:
            fid = f"{adj.network}.{adj.station}.{adj.component}.adj"
            adj.write(filename=os.path.join(path, fid), format="SPECFEM",
                      time_offset=offset
                      )

    def write_additional_adjoint_files(self, path, inv, component_list):
        """
        Generates the STATIONS_ADJOINT file expected by the SPECFEM 
        adjoint simulation, and blank adjoint sources.

        :type path: str
        :param path: path to the current SPECFEM working directory
        :type inv: obspy.core.inventory.Inventory
        :param inv: Inventory created from the SPECFEM STATIONS file, to be 
            used for checking station names
        :type component_list: list of str
        :param component_list: list of components that are being used, for 
            checking against adjoint sources
        """
        # A line template for the STATIONS_ADJOINT file
        tmplt = "{sta:>6}{net:>6}{lat:12.4f}{lon:12.4f}{elv:11.1f}{bur:11.1f}\n"

        # Check that adjoint sources have been written by prepare_eval_grad()
        adj_path = os.path.join(path, "traces", "adj")
        adjoint_traces = glob(os.path.join(adj_path, "*"))

        # Open up the stations adjoint file to be written to
        with open(os.path.join(path, "DATA", "STATIONS_ADJOINT"), "w") as f:

            # If no adjoint traces were written, will create empty file
            if not adjoint_traces:
                return

            # Create an zeroed adjoint source trace for filling blanks
            example_trace = np.loadtxt(adjoint_traces[0])
            example_trace[:, 1] = 0

            # Check for the existence of adjoint sources for each station
            for net in inv:
                for sta in net:
                    fid = os.path.join(adj_path, 
                                       f"{net.code}.{sta.code}.??{'{}'}.adj"
                                       )
                    # Check for the existence of any adjoint sources
                    if glob(fid.format("?")):
                        for comp in component_list:
                            if not os.path.exists(fid.format(comp)):
                                # Write a blank adjoint source for empty comps
                                np.savetxt(fid.format(comp), example_trace)

                        # Write the station into the STATIONS_ADJOINT file
                        f.write(tmplt.format(sta=sta.code, net=net.code,
                                             lat=sta.latitude, 
                                             lon=sta.longitude, 
                                             elv=sta.elevation,
                                             bur=0)
                        )

    def finalize(self):
        """
        Run some finalization functions that will assist in assessing misfit
        behavior 
        """

