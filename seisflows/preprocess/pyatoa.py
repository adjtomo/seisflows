#!/usr/bin/env python
"""
This is the base class seisflows.preprocess.Pyatoa

This is a main Seisflows class, it controls the preprocessing.
This class uses the Python package Pyatoa to perform preprocessing, and
misfit measurement.
"""
import os
import sys
import pyatoa
import logging
import numpy as np
from glob import glob
from pyasdf import ASDFDataSet
from seisflows.tools import unix
from seisflows.tools.err import ParameterError
from pyatoa.utils.asdf.clean import clean_dataset
from pyatoa.utils.images import imgs_to_pdf

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]


class Pyatoa:
    """
    Data preprocessing class using the Pyatoa package
    """
    def __init__(self):
        """
        These parameters should not be set by __init__!
        Attributes are just initialized as NoneTypes for clarity and docstrings

        :type data: str
        :param data: directory where data from the preprocessing is stored
        :type figures: str
        :param figures: directory where figures are stored
        :type config: pyatoa.core.Config
        :param config: a general config object that will be parsed into
            the preprocessing workflow
        """
        self.data = None
        self.figures = None
        self.config = None

    def check(self):
        """ 
        Checks Parameter and Path files, will be run at the start of a Seisflows
        workflow to ensure that things are set appropriately.
        """
        # Check the path requirements
        if "PREPROCESS" not in PATH:
            setattr(PATH, "PREPROCESS", 
                    os.path.join(PATH.SCRATCH, "preprocess"))

        if "DATA" not in PATH:
            setattr(PATH, "DATA", None)

        if "RESPONSE" not in PATH:
            setattr(PATH, "RESPONSE", None)

        # Check the existence of required parameters
        required_parameters = ["COMPONENTS", "UNIT_OUTPUT", "MIN_PERIOD",
                               "MAX_PERIOD", "CORNERS", "CLIENT", "ROTATE",
                               "ADJ_SRC_TYPE", "PYFLEX_PRESET",
                               "FIX_WINDOWS", "PLOT", "FORMAT"
                               ]
        for req in required_parameters:
            if req not in PAR:
                raise ParameterError(PAR, req)

        # Check specific parameter requirements
        if PAR.FORMAT != "ascii":
            raise ValueError("Pyatoa preprocess currently only works with "
                             "the 'ascii' format")

        # Set default values parameters for any non-set parameters
        if "PLOT" not in PAR:
            setattr(PAR, "PLOT", True)

        if "LOGGING" not in PAR:
            setattr(PAR, "LOGGING", "DEBUG")

        if "MAP_CORNERS" not in PAR:
            setattr(PAR, "MAP_CORNERS", None)

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

    def setup(self):
        """
        Sets up data preprocessing machinery by establishing an internally
        defined directory structure that will be used to store the outputs 
        of the preprocessing workflow

        Akin to an __init__ class, but to be called externally by the workflow.
        """
        # Late import because preprocess is loaded before optimize
        solver = sys.modules["seisflows_solver"]

        self.data = os.path.join(PATH.PREPROCESS, "data")
        self.figures = os.path.join(PATH.PREPROCESS, "figures")

        # Make data and figure directories for each source
        unix.mkdir(self.data)  
        for source_name in solver.source_names:
            unix.mkdir(os.path.join(self.figures, source_name))

    def prepare_eval_grad(self, path, cwd, source_name):
        """
        Prepare the gradient evaluation by gathering, preprocessing waveforms, 
        and measuring misfit between observations and synthetics using Pyatoa.
        
        This is a process specific task and intended to be run in parallel

        :type path: str
        :param path: path to the current function evaluation for saving residual
        :type cwd: str
        :param cwd: the path to the current Specfem working directory
        :type source_name: str
        :param source_name: the event id to be used for tagging and data lookup
        """

        # Set the logging level, which will be outputted to stdout
        for log in ["pyatoa", "pyflex", "pyadjoint"]:
            logging.getLogger(log).setLevel(PAR.LOGGING.upper())

        # Generate auxiliary information necessary for misfit assessment
        config = self.create_config(cwd=cwd, source_name=source_name)
        inv = pyatoa.read_stations(os.path.join(cwd, "DATA", "STATIONS"))

        # Misfit assessment for a given event and set of stations
        misfit, nwin = self.process_event(config, inv, cwd)

        # Generate the necessary files to continue the inversion
        if misfit:
            # Event misfit defined by Tape et al. (2010)
            self.write_residuals(path=path, scaled_misfit=0.5 * misfit / nwin,
                                 source_name=source_name)

            # Create blank adjoint sources and STATIONS_ADJOINT
            self.write_additional_adjoint_files(cwd=cwd, inv=inv)

    def finalize(self):
        """
        Run some serial finalization tasks specific to Pyatoa, which will help
        aggregate the collection of output information:
            Aggregate misfit windows using the Inspector class
            Generate PDFS of waveform figures for easy access
        """
        insp = pyatoa.Inspector(PAR.TITLE, verbose=False)
        insp.discover(path=os.path.join(self.data))
        insp.save(path=self.data) 

        self.make_pdfs()

    def write_residuals(self, path, scaled_misfit, source_name):
        """
        Computes residuals and saves them to a text file in the appropriate path

        :type path: str        
        :param path: scratch directory path, e.g. PATH.GRAD or PATH.FUNC
        :type scaled_misfit: float
        :param scaled_misfit: the summation of misfit from each 
            source-receiver pair calculated by prepare_eval_grad()
        :type source_name: str
        :param source_name: name of the source related to the misfit, used
            for file naming
        """
        residuals_dir = os.path.join(path, "residuals")        

        if not os.path.exists(residuals_dir):
            unix.mkdir(residuals_dir)
        
        event_residual = os.path.join(residuals_dir, source_name)        
     
        np.savetxt(event_residual, [scaled_misfit], fmt="%11.6e")

    def sum_residuals(self, files):
        """
        Averages the event misfits and returns the total misfit.

        Total misfit defined by Tape et al. (2010)

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

    def check_fixed_windows(self, config):
        """
        Determine how to address fixed time windows based on the user parameter
        as well as the current iteration and step count in relation to the 
        inversion location.
        
        Options:
            True: Always fix windows except for i01s00 because we don't have any
                  windows for the first function evaluation
            False: Don't fix windows, always choose a new set of windows
            Iter: Pick windows only on the initial step count (0th) for each 
                  iteration. WARNING - does not work well with Thrifty Inversion 
                  because the 0th step count is usually skipped
            Once: Pick new windows on the first function evaluation and then fix
                  windows. Useful for when parameters have changed, e.g. filter
                  bounds

        :type config: pyatoa.core.config.Config
        :param config: Config object that should contain current iteration
            and step_count values
        :rtype: bool
        :return: flag to denote whether or not to pick new windows in the 
            processing workflow
        """
        # First function evaluation never fixes windows 
        if config.iteration == 1 and config.step_count == 0:
            fix_windows = False
        # By 'iter'ation only pick new windows on the first step count
        elif PAR.FIX_WINDOWS.upper() == "ITER":
            if config.step_count == 0:
                fix_windows = False
            else:
                fix_windows = True
        # 'Once' picks windows only for the first function evaluation of the 
        # current set of iterations.
        elif PAR.FIX_WINDOWS.upper() == "ONCE":
            if config.iteration == PAR.BEGIN and config.step_count == 0:
                fix_windows = False 
            else:
                fix_windows = True
        # Bool fix windows simply sets the parameter
        else:
            fix_windows = PAR.FIX_WINDOWS

        return fix_windows

    def create_config(self, cwd, source_name):
        """
        Creates the Pyatoa Configuration object using the internal Seisflows
        parameters and some unique identifiers from the preprocess module.

        :type cwd: str
        :param cwd: the path to the current Specfem working directory
        :type source_name: str
        :param source_name: the event id to be used for tagging and data lookup
        """
        # Late import because preprocess is loaded before optimize,
        # Optimize required to know which iteration/step_count we are at
        optimize = sys.modules["seisflows_optimize"]

        # Only query FDSN for i00s00, else turn off by setting client to None
        # Dont fix windows for the first function evaluation
        if optimize.iter == 1 and optimize.line_search.step_count == 0:
            client = PAR.CLIENT
        else:
            client = None

        # Establish the Pyatoa Configuration object using Seisflows parameters
        return pyatoa.Config(
                        event_id=source_name,
                        iteration=optimize.iter,
                        step_count=optimize.line_search.step_count,
                        synthetics_only=bool(PAR.CASE.lower() == "synthetic"),
                        component_list=list(PAR.COMPONENTS),
                        rotate_to_rtz=PAR.ROTATE,
                        min_period=PAR.MIN_PERIOD,
                        max_period=PAR.MAX_PERIOD,
                        filter_corners=PAR.CORNERS,
                        unit_output=PAR.UNIT_OUTPUT,
                        client=client,
                        start_pad=PAR.START_PAD,
                        end_pad=PAR.END_PAD,
                        adj_src_type=PAR.ADJ_SRC_TYPE,
                        pyflex_preset=PAR.PYFLEX_PRESET,
                        paths={
                            "waveforms": [os.path.join(PATH.DATA, "mseeds"),
                                          os.path.join(cwd, "traces", "obs")],
                            "synthetics": [os.path.join(cwd, "traces", "syn")],
                            "responses": [os.path.join(PATH.DATA, "seed")],
                            }
                        )

    def process_event(self, config, inv, cwd):
        """
        Run through the Pyatoa internal processing workflow to gather, process
        data and output misfit windows and adjoint sources

        :type config: pyatoa.core.config.Config
        :param config: Configuration object for the given event and function
            evaluation. Created by create_config()
        :type inv: obspy.core.inventory.Inventory
        :param inv: an inventory containing station names to process for a given
            source. Only network and station codes are required in the inventory
        :type cwd: str
        :param cwd: the current SPECFEM solver directory, used only to determine
            where to write adjoint traces.
        :return:
        """
        # Track the total misfit and number of windows
        misfit, nwin = 0, 0
        # Track the number of successes and fails for a summary statement
        _stations, _processed, _exceptions = 0, 0, 0
    

        # Begin the Pyatoa processing workflow
        with ASDFDataSet(os.path.join(self.data,
                                      f"{config.event_id}.h5")) as ds:
            clean_dataset(ds, iteration=config.iteration,
                          step_count=config.step_count)

            config.write(write_to=ds)
            mgmt = pyatoa.Manager(ds=ds, config=config)
            for net in inv:
                for sta in net:
                    _stations += 1
                    pyatoa.logger.info(
                        f"\n{'=' * 80}\n\n{net.code}.{sta.code}\n\n{'=' * 80}")
                    mgmt.reset()

                    # Gather data; if fail, move onto the next station
                    try:
                        mgmt.gather(code=f"{net.code}.{sta.code}.*.HH*")
                    except pyatoa.ManagerError as e:
                        pyatoa.logger.warning(e)
                        continue

                    # Process data; if fail, move onto waveform plotting
                    try:
                        mgmt.flow(fix_windows=self.check_fixed_windows(config))

                        self.write_adjoint_traces(
                            path=os.path.join(cwd, "traces", "adj"),
                            adjsrcs=mgmt.adjsrcs.values(),
                            offset=mgmt.stats.time_offset_sec
                        )

                        misfit += mgmt.stats.misfit
                        nwin += mgmt.stats.nwin
                        _processed += 1
                    except pyatoa.ManagerError as e:
                        pyatoa.logger.warning(e)
                        pass
                    except Exception as e:
                        # Uncontrolled exceptions warrant lengthier error msg
                        pyatoa.logger.warning(e, exc_info=True)
                        _exceptions += 1
                        pass

                    if PAR.PLOT:
                        # Save the data with a specific tag
                        # e.g. path/to/figures/i01s00_NZ_BFZ.png
                        mgmt.plot(corners=PAR.MAP_CORNERS, show=False,
                                  save=os.path.join(self.figures,
                                                    config.event_id,
                                                    f"{config.iter_tag}"
                                                    f"{config.step_tag}_"
                                                    f"{net.code}_{sta.code}.png"
                                                    )
                                  )
        # Record summary information at the end of the Pyatoa log file
        pyatoa.logger.info(f"\n{'=' * 80}\n\nSUMMARY\n\n{'=' * 80}\n"
                           f"SOURCE NAME: {config.event_id}\n"
                           f"STATIONS: {_processed} / {_stations}\n"
                           f"WINDOWS: {nwin}\n"
                           f"RAW MISFIT: {misfit:.2f}\n"
                           f"UNEXPECTED ERRORS: {_exceptions}"
                           )

        return misfit, nwin

    def write_adjoint_traces(self, path, adjsrcs, offset):
        """
        Writes adjoint sources required for gradient computation using the 
        functionality contained in the Pyadjint AdjointSource object

        :type path: str
        :param path: location "adjoint traces" will be written
        :type adjsrcs: dict of pyadjoint.AdjointSource's
        :param adjsrcs: adjoint source objects that contain a write function
        :type offset: float
        :param offset: required time offset that is set by Specfem and defined
            by Pyatoa
        """
        for adj in adjsrcs:
            fid = f"{adj.network}.{adj.station}.{adj.component}.adj"
            adj.write(filename=os.path.join(path, fid), format="SPECFEM",
                      time_offset=offset
                      )

    def write_additional_adjoint_files(self, cwd, inv):
        """
        Generates the STATIONS_ADJOINT file expected by the SPECFEM 
        adjoint simulation, and blank adjoint sources.

        :type cwd: str
        :param cwd: path to the current SPECFEM working directory
        :type inv: obspy.core.inventory.Inventory
        :param inv: Inventory created from the SPECFEM STATIONS file, to be 
            used for checking station names
        """
        # A line template for the STATIONS_ADJOINT file
        tmplt = "{sta:>6}{net:>6}{lat:12.4f}{lon:12.4f}{elv:11.1f}{bur:11.1f}\n"

        # Check that adjoint sources have been written by prepare_eval_grad()
        unix.cd(os.path.join(cwd, "traces", "adj"))
        adjoint_traces = glob("*")

        # Open up the stations adjoint file to be written to
        with open(os.path.join(cwd, "DATA", "STATIONS_ADJOINT"), "w") as f:

            # If no adjoint traces were written, will create empty file
            if not adjoint_traces:
                return

            # Create a zeroed adjoint source trace for filling blanks
            example_trace = np.loadtxt(adjoint_traces[0])
            example_trace[:, 1] = 0

            # Check for the existence of adjoint sources for each station
            for net in inv:
                for sta in net:
                    fid = f"{net.code}.{sta.code}.{'{}'}.adj"
                    # Check for the existence of any adjoint sources
                    check_station = glob(fid.format("???"))
                    if check_station:
                        _, _, channel, _ = check_station[0].split(".")
                        for comp in PAR.COMPONENTS:
                            # Take channel name from the file name
                            channel = channel[:-1] + comp
                            fid_out = fid.format(channel.upper())
                            if not os.path.exists(fid_out):
                                # Write a blank adjoint source for empty comps
                                np.savetxt(fid_out, example_trace)

                        # Write the station into the STATIONS_ADJOINT file
                        f.write(tmplt.format(sta=sta.code, net=net.code,
                                             lat=sta.latitude, 
                                             lon=sta.longitude, 
                                             elv=sta.elevation,
                                             bur=0)
                                )

    def make_pdfs(self):
        """
        Utility function to make PDFs of waveform images for easier transport
        and access during an inversion. Checks tags based on the file names
        specified during process_event(). Removes png files after pdf made.
        """
        unix.cd(self.figures)
        sources = [os.path.abspath(_) for _ in glob("*") if os.path.isdir(_)]
        for source in sources:
            event = os.path.basename(source)
            unix.cd(source)
            all_imgs = glob("*.png")
            img_tags = set([_.split("_")[0] for _ in all_imgs])
            for tag in img_tags:
                fids = glob(f"{tag}_*.png")
                imgs_to_pdf(fids=fids,
                            fid_out=os.path.join(self.figures,
                                                 f"{tag}_{event}.pdf"))
            unix.rm(glob("*.png"))
