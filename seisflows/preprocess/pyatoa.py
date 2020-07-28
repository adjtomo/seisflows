#!/usr/bin/env python
"""
This is the base class seisflows.preprocess.Pyatoa

This is a main Seisflows class, it controls the preprocessing.
This class uses the Python package Pyatoa to perform preprocessing, and
misfit measurement.
"""
import os
import sys
import copy
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

    @staticmethod
    def check():
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
       
        # Wee bit hacky, but make the 'residuals' directories for saving misfit
        for path_ in [os.path.join(PATH.GRAD, "residuals"),
                      os.path.join(PATH.FUNC, "residuals")]:
            unix.mkdir(path_)

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
        # Late import because preprocess is loaded before optimize, 
        # Optimize required to know which iteration/step_count we are at
        optimize = sys.modules["seisflows_optimize"]

        # Set the logging level, which will be outputted to stdout
        for log in ["pyatoa", "pyflex", "pyadjoint"]:
            logging.getLogger(log).setLevel(PAR.LOGGING.upper())

        # Only query FDSN for i00s00, else turn off by setting client to None
        if optimize.iter == 1 and optimize.line_search.step_count == 0:
            client = PAR.CLIENT
        else:   
            client = None

        # Establish the Pyatoa Configuration object using Seisflows parameters
        config = pyatoa.Config(
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
                        cfgpaths={
                            "waveforms": [os.path.join(PATH.DATA, "mseeds"),
                                          os.path.join(cwd, "traces", "obs")],
                            "synthetics": [os.path.join(cwd, "traces", "syn")],
                            "responses": [os.path.join(PATH.DATA, "seed")],
                            }
                        )

        # Begin processing using Pyatoa
        misfit, nwin = 0, 0
        inv = pyatoa.read_stations(os.path.join(cwd, "DATA", "STATIONS"))

        with ASDFDataSet(os.path.join(self.data, f"{source_name}.h5")) as ds:
            clean_dataset(ds, iteration=optimize.iter, 
                          step_count=optimize.line_search.step_count
                          )

            config.write(write_to=ds)
            mgmt = pyatoa.Manager(ds=ds, config=config)
            for net in inv:
                for sta in net:
                    pyatoa.logger.info(
                        f"\n{'='*80}\n\n{net.code}.{sta.code}\n\n{'='*80}")
                    mgmt.reset()

                    # Gather data; if fail, move onto the next station
                    try:
                        mgmt.gather(station_code=f"{net.code}.{sta.code}.*.HH*")
                    except pyatoa.ManagerError as e:
                        pyatoa.logger.warning(e)
                        continue

                    # Process data; if fail, move onto waveform plotting
                    try:
                        mgmt.flow(fix_windows=PAR.FIX_WINDOWS)

                        self.write_adjoint_traces(
                                       path=os.path.join(cwd, "traces", "adj"),
                                       adjsrcs=mgmt.adjsrcs.values(),
                                       offset=mgmt.stats.time_offset_sec
                                       )

                        misfit += mgmt.stats.misfit
                        nwin += mgmt.stats.nwin
                    except pyatoa.ManagerError as e:
                        pyatoa.logger.warning(e)
                        pass
                    except Exception as e:
                        # Uncontrolled exceptions warrant lengthier error msg
                        pyatoa.logger.warning(e, exc_info=True)
                        pass

                    if PAR.PLOT:
                        # Save the data with a specific tag
                        # e.g. path/to/figures/i01s00_NZ_BFZ.png
                        mgmt.plot(corners=PAR.MAP_CORNERS, show=False, 
                                  save=os.path.join(self.figures, source_name, 
                                                    f"{config.iter_tag}"
                                                    f"{config.step_tag}_"
                                                    f"{net.code}_{sta.code}.png"
                                                    )
                                  )

        # Generate the necessary files to continue the inversion
        if misfit:
            # Write the event misfit a la Tape et al. (2010)
            np.savetxt(os.path.join(path, "residuals", source_name),
                       [0.5 * misfit / nwin], fmt="%11.6e")

            # Create blank adjoint sources and STATIONS_ADJOINT
            self.write_additional_adjoint_files(cwd=cwd, inv=inv)

    def finalize(self):
        """
        Run some serial finalization tasks specific to Pyatoa, which will help
        aggregate the collection of output information:
            Aggregate misfit windows using the Inspector class
            Generate PDFS of waveform figures for easy access
        """
        insp = pyatoa.Inspector(PAR.TITLE)
        insp.discover(path=os.path.join(self.data))
        insp.save(path=self.data) 

        # Combine images into a PDF, delete all pngs afterwards
        unix.cd(self.figures)
        for source in glob("*"):
            if os.path.isdir(source):
                unix.cd(source)
                all_imgs = glob("*.png")
                img_tags = set([_.split("_")[0] for _ in all_imgs])
                for tag in img_tags:
                    fids = glob(f"{tag}_*.png")
                    imgs_to_pdf(fids=fids, 
                                fid_out=os.path.join(self.figures, 
                                                     f"{tag}_{source}.pdf"))
                unix.rm(glob("*.png"))


    @staticmethod
    def sum_residuals(files):
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

    @staticmethod
    def write_adjoint_traces(path, adjsrcs, offset):
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

    @staticmethod
    def write_additional_adjoint_files(cwd, inv):
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
        adj_path = os.path.join(cwd, "traces", "adj")
        adjoint_traces = glob(os.path.join(adj_path, "*"))

        # Open up the stations adjoint file to be written to
        with open(os.path.join(cwd, "DATA", "STATIONS_ADJOINT"), "w") as f:

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
                        for comp in PAR.COMPONENTS:
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
