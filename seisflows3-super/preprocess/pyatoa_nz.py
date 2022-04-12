#!/usr/bin/env python3
"""
This is the sub class seisflows.preprocess.PyatoaMaui

Slightly altered processing function for the New Zealand tomography scenario
"""
import sys
import pyatoa
from pyasdf import ASDFDataSet
from seisflows3.config import custom_import
from pyatoa.utils.read import read_station_codes

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]


def process_event(self, source_name, codes=None, **kwargs):
    """
    .. note::
        MONKEY PATCH: This patch of the Pyaflowa function provides a slight
        alteration to the Pyaflowa source code by treating some data differently
        in the processing workflow
    
    The main processing function for Pyaflowa misfit quantification.

    Processes waveform data for all stations related to a given event,
    produces waveform and map plots during the processing step, saves data
    to an ASDFDataSet and writes adjoint sources and STATIONS_ADJOINT file,
    required by SPECFEM3D's adjoint simulations, to disk.

    Kwargs passed to pyatoa.Manager.flow() function.

    :type source_name: str
    :param source_name: event id to be used for data gathering, processing
    :type codes: list of str
    :param codes: list of station codes to be used for processing. If None,
        will read station codes from the provided STATIONS file
    :rtype: float
    :return: the total scaled misfit collected during the processing chain
    """
    # Create the event specific configurations and attribute container (io)
    io = self.setup(source_name)

    # Allow user to provide a list of codes, else read from station file
    if codes is None:
        codes = read_station_codes(io.paths.stations_file,
                                   loc="*",cha="HH?")

    # Open the dataset as a context manager and process all events in serial
    with ASDFDataSet(io.paths.ds_file) as ds:
        mgmt = pyatoa.Manager(ds=ds, config=io.config)
        for code in codes:
            net, sta, loc, cha = code.split(".")
            # Don't remove response for temp networks as they are already
            # in physical units. Redundant bool for clarity
            rem_resp = bool(net.upper() not in ["ZX", "Z8"])
            mgmt_out, io = self.process_station(mgmt=mgmt, code=code,
                                                io=io, remove_response=rem_resp,
                                                **kwargs)

    scaled_misfit = self.finalize(io)

    return scaled_misfit


class PyatoaNz(custom_import("preprocess", "pyatoa")):
    """
    Data preprocessing class using the Pyatoa package with a custom processing
    function to deal with NZ data
    """
    def prepare_eval_grad(self, source_name, **kwargs):
        """
        Prepare the gradient evaluation by gathering, preprocessing waveforms,
        and measuring misfit between observations and synthetics using Pyatoa.

        This is a process specific task and intended to be run in parallel

        :type path: str
        :param path: path to the current function evaluation for saving residual
        :type source_name: str
        :param source_name: the event id to be used for tagging and data lookup
        """
        # Late import because preprocess is loaded before optimize,
        # Optimize required to know which iteration/step_count we are at
        optimize = sys.modules["seisflows_optimize"]

        # Apply the monkey patch before initiating
        pyatoa.Pyaflowa.process_event = process_event

        # Inititate the Pyaflowa class which abstracts processing functions
        # Communicate to Pyaflowa the current iteration and step count
        pyaflowa = pyatoa.Pyaflowa(structure="seisflows", sfpaths=PATH,
                                   sfpar=PAR, iteration=optimize.iter,
                                   step_count=optimize.line_search.step_count)

        # Process all the stations for a given event using Pyaflowa
        misfit = pyaflowa.process_event(source_name, 
                                        fix_windows=PAR.FIX_WINDOWS,
                                        event_id_prefix=PAR.SOURCE_PREFIX)

        # Generate the necessary files to continue the inversion
        cwd = pyaflowa.path_structure.cwd.format(source_name=source_name) 

        # Event misfit defined by Tape et al. (2010) written to solver dir. 
        self.write_residuals(path=cwd, scaled_misfit=misfit)  

