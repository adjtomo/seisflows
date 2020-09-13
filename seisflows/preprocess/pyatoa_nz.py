#!/usr/bin/env python
"""
This is the sub class seisflows.preprocess.PyatoaMaui

Slightly altered processing function for the New Zealand tomography scenario
"""
import os
import sys
import pyatoa
from glob import glob
from pyasdf import ASDFDataSet
from seisflows.tools import unix
from seisflows.config import custom_import
from pyatoa.utils.images import merge_pdfs
from pyatoa.utils.asdf.clean import clean_dataset

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]


class PyatoaNz(custom_import("preprocess", "pyatoa")):
    """
    Data preprocessing class using the Pyatoa package with a custom processing
    function to deal with NZ data
    """
    def process_event(self, config, inv, cwd):
        """
        Run through the Pyatoa internal processing workflow to gather, process
        data and output misfit windows and adjoint sources.

        Ignores instrument response removal for certain network codes

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
                # ADDITION:
                # Bannister network data has already had the instrument
                # response removed, so we will ignore that.
                if net.code in ["Z8", "ZX"]:
                    remove_response = False
                else:
                    remove_response = True

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

                    # ADDITION:
                    # Occasionally multiple traces are available for a single
                    # component which is not ideal. Ignore these, they don't 
                    # occur often enough to worry about
                    if mgmt.stats.len_obs > 3:
                        pyatoa.logger.info("Too many observations, ignoring...")
                        continue

                    # Process data; if fail, move onto waveform plotting
                    try:
                        mgmt.flow(fix_windows=self.check_fixed_windows(config),
                                  remove_response=remove_response)

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

                    # Save the data with a specific tag
                    # e.g. path/to/figures/i01s00_NZ_BFZ.png
                    if PAR.PLOT:
                        # Save the data with a specific tag
                        # e.g. path/to/figures/i01s00_NZ_BFZ.pdf
                        mgmt.plot(corners=PAR.MAP_CORNERS, show=False,
                                  save=os.path.join(self.figures,
                                                    config.event_id,
                                                    f"{config.iter_tag}"
                                                    f"{config.step_tag}_"
                                                    f"{net.code}_{sta.code}.pdf"
                                                    )
                                  )

        # To reduce file count, merge all output pdfs into a single document
        all_created_pdfs = glob(os.path.join(self.figures, config.event_id,
                                             "*.pdf"))
        merge_pdfs(fids=sorted(all_created_pdfs),
                   fid_out=os.path.join(self.figures, f"{config.iter_tag}"
                                                      f"{config.step_tag}_"
                                                      f"{config.event_id}.pdf")
                   )
        unix.rm(all_created_pdfs)

        # Record summary information at the end of the Pyatoa log file
        pyatoa.logger.info(f"\n{'=' * 80}\n\nSUMMARY\n\n{'=' * 80}\n"
                           f"SOURCE NAME: {config.event_id}\n"
                           f"STATIONS: {_processed} / {_stations}\n"
                           f"WINDOWS: {nwin}\n"
                           f"RAW MISFIT: {misfit:.2f}\n"
                           f"UNEXPECTED ERRORS: {_exceptions}"
                           )

        return misfit, nwin
