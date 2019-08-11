#
# This is Seisflows
#
# See LICENCE file
#
# Functions to write signals to files  (using Obspy)
#
# SeisFlows uses obspy stream objects for holding and processing seismic data.
# In some cases, obspy.read doesn't provide the desired behavior, so we
# introduce an additonal level of indirection
#
# used by the PREPROCESS class and specified by the WRITER parameter
#
###############################################################################

# Import system module
import sys

# Import numpy
import numpy as np

PAR = sys.modules['seisflows_parameters']


def su(stream, path, filename):
    """ Write Seismic Unix files.
        Function writeBigSuFile is a hack to write a .su file when the number
        of samples per trace is two big.
        In the su format only 2 bytes per trace are dedicated to encoding for
        the number of samples (as signed int, see:
        http://lists.swapbytes.de/archives/obspy-users/2017-March/002359.html).
        Even if it's an old format it's still extremely stupid.
        This proove the lack of vision the designer of this format had at that
        time. They could have chosen 8 bytes or 16 bytes it was no big deal...
        They've cost me a day's work.
        But let us forget about the past. This limits the size of the
        traces in the header to maximum 32768.
        We use Obspy to write the file with dummy values there instead of the
        real number of sample (that we now anyway : it is PAR.NT).
        We thus rewrote a quick version of this function from Obspy replacing
        the number of point by PAR.NT

    """
    for t in stream:
        # work around obspy data type conversion
        t.data = t.data.astype(np.float32)

    max_npts = 32767
    max_delta = 0.065535
    dummy_delta = max_delta

    if stream[0].stats.delta > max_delta:
        for t in stream:
            t.stats.delta = dummy_delta

    # write data to file
    if PAR.NT < max_npts:
        stream.write(path+'/'+filename, format='SU')
    else:
        writeBigSuFile(stream, path+'/'+filename)


def ascii(stream, path, filenames):
    """ Write ascii signal file
    """
    for ir, tr in enumerate(stream):
        nt = tr.stats.npts
        t1 = float(tr.stats.starttime)
        t2 = t1 + tr.stats.npts*tr.stats.sampling_rate
        print nt, t1, t2

        t = np.linspace(t1, t2, nt)
        w = tr.data

        print path + '/' + tr.stats.filename
        print times.shape, tr.data.shape
        np.savetxt(path + '/' + tr.stats.filename,
                   np.column_stack((t, w)))


def writeBigSuFile(stream, path, byteorder='<'):
    """ This function is a hack to write a .su file when the number
        of samples per trace is two big.
        In the su format only 2 bytes per trace are dedicated to encoding for
        the number of samples (as signed int, see:
        http://lists.swapbytes.de/archives/obspy-users/2017-March/002359.html).
        Even if it's an old format it's still extremely stupid.
        This proove the lack of vision the designer of this format had at that
        time. They could have chosen 8 bytes or 16 bytes it was no big deal...
        They've cost me a day's work.
        But let us forget about the past. This limits the size of the
        traces in the header to maximum 32768.
        We use Obspy to write the file with dummy values there instead of the
        real number of sample (that we now anyway : it is PAR.NT).
        We thus rewrote a quick version of this function from Obspy replacing
        the number of point by PAR.NT
        This is mostly copy-pastes from Obspy source code
    """

    from obspy.core.utcdatetime import UTCDateTime
    from obspy.core.util import AttribDict
    from obspy.io.segy.core import SUFile
    from obspy.io.segy.segy import SEGYWritingError
    from obspy.io.segy.segy import SEGYTrace
    from obspy.io.segy.header import TRACE_HEADER_FORMAT
    from obspy.io.segy.header import DATA_SAMPLE_FORMAT_PACK_FUNCTIONS

    dummy_npts = 9999
    su_file = SUFile()

    # Add all traces
    for trace in stream:
        new_trace = SEGYTrace()
        new_trace.data = trace.data
        # Use header saved in stats if one exists.
        if hasattr(trace.stats, 'su') and \
           hasattr(trace.stats.su, 'trace_header'):
            this_trace_header = trace.stats.su.trace_header
        else:
            this_trace_header = AttribDict()
        new_trace_header = new_trace.header
        # Again loop over all field of the trace header and if they exists, set
        # them. Ignore all additional attributes.
        for _, item, _, _ in TRACE_HEADER_FORMAT:
            if hasattr(this_trace_header, item):
                setattr(new_trace_header, item,
                        getattr(this_trace_header, item))
        starttime = trace.stats.starttime
        # Set some special attributes, e.g. the sample count and other stuff.
        new_trace_header.number_of_samples_in_this_trace = trace.stats.npts
        new_trace_header.sample_interval_in_ms_for_this_trace = \
            int(round((trace.stats.delta * 1E6)))
        # Set the date of the Trace if it is not UTCDateTime(0).
        if starttime == UTCDateTime(0):
            new_trace.header.year_data_recorded = 0
            new_trace.header.day_of_year = 0
            new_trace.header.hour_of_day = 0
            new_trace.header.minute_of_hour = 0
            new_trace.header.second_of_minute = 0
        else:
            new_trace.header.year_data_recorded = starttime.year
            new_trace.header.day_of_year = starttime.julday
            new_trace.header.hour_of_day = starttime.hour
            new_trace.header.minute_of_hour = starttime.minute
            new_trace.header.second_of_minute = starttime.second
        # Set the data encoding and the endianness.
        new_trace.endian = byteorder
        # Add the trace to the SEGYFile object.
        su_file.traces.append(new_trace)

    # Write the file
    file = open(path, 'wb')
    for trace in su_file.traces:
        trace.header.number_of_samples_in_this_trace = dummy_npts
        endian = trace.endian
        data_encoding = 5
        # Write the header.
        trace.header.write(file, endian=endian)
        # Write the data.
        if trace.data is None:
            msg = "No data in the SEGYTrace."
            raise SEGYWritingError(msg)
        DATA_SAMPLE_FORMAT_PACK_FUNCTIONS[data_encoding](file, trace.data,
                                                         endian=endian)
    file.close()
