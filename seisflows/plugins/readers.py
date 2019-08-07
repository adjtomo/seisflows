#
# This is Seisflows
#
# See LICENCE file
#
# Functions to write signals to files (using Obspy)
#
# SeisFlows uses obspy stream objects for holding and processing seismic data.
# In some cases, obspy.read doesn't  provide the desired behavior, so we
# introduce an additonal level of indirection
# Used by the PREPROCESS class and specified by the READER parameter
#
###############################################################################

# Import system module
import sys
import os

PAR = sys.modules['seisflows_parameters']


def su(path, filename):
    """ Reads Seismic Unix files.
        Hardwired ''

        Function readBigSuFile is a hack to read su file containing too many
        samples per trace.
        In the su format only 2 bytes per trace are dedicated to encoding for
        the number of samples (as signed int, see:
        http://lists.swapbytes.de/archives/obspy-users/2017-March/002359.html).
        Even if it's an old format it's still extremely stupid.
        This proove the lack of vision the designer of this format had at that
        time. They could have chosen 8 bytes or 16 bytes it was no big deal...
        They've cost me a day's work.
        But let us forget about the past. This limits the size of the
        traces to 32768 samples (NSTEP beween -32768 to 32768). Let us now
        suppose that we have NSTEP = 40000 samples per trace. We still want to
        use Obspy. The problem is that the NSTEP written in the .su file does
        not make any sense anymore and it is read by the obspy.read function!
        We thus rewrote a quick version of this function from Obspy replacing
        the number of point by PAR.NT
    """
    import obspy
    if PAR.NT < 32768:
        stream = obspy.read(path + '/' + filename,
                            format='SU',
                            byteorder='<')
    else:
        stream = readBigSuFile(path + '/' + filename, PAR.NT,
                               format='SU',
                               byteorder='<')
    return stream


def ascii(path, filenames):
    """ Reads SPECFEM3D-style ascii data
    """
    from numpy import loadtxt
    from obspy.core import Stream, Stats, Trace

    stream = Stream()
    for filename in filenames:
        stats = Stats()
        data = loadtxt(path + '/' + filename)

        stats.filename = filename
        stats.starttime = data[0, 0]
        stats.sampling_rate = data[0, 1] - data[0, 0]
        stats.npts = len(data[:, 0])

        try:
            parts = filename.split('.')
            stats.network = parts[0]
            stats.station = parts[1]
            stats.channel = temp[2]
        except:
            pass

        stream.append(Trace(data=data[:, 1], header=stats))

    return stream


def readBigSuFile(nameOfFile, nt, format='SU', byteorder='<'):
    """ This function is a hack to read .su file containing too many samples
        per traces.
        In the su format only 2 bytes per trace are dedicated to encoding for
        the number of samples (as signed int, see:
        http://lists.swapbytes.de/archives/obspy-users/2017-March/002359.html).
        Even if it's an old format it's still extremely stupid.
        This proove the lack of vision the designer of this format had at that
        time. They could have chosen 8 bytes or 16 bytes it was no big deal...
        They've cost me a day's work.
        But let us forget about the past. This limits the size of the
        traces to 32768 samples (NSTEP beween -32768 to 32768). Let us now
        suppose that we have NSTEP = 80000 samples per trace. We still want to
        use Obspy. The problem is that the NSTEP written in the .su file does
        not make any sense anymore and it is read by the obspy.read function!
        We thus rewrote a quick version of this function replacing the number
        of point by PAR.NT
        It is mainly copy-pasted from Obspy source code.
    """

    from obspy.core import Stream, Trace
    from obspy.core.utcdatetime import UTCDateTime
    from obspy.core.util import AttribDict
    from obspy.io.segy.core import LazyTraceHeaderAttribDict
    from obspy.io.segy.segy import SEGYTraceReadingError
    from obspy.io.segy.segy import SEGYTraceHeaderTooSmallError
    from obspy.io.segy.segy import SEGYTraceHeader
    from obspy.io.segy.header import DATA_SAMPLE_FORMAT_UNPACK_FUNCTIONS
    from obspy.io.segy.header import DATA_SAMPLE_FORMAT_SAMPLE_SIZE

    file_object = open(nameOfFile, 'r')
    endian = byteorder
    datas = []
    headers = []
    data_encoding = 5
    # Big loop to read all data traces.
    while True:
        # Read and as soon as the trace header is too small abort.
        try:
            # Always unpack with IEEE
            filesize = os.fstat(file_object.fileno())[6]
            trace_header = file_object.read(240)
            # Check if it is smaller than 240 byte.
            if len(trace_header) != 240:
                msg = 'The trace header needs to be 240 bytes long'
                raise SEGYTraceHeaderTooSmallError(msg)
            header = SEGYTraceHeader(trace_header, endian=endian)
            header.number_of_samples_in_this_trace = nt
            # The number of samples in the current trace.
            npts = header.number_of_samples_in_this_trace
            # Do a sanity check if there is enough data left.
            pos = file_object.tell()
            data_left = filesize - pos
            data_needed = DATA_SAMPLE_FORMAT_SAMPLE_SIZE[data_encoding] * npts
            if npts < 1 or data_needed > data_left:
                msg = """
                      Too little data left in the file to unpack it according
                      to its trace header. This is most likely either due to a
                      wrong byte order or a corrupt file.
                      """.strip()
                raise SEGYTraceReadingError(msg)
            else:
                # Unpack the data
                dsfuf = DATA_SAMPLE_FORMAT_UNPACK_FUNCTIONS
                data = dsfuf[data_encoding](file_object, npts, endian=endian)
            datas.append(data)
            headers.append(header)
        except SEGYTraceHeaderTooSmallError:
            break

    file_object.close()
    # Create the stream object.
    stream = Stream()

    # Loop over all traces.
    for idx, data in enumerate(datas):
        # Create new Trace object for every segy trace and append to the Stream
        # object.
        trace = Trace()
        stream.append(trace)
        trace.data = data
        trace.stats.su = AttribDict()

        # Add the trace header as a new lazy attrib dictionary.
        header = LazyTraceHeaderAttribDict(headers[idx].unpacked_header,
                                           endian)
        trace.stats.su.trace_header = header
        # Also set the endianness.
        trace.stats.su.endian = endian
        # The sampling rate should be set for every trace. It is a sample
        # interval in microseconds. The only sanity check is that is should be
        # larger than 0.
        tr_header = trace.stats.su.trace_header
        if tr_header.sample_interval_in_ms_for_this_trace > 0:
            trace.stats.delta = \
                float(headers[idx].sample_interval_in_ms_for_this_trace) / 1E6
        # If the year is not zero, calculate the start time. The end time is
        # then calculated from the start time and the sampling rate.
        # 99 is often used as a placeholder.
        if tr_header.year_data_recorded > 0:
            year = tr_header.year_data_recorded
            # The SEG Y rev 0 standard specifies the year to be a 4 digit
            # number.  Before that it was unclear if it should be a 2 or 4
            # digit number. Old or wrong software might still write 2 digit
            # years. Every number <30 will be mapped to 2000-2029 and every
            # number between 30 and 99 will be mapped to 1930-1999.
            if year < 100:
                if year < 30:
                    year += 2000
                else:
                    year += 1900
            julday = tr_header.day_of_year
            julday = tr_header.day_of_year
            hour = tr_header.hour_of_day
            minute = tr_header.minute_of_hour
            second = tr_header.second_of_minute
            trace.stats.starttime = UTCDateTime(
                year=year, julday=julday, hour=hour, minute=minute,
                second=second)

    # set _format identifier for each element
    for trace in stream:
        trace.stats._format = format

    return stream
