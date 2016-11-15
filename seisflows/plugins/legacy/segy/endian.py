"""
from obspy.io.segy

ObsPy is licensed under the LGPL v3.0, i.e. it is licensed with the GPL v3.0
and the additional set of permissions granted by the LGPL v3.0 license.
"""

import io
import os
from struct import pack, unpack

import numpy as np


from obspy.segy.header import (
    BINARY_FILE_HEADER_FORMAT,
    DATA_SAMPLE_FORMAT_PACK_FUNCTIONS,
    DATA_SAMPLE_FORMAT_SAMPLE_SIZE,
    DATA_SAMPLE_FORMAT_UNPACK_FUNCTIONS, ENDIAN,
    TRACE_HEADER_FORMAT, TRACE_HEADER_KEYS
    )

from obspy.segy.unpack import OnTheFlyDataUnpacker
from obspy.segy.util import unpack_header_value


def endian(file):
    """
    Takes an open file and tries to determine the endianness of a Seismic
    Unix data file by doing some sanity checks with the unpacked header values.

    Returns False if the sanity checks failed and the endianness otherwise.

    It is assumed that the data is written as 32bit IEEE floating points in
    either little or big endian.

    The test currently can only identify SU files in which all traces have the
    same length. It basically just makes a sanity check for various fields in
    the Trace header.
    """
    pos = file.tell()
    if isinstance(file, io.BytesIO):
        file.seek(0, 2)
        size = file.tell()
        file.seek(pos, 0)
    else:
        size = os.fstat(file.fileno())[6]
    if size < 244:
        return False
    # Also has to be a multiple of 4 in length because every header is 400 long
    # and every data value 4 byte long.
    elif (size % 4) != 0:
        return False
    # Jump to the number of samples field in the trace header.
    file.seek(114, 0)
    sample_count = file.read(2)
    interval = file.read(2)
    # Jump to the beginning of the year fields.
    file.seek(156, 0)
    year = file.read(2)
    jul_day = file.read(2)
    hour = file.read(2)
    minute = file.read(2)
    second = file.read(2)
    # Jump to previous position.
    file.seek(pos, 0)
    # Unpack in little and big endian.
    le_sample_count = unpack(b'<h', sample_count)[0]
    be_sample_count = unpack(b'>h', sample_count)[0]
    # Check if both work.
    working_byteorders = []
    if le_sample_count > 0:
        length = 240 + (le_sample_count * 4)
        if (size % length) == 0:
            working_byteorders.append('<')
    if be_sample_count > 0:
        length = 240 + (be_sample_count * 4)
        if (size % length) == 0:
            working_byteorders.append('>')
    # If None works return False.
    if len(working_byteorders) == 0:
        return False
    # Check if the other header values make sense.
    still_working_byteorders = []
    for bo in working_byteorders:
        fmt = ("%sh" % bo).encode('ascii', 'strict')
        this_interval = unpack(fmt, interval)[0]
        this_year = unpack(fmt, year)[0]
        this_julday = unpack(fmt, jul_day)[0]
        this_hour = unpack(fmt, hour)[0]
        this_minute = unpack(fmt, minute)[0]
        this_second = unpack(fmt, second)[0]
        # Make a sanity check for each.
        # XXX: The arbitrary maximum of the sample interval is 10 seconds.
        if this_interval <= 0 or this_interval > 10E7:
            continue
        # Some programs write two digit years.
        if this_year != 0 and (this_year < 1930 or this_year >= 2030) and \
                (this_year < 0 or this_year >= 100):
            continue
        # 9999 is often used as a placeholder
        if (this_julday > 366 or this_julday < 0) and this_julday != 9999:
            continue
        if this_hour > 24 or this_hour < 0:
            continue
        if this_minute > 60 or this_minute < 0:
            continue
        if this_second > 60 or this_second < 0:
            continue
        still_working_byteorders.append(bo)
    length = len(still_working_byteorders)
    if not length:
        return False
    elif length == 1:
        return still_working_byteorders[0]
    else:
        # XXX: In the unlikely case both byte orders pass the sanity checks
        # something else should be checked. Currently it is not.
        msg = """
            Both possible byte orders passed all sanity checks. Please contact
            the ObsPy developers so they can implement additional tests.
            """.strip()
        raise Exception(msg)


