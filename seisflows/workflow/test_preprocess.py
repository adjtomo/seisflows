
from glob import glob
from os import basename, exists

from seisflows.config import   \
    ParameterError

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

import preprocess


class test_preprocess(object):
    """ Preprocess integration test

        Not yet implemented. The following is just a sketch. None of the 
        methods work yet.
    """

    def check(self):
        """ Checks parameters and paths
        """
        if 'DATA' not in PATH:
            raise Exception

        if not exists(PATH.DATA):
            raise Exception

        if 'OUTPUT' in PATH:
            assert exists(PATH.OUTPUT)


    def main(self):
        """ Tests data processing methods
        """
        data = {}

        try:
            preprocess.setup()
        except:
            print 'setup failed'
        else:
            print 'setup succeeded'

        # test reader
        try:
            for channel in self.channels():
                data[channel] = preprocess.reader(PATH.DATA, channel)
        except:
            print 'reader failed'
        else:
            print 'reader succeeded'

        # test processing
        try:
            for channel in self.channels():
                data[channel] = preprocess.apply_filter(data[channel])
        except:
            print 'processing failed'
        else:
            print 'processing succeeded'

        try:
            for channel in self.channels():
                preprocess.writer(data[channel], PATH.OUTPUT, channel)
        except:
            print 'writer failed'
        else:
            print 'writer succeeded'


    @property
    def channels(self):
        channels = []
        for fullname in glob(PATH.DATA):
            channels += [basename(fullname)]

        return channels


