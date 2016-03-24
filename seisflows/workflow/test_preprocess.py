
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import preprocess


class test_preprocess(object):
    """ Preprocess integration test

        Not yet implemented. The following is just a sketch. None of the 
        methods work yet.
    """

    def check(self):
        """ Checks parameters and paths
        """
        #raise NotImplementedError

        # mute settings
        if 'MUTE' not in PAR:
            setattr(PAR, 'MUTE', False)

        if 'MUTESLOPE' not in PAR:
            setattr(PAR, 'MUTESLOPE', 0.)

        if 'MUTECONST' not in PAR:
            setattr(PAR, 'MUTECONST', 0.)

        # filter settings
        if 'BANDPASS' not in PAR:
            setattr(PAR, 'BANDPASS', False)

        if 'FREQLO' not in PAR:
            setattr(PAR, 'FREQLO', 0.)

        if 'FREQHI' not in PAR:
            setattr(PAR, 'FREQHI', 0.)

        # check paths
        if 'OBSERVATIONS' not in PATH:
            raise Exception

        if 'SYNTHETICS' not in PATH:
            raise Exception

        if 'OUTPUT' not in PATH:
            raise Exception


    def main(self):
        """ Tests data processing methods
        """

    try:
        preprocess.setup()
    except:
        print 'SETUP failed'
    else:
        print 'SETUP succeeded'

    try:
        d, h = preprocess.load(prefix=PATH.OBSERVATIONS)
        s, h = preprocess.load(prefix=PATH.SYNTHETICS)
    except:
        print 'LOAD failed'
    else:
        print 'LOAD succeeded'

    try:
        d = preprocess.process_traces([d], [h]) 
        s = preprocess.process_traces([s], [h]) 
    except:
        print 'PROCESS_TRACES failed'
    else:
        print 'PROCESS_TRACES succeeded'


    try:
        preprocess.save(d, h, prefix=PATH.OBSERVATIONS_PRE)
        preprocess.save(s, h, prefix=PATH.SYNTHETICS_PRE)
    except:
        print 'OUTPUT_TRACES failed'
    else:
        print 'OUTPUT_TRACES succeeded'
