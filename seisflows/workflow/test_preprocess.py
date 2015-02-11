from seisflows.tools.config import ConfigObj, ParameterObj

OBJ = ConfigObj('SeisflowsObjects')
PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import preprocess


class test_preprocess(object):
    """ Preprocess integration test

        Not yet implemented. The following is just a sketch. None of the 
        methods work yet.
    """

    def check(self):
        """ Checks parameters and paths
        """
        raise NotImplementedError

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
        d = preprocess.load(prefix=PAR.OBSERVATIONS)
        s = preprocess.load(prefix=PAR.SYNTHETICS)
    except:
        print 'LOAD failed'
    else:
        print 'LOAD succeeded'

    try:
        d = preprocess.process_traces(d)
        s = preprocess.process_traces(s)
    except:
        print 'PROCESS_TRACES failed'
    else:
        print 'PROCESS_TRACES succeeded'




