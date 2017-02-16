
import sys

from glob import glob
from os.path import basename
from seisflows.config import ParameterError

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

solver = sys.modules['seisflows_solver']
preprocess = sys.modules['seisflows_preprocess']


class test_preprocess_adj(object):
    """ Preprocess integration test

        Not yet implemented. The following is just a sketch. None of the 
        methods work yet.
    """

    def check(self):
        """ Checks parameters and paths
        """
        if 'SOLVER' not in PAR:
            raise ParameterError()

        # check paths
        if 'OBS' not in PATH:
            raise Exception

        if 'SYN' not in PATH:
            raise Exception

        if 'OUTPUT' not in PATH:
            raise Exception


    def main(self):
        """ Tests data processing methods
        """
        try:
            preprocess.setup()
        except:
            print 'setup failed'
        else:
            print 'setup succeeded'

        for name in self.data_filenames:
            obs = preprocess.reader(PATH.OBS, name)
            syn = preprocess.reader(PATH.SYN, name)

        for name in self.data_filenames:

            try:
                obs = preprocess.reader(PATH.OBS, name)
                syn = preprocess.reader(PATH.SYN, name)
            except:
                print 'reader failed'

            preprocess.write_adjoint_traces('.', syn, obs, name)

            try:
                preprocess.write_adjoint_traces('.', syn, obs, name)
            except:
                print 'processing failed'
            else:
                print 'processing succeeded'


    @property
    def data_filenames(self):
        filenames = []
        for fullname in glob(PATH.OBS +'/*'):
            filenames += [basename(fullname)]
        obs = filenames

        filenames = []
        for fullname in glob(PATH.SYN +'/*'):
            filenames += [basename(fullname)]
        syn = filenames

        if obs != syn:
            raise Exception()

        if PAR.FORMAT in ['ASCII', 'ascii']:
            return [filenames]
        else:
            return filenames



