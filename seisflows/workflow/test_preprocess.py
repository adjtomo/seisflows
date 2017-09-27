
import sys
import traceback
import numpy as np

from glob import glob
from os.path import abspath, basename, dirname, exists
from seisflows.config import ParameterError
from seisflows.workflow.base import base

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

preprocess = sys.modules['seisflows_preprocess']


class test_preprocess(base):
    """ Signal processing integration test
    """

    def check(self):
        """ Checks parameters and paths
        """
        # data file format
        if 'FORMAT' not in PAR:
            raise ParameterError(PAR, 'FORMAT')

        # data normalization option
        if 'NORMALIZE' not in PAR:
            setattr(PAR, 'NORMALIZE', None)

        # data muting option
        if 'MUTE' not in PAR:
            setattr(PAR, 'MUTE', None)

        # data filtering option
        if 'FILTER' not in PAR:
            setattr(PAR, 'FILTER', None)


        if 'DATA' not in PATH:
            raise Exception

        if not exists(PATH.DATA):
            raise Exception

        if 'SYNTHETICS' not in PATH:
            setattr(PATH, 'SYNTHETICS', '')

        if PATH.SYNTHETICS:
            assert exists(PATH.SYNTHETICS)

        if 'WORKDIR' not in PATH:
            setattr(PATH, 'WORKDIR', abspath('.'))


    def main(self):
        """ Tests data processing methods
        """

        print 'testing reader...'
        data = self.test_reader()


        print 'testing writer...'
        self.test_writer(data)


        if PAR.NORMALIZE:
            print 'testing normalizing...'
            self.test_normalize(data)


        if PAR.FILTER:
            print 'testing filtering...'
            self.test_filter(data)


        if PAR.MUTE:
            print 'testing muting...'
            self.test_mute(data)


        if PAR.MISFIT and \
           PATH.DATA and \
           PATH.SYNTHETICS:

            dat = preprocess.reader(dirname(PATH.DATA),
                basename(PATH.DATA))

            syn = preprocess.reader(dirname(PATH.SYNTHETICS),
                basename(PATH.SYNTHETICS))

            print 'testing misfit...'
            self.test_misfit(dat, syn)

            print 'testing adjoint...'
            self.test_adjoint(dat, syn)

        print 'SUCCESS\n'



    def test_reader(self):
        try:
            preprocess.setup()

        except Exception,e:
            print 'setup FAILED\n'
            sys.exit(-1)

        try:
            data = preprocess.reader(dirname(PATH.DATA),
                basename(PATH.DATA))

        except Exception,e:
            print 'reader FAILED'
            sys.exit(-1)

        else:
            print ''
            return data



    def test_writer(self, data):
        try:
            if PAR.FORMAT in ['SU', 'su']:
                extension = '.su'
            else:
                extension = ''

            preprocess.writer(data, PATH.WORKDIR, 'output_data'+extension)


        except Exception,e:
            print 'writer FAILED\n'
            print e.message
            print e.__class__.__name__
            traceback.print_exc(e)
            sys.exit(-1)

        else:
            print ''



    def test_normalize(self, dat):
        try:
            out = preprocess.apply_normalize(dat)

        except Exception,e:
            print 'normalization FAILED\n'
            print e.message
            print e.__class__.__name__
            traceback.print_exc(e)
            sys.exit(-1)

        else:
            self.save(out, 'output_data_normalized')
            print ''



    def test_filter(self, dat):
        try:
            out = preprocess.apply_filter(dat)

        except Exception,e:
            print 'filtering FAILED\n'
            print e.message
            print e.__class__.__name__
            traceback.print_exc(e)
            sys.exit(-1)

        else:
            self.save(out, 'output_data_filtered')
            print ''


    def test_mute(self, dat):
        try:
            out = preprocess.apply_mute(dat)

        except Exception,e:
            print 'muting FAILED\n'
            print e.message
            print e.__class__.__name__
            traceback.print_exc(e)
            sys.exit(-1)

        else:
            self.save(out, 'output_data_muted')
            print ''


    def test_misfit(self, dat, syn):
        nt, dt, _ = preprocess.get_time_scheme(syn)
        nn, _ = preprocess.get_network_size(syn)

        rsd = []
        for ii in range(nn):
            rsd.append(preprocess.misfit(syn[ii].data, dat[ii].data, nt, dt))


        filename = PATH.WORKDIR+'/'+'output_misfit'
        np.savetxt(filename, rsd)

        print ''


    def test_adjoint(self, dat, syn):
        nt, dt, _ = preprocess.get_time_scheme(syn)
        nn, _ = preprocess.get_network_size(syn)

        adj = syn
        for ii in range(nn):
            adj[ii].data = preprocess.adjoint(syn[ii].data, dat[ii].data, nt, dt)

        self.save(adj, 'output_adjoint')

        print ''


    def save(self, data, filename):
        if PAR.FORMAT in ['SU', 'su']:
            extension = '.su'
        else:
            extension = ''

        preprocess.writer(data, PATH.WORKDIR, filename+extension)

