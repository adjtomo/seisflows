
import random

import numpy as np

from seisflows.tools import unix
from seisflows.tools.codetools import Struct, exists
from seisflows.tools.configtools import loadclass, ConfigObj, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')
OBJ = ConfigObj('SeisflowsObjects')


class FwiSourceEncoding(loadclass('workflow','inversion')):
    """ Source encoding subclass
    """

    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        super(FwiSourceEncoding,self).check()

        # check source encoding parameters
        if 'ENCODING' not in PAR:
            setattr(PAR,'ENCODING',1)

        if 'FAILRATE' not in PAR:
            setattr(PAR,'FAILRATE',0.)

        if 'SHIFT' not in PAR and PAR.ENCODING in [3,4]:
            raise Exception

        if PAR.ENCODING in [3,4]:
            PAR.NT_PADDED = PAR.NT + (PAR.NSRC-1)*PAR.SHIFT
        else:
            PAR.NT_PADDED = PAR.NT


        # check dependencies
        if 'preprocess' not in OBJ:
            raise Exception

        if 'system' not in OBJ:
            raise Excpetion

        global preprocess
        import preprocess

        global solver
        import solver


        # assertions
        assert('SourceEncoding' in PAR.SOLVER)
        assert(PAR.PREPROCESS in ['default'])


    def setup(self):
        """ Lays groundwork for inversion
        """
        super(FwiSourceEncoding,self).setup()
        self.multiload(path=PATH.DATA,tag='obs')


    def initialize(self):
        """ Prepares for next model update iteration
        """
        # collect headers
        sinfo = self.prepare_sources()
        rinfo = self.prepare_receivers()

        h = self.headers()[0]
        h.dt = PAR.DT
        h.nt = PAR.NT_PADDED

        zeros = np.zeros(h.nr)
        h = Struct(h.items() + \
            [['sx',zeros],['sy',zeros],['sz',zeros]])

        # combine observations into single 'supergather'
        self.combine(h,sinfo,rinfo,tag='obs')

        # update input files
        solver.write_receivers()

        solver.write_sources(
            sinfo=sinfo,
            mapping=lambda _:range(PAR.NSRC))

        # generate synthetics
        super(FwiSourceEncoding,self).initialize()


    def prepare_sources(self):
        """ Generates source encoding factors
        """
        s = self.headers()

        if PAR.ENCODING == 0:
            ts = np.zeros(PAR.NSRC)
            fs = np.ones(PAR.NSRC)

        elif PAR.ENCODING == 1:
            # binary weights
            ts = np.zeros(PAR.NSRC)
            fs = np.sign(np.random.randn(PAR.NSRC))

        elif PAR.ENCODING == 2:
            # Gaussian weights
            ts = np.zeros(PAR.NSRC)
            fs = np.random.randn(PAR.NSRC)

        elif PAR.ENCODING == 3:
            # "staggered" shifts
            ts = np.arange(0,PAR.NSRC)*PAR.SHIFT*PAR.DT
            random.shuffle(ts)
            fs = np.ones(PAR.NSRC)

        elif PAR.ENCODING == 4:
            # random shifts
            ts = ((PAR.NSRC-1)*PAR.SHIFT*PAR.DT)*np.random.rand(PAR.NSRC)
            fs = np.ones(PAR.NSRC)

        # collect factors
        for i in range(PAR.NSRC):
            s[i].ts = ts[i]
            s[i].fs = fs[i]

        return s


    def prepare_receivers(self):
        """ Generates receiver factors
        """
        if self.iter==1:
            # generate receiver factors
            if PAR.FAILRATE==0:
                rs = np.ones((PAR.NREC,PAR.NSRC))
            else:
                rs = np.random.rand(PAR.NREC,PAR.NSRC)
                rs = (rs > PAR.FAILRATE).astype(int)
            np.savetxt(PATH.GLOBAL+'/'+'rs',rs,'%3d')

        # collect factors
        r = []
        rs = np.loadtxt(PATH.GLOBAL+'/'+'rs')
        for i in range(PAR.NSRC):
            r.append(rs[:,i])
        return r


    def solver_status(self):
        """ Because source encoding factors change each iteration, the adjoint
            simulation is never "ready to go" at the start of a new iteration
        """
        isready = False
        return isready


    ### data processing utilities

    def combine(self,h,sinfo,rinfo,tag='obs',inplace=0):
        """ Combines multiple source gathers into a single supergather
        """
        self.multiload(PATH.DATA,tag)
        obj = globals()[tag]
        sum = Struct()

        # allocate arrays
        for channel in solver.channels:
            sum[channel] = np.zeros((h.nt,h.nr))

        # combine gathers
        for i,key in enumerate(obj.keys):
            sum = preprocess.apply(
                      self.combine_traces,[sum,obj.d[key]],[sinfo[i],rinfo[i]])

        # save results
        preprocess.save(
            sum,h,prefix=PATH.SOLVER+'/'+'000000/traces/'+tag,suffix='.bin')


    def combine_traces(self,sum,d,sinfo,rinfo):
        """ Applies source encoding factors to one source gather and adds it to
            another source gather
        """
        if PAR.FAILRATE > 0:
            # apply "receiver factors"
            for j in range(PAR.NREC):
                d[:,j] = rinfo[j]*d[:,j]

        # apply "source factors"
        imin = int(sinfo.ts/PAR.DT)
        imax = imin+PAR.NT
        sum[imin:imax,:] = sum[imin:imax,:] + sinfo.fs*d[:,:]

        return sum


    def multiload(self,path='',tag='obs',inplace=0,debug=0):
        """ Loads data from multiple sources and keeps it in memory
        """
        if tag in globals():
            if inplace:
                pass
            else:
                return
        else:
            obj=Struct()
            obj.d={}
            obj.h={}
            obj.keys=[]

        # generate keys
        keys = unix.ls(path)
        keys.sort()
        obj.keys = keys[0:PAR.NSRC]

        # load data
        for key in obj.keys:
            obj.d[key],obj.h[key] = preprocess.load(path+'/'+key)
            if debug:
                print key

        globals()[tag] = obj


    def headers(self,tag='obs'):
        obj = globals()[tag]
        return [obj.h[key] for key in obj.keys]
