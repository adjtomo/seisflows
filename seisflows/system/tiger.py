
from seisflows.tools.config import loadclass
from seisflows.tools.config import ParameterError, SeisflowsParameters, SeisflowsPaths

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


# ensure number of processers per source is defined
if 'NPROC' not in PAR:
    raise Exception

# there are 16 processers per node on tiger
if 'NODESIZE' in PAR:
    assert(PAR.NODESIZE == 16)
else:
    PAR.NODESIZE = 16

# if nproc per source exceeds nproc per node, use tiger_lg
# otherwise, use tiger_sm
if PAR.NPROC > PAR.NODESIZE:
    tiger = loadclass('system','tiger_lg')
else:
    tiger = loadclass('system','tiger_sm')

