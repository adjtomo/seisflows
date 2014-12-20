
from seisflows.tools.config import ConfigObj, ParameterObj, loadclass

OBJ = ConfigObj('SeisflowsObjects')
PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


# check NPROC is defined (NPROC denotes number of processers per source)
if 'NPROC' not in PAR:
    raise Exception

# on the tiger cluster, there are 16 processers per node
if 'NPROC_PER_NODE' in PAR:
    assert(PAR.NPROC_PER_NODE == 16)
else:
    PAR.NPROC_PER_NODE = 16

# if number of processers per source exceeds number of processers per node, 
# use tiger_lg_job; otherwise, use tiger_sm_job
if PAR.NPROC > PAR.NPROC_PER_NODE:
    tiger = loadclass('system','tiger_lg_job')
else:
    tiger = loadclass('system','tiger_sm_job')

