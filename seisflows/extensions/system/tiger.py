
from seisflows.tools.config import ConfigObj, ParameterObj, loadclass

OBJ = ConfigObj('SeisflowsObjects')
PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


# ensure that NPROC (number of processers per source) is defined
if 'NPROC' not in PAR:
    raise Exception

# on tiger.princeton.edu, there are 16 processers per node
if 'NPROC_PER_NODE' in PAR:
    assert(PAR.NPROC_PER_NODE == 16)
else:
    PAR.NPROC_PER_NODE = 16

# if processers per source exceeds processers per node, use tiger_lg_job, 
# otherwise, use tiger_sm_job
if PAR.NPROC > PAR.NPROC_PER_NODE:
    tiger = loadclass('system','tiger_lg_job')
else:
    tiger = loadclass('system','tiger_sm_job')

