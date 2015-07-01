#!/bin/env python

from seisflows.tools.config import SeisflowsObjects, SeisflowsParameters, SeisflowsPaths, \
    loadclass

SeisflowsParameters().load()
SeisflowsPaths().load()

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


# run test
if __name__ == '__main__':

    if 'SYSTEM' not in PAR:
        PAR.SYSTEM = 'serial'

    if 'PREPROCESS' not in PAR:
        PAR.PREPROCESS = None

    if 'SOLVER' not in PAR:
        PAR.SOLVER = None

    if 'POSTPROCESS' not in PAR:
        PAR.POSTPROCESS = None

    if 'OPTIMIZE' not in PAR:
        PAR.OPTIMIZE = None

    if 'WORKFLOW' not in PAR:
        PAR.WORKFLOW = 'test_system'

    SeisflowsObjects().load()
    import system
    import workflow

    system.submit(workflow)

