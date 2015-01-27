#!/bin/env python

from seisflows.tools.config import loadclass, loadvars, ConfigObj, ParameterObj, Null

OBJ = ConfigObj('SeisflowsObjects')
PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


# run test
if __name__ == '__main__':

    PAR.update(loadvars('parameters','.'))
    PATH.update(loadvars('paths','.'))

    register=OBJ.register

    system = loadclass('system',PAR.SYSTEM)()
    register('system',system)

    preprocess = Null()
    register('preprocess',preprocess)

    solver = Null()
    register('solver',solver)

    postprocess = Null()
    register('postprocess',postprocess)

    optimize = Null()
    register('optimize',optimize)

    workflow = loadclass('workflow','test_system')()
    register('workflow',workflow)

    system.check()
    workflow.check()
    solver.check()
    optimize.check()
    preprocess.check()
    postprocess.check()

    system.submit(workflow)

