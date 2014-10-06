
from seisflows.tools.configure import getclass, ParameterObject

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')

system = getclass('system',PAR.SYSTEM)()


class specfem3d_Thomsen(getclass('extensions.solver','specfem3d_legacy')):

    # model parameters expected by solver
    model_parameters = []
    model_parameters += ['rho']
    model_parameters += ['vp']
    model_parameters += ['vs']
    model_parameters += ['epsilon']
    model_parameters += ['delta']
    model_parameters += ['gamma']
    model_parameters += ['theta']
    model_parameters += ['azimuth']

    # model parameters included in inversion
    inversion_parameters = []
    model_parameters += ['rho']
    inversion_parameters += ['vp']
    inversion_parameters += ['vs']
    inversion_parameters += ['epsilon']
    inversion_parameters += ['delta']
    inversion_parameters += ['gamma']
    inversion_parameters += ['theta']
    inversion_parameters += ['azimuth']

    # data channels
    channels = []
    channels += ['x']

    kernel_map = {
        'rho':'rho_kernel',
        'vp':'alpha_kernel',
        'vs':'beta_kernel',
        'epsilon':'epsilon_kernel',
        'delta':'delta_kernel',
        'gamma':'gamma_kernel',
        'theta':'theta_kernel',
        'azimuth':'azimuth_kernel'}
