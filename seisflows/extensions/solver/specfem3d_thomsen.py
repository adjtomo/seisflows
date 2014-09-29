
from seisflows.tools.configure import getclass, ParameterObject

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')

system = getclass('system',PAR.SYSTEM)()


class specfem3d_thomsen(getclass('extensions.solver','specfem3d_legacy')):

  def configure_model(self):

    # model parameters expected by solver
    model_parameters = []
    model_parameters += ['rho']
    model_parameters += ['vp']
    model_parameters += ['vs']
    model_parameters += ['epsilon']
    model_parameters += ['delta']
    model_parameters += ['gamma']
    self.model_parameters = model_parameters

    # model parameters included in inversion
    inversion_parameters = []
    inversion_parameters += ['vp']
    inversion_parameters += ['vs']
    inversion_parameters += ['epsilon']
    inversion_parameters += ['delta']
    inversion_parameters += ['gamma']
    self.inversion_parameters = inversion_parameters

    self.kernel_dict = {
        'rho':'rho_kernel',
        'vp':'alpha_kernel',
        'vs':'beta_kernel',
        'epsilon':'epsilon_kernel',
        'delta':'delta_kernel',
        'gamma':'gamma_kernel'}

