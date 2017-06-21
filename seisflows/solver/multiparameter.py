
import sys

from os.path import join

from seisflows.plugins import io
from seisflows.tools.seismic import ModelDict

from seisflows.tools import unix
from seisflows.tools.tools import exists
from seisflows.config import ParameterError, custom_import

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class multiparameter(custom_import('solver', 'base')):
    """ Adds multiparameter inversion machinery

      SOLVER.BASE does already offer some multiparameter capabilites, but these
      are limited to isotropic elastic inversions parameterized in terms of 
      rho,vp,vs or other very simple cases.

      SOLVER.MULITPARAMETER provides a much more flexible framework able to
      accomadate
         - VTI,HTI,TTI inversions
         - fully anisotropic elastic inversions
         - alternative isotropic elastic parameterizations
    """

    # NOT YET TESTED, PROVIDED FOR INFORMATIONAL PURPOSES ONLY
    raise NotImplementedError

    assert PAR.DENSITY in [
        'Constant'
        'Variable'
        'Gardner']

    assert PAR.MATERIALS in [
        'Acoustic'
        'LegacyAcoustic'
        'Elastic'
        'phi_beta',
        'kappa_mu'
        'lambda_mu',
        'alpha_beta'
        'phi_beta_gardner',
        'kappa_mu_gardner'
        'lambda_mu_gardner',
        'alpha_beta_gardner',
        'phi_beta_gardner2',
        'kappa_mu_gardner2'
        'lambda_mu_gardner2',
        'alpha_beta_gardner2',
        'ChenTromp2d',
        'Voigt2d',
        'Thomsen2d']

    if PAR.MATERIALS == 'phi_beta':
        from seisflows.plugins.materials import phi_beta_forward as forward
        from seisflows.plugins.materials import phi_beta_inverse as inverse
        model_parameters = []
        model_parameters += ['rho']
        model_parameters += ['vp']
        model_parameters += ['vs']
        kernel_parameters = []
        kernel_parameters += ['bulk_c']
        kernel_parameters += ['bulk_beta']

    if PAR.MATERIALS == 'kappa_mu':
        from seisflows.plugins.materials import kappa_mu_forward as forward
        from seisflows.plugins.materials import kappa_mu_inverse as inverse
        model_parameters = []
        model_parameters += ['rho']
        model_parameters += ['vp']
        model_parameters += ['vs']
        kernel_parameters = []
        kernel_parameters += ['kappa']
        kernel_parameters += ['mu']

    if PAR.MATERIALS == 'lambda_mu':
        from seisflows.plugins.materials import lambda_mu_forward as forward
        from seisflows.plugins.materials import lambda_mu_inverse as inverse
        model_parameters = []
        model_parameters += ['rho']
        model_parameters += ['vp']
        model_parameters += ['vs']
        kernel_parameters = []
        kernel_parameters += ['lame1']
        kernel_parameters += ['lame2']

    if PAR.MATERIALS == 'alpha_beta':
        from seisflows.plugins.materials import vp_vs_forward as forward
        from seisflows.plugins.materials import vp_vs_inverse as inverse
        model_parameters = []
        model_parameters += ['rho']
        model_parameters += ['vp']
        model_parameters += ['vs']
        kernel_parameters = []
        kernel_parameters += ['vp']
        kernel_parameters += ['vs']


    if PAR.MATERIALS == 'phi_beta_gardner':
        from seisflows.plugins.materials import phi_beta_gardner_forward as forward
        from seisflows.plugins.materials import phi_beta_gardner_inverse as inverse
        model_parameters = []
        model_parameters += ['rho']
        model_parameters += ['vp']
        model_parameters += ['vs']
        kernel_parameters = []
        kernel_parameters += ['bulk_c']
        kernel_parameters += ['bulk_beta']

    if PAR.MATERIALS == 'kappa_mu_gardner':
        from seisflows.plugins.materials import kappa_mu_gardner_forward as forward
        from seisflows.plugins.materials import kappa_mu_gardner_inverse as inverse
        model_parameters = []
        model_parameters += ['rho']
        model_parameters += ['vp']
        model_parameters += ['vs']
        kernel_parameters = []
        kernel_parameters += ['kappa']
        kernel_parameters += ['mu']

    if PAR.MATERIALS == 'lambda_mu_gardner':
        from seisflows.plugins.materials import lambda_mu_garnder_forward as forward
        from seisflows.plugins.materials import lambda_mu_gardner_inverse as inverse
        model_parameters = []
        model_parameters += ['rho']
        model_parameters += ['vp']
        model_parameters += ['vs']
        kernel_parameters = []
        kernel_parameters += ['lame1']
        kernel_parameters += ['lame2']

    if PAR.MATERIALS == 'alpha_beta_gardner':
        from seisflows.plugins.materials import vp_vs_gardner_forward as forward
        from seisflows.plugins.materials import vp_vs_gardner_inverse as inverse
        model_parameters = []
        model_parameters += ['rho']
        model_parameters += ['vp']
        model_parameters += ['vs']
        kernel_parameters = []
        kernel_parameters += ['vp']
        kernel_parameters += ['vs']


    if PAR.MATERIALS == 'ChenTromp2d':
        from seisflows.plugins.maps import voigt_chentromp_2d as forward
        from seisflows.plugins.maps import chentromp_voigt_2d as inverse
        model_parameters = []
        model_parameters += ['c11']
        model_parameters += ['c13']
        model_parameters += ['c15']
        model_parameters += ['c33']
        model_parameters += ['c35']
        model_parameters += ['c55']
        kernel_parameters = []
        kernel_parameters += ['A']
        kernel_parameters += ['C']
        kernel_parameters += ['N']
        kernel_parameters += ['L']
        kernel_parameters += ['F']

    if PAR.MATERIALS == 'Voigt2d':
        from seisflows.plugins.maps import voigt_voigt_2d as forward
        from seisflows.plugins.maps import voigt_voigt_2d as inverse
        model_parameters = []
        model_parameters += ['c11']
        model_parameters += ['c13']
        model_parameters += ['c15']
        model_parameters += ['c33']
        model_parameters += ['c35']
        model_parameters += ['c55']
        kernel_parameters = []
        kernel_parameters += ['c11']
        kernel_parameters += ['c13']
        kernel_parameters += ['c15']
        kernel_parameters += ['c33']
        kernel_parameters += ['c35']
        kernel_parameters += ['c55']

    if PAR.MATERIALS == 'Thomsen2d':
        from seisflows.plugins.maps import voigt_thomsen_2d as forward
        from seisflows.plugins.maps import thomsen_voigt_2d as inverse
        model_parameters = []
        model_parameters += ['c11']
        model_parameters += ['c13']
        model_parameters += ['c15']
        model_parameters += ['c33']
        model_parameters += ['c35']
        model_parameters += ['c55']
        kernel_parameters = []
        kernel_parameters += ['vp']
        kernel_parameters += ['vs']
        kernel_parameters += ['epsilon']
        kernel_parameters += ['delta']
        kernel_parameters += ['gamma']
        kernel_parameters += ['theta']


    if PAR.DENSITY == 'Variable':
        kernel_parameters += ['rho']



    def load(self, path, parameters=[], prefix='', suffix=''):
        self.check_parameters(parameters)
        dict = ModelDict()
        if parameters == self.kernel_parameters:
            for iproc in range(self.mesh_properties.nproc):
                dict[key] = self.read_slice(path, prefix+key+suffix, iproc)
        else:
            for iproc in range(self.mesh_properties.nproc):
                keys = [prefix+key+suffix for key in self.model_parameters]
                vals = self.read_slice(path, keys, iproc)
                mapped = self.forward(keys, vals)
                for key, val in mapped.items():
                    dict[key] = val
        return dict


    def save(self, dict, path, parameters=[], prefix='', suffix=''):
        self.check_parameters(parameters)
        unix.mkdir(path)
        if parameters == self.kernel_parameters:
            for iproc in range(self.mesh_properties.nproc):
                for key in dict.keys():
                    self.write_slice(
                         dict[key][iproc], path, prefix+key+suffix, iproc)
        else:
            for iproc in range(self.mesh_properties.nproc):
                keys = dict.keys()
                vals = [dict[key][iproc] for key in keys]

                if PAR.DENSITY in ['Constant']:
                      keys += ['rho']
                      vals += self.read_slice(
                          PATH.MODEL_INIT, prefix+'rho'+suffix, iproc)

                mapped = self.inverse(keys, vals)
                for key, val in mapped.items():
                    self.write_slice(val, path, prefix+key+suffix, iproc)


    def export_model(self, path, parameters=None)
        if not parameters:
            parameters = self.model_parameters
        super(multiparameter, self).export_model(path, parameters)


    @property
    def parameters(self):
        return self.kernel_parameters
            

    def check_parameters(self, parameters)
        if parameters: assert \
            parameters == self.kernel_parameters or \
            parameters == self.model_parameters 

