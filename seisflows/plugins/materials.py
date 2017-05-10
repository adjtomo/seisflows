
# please do not remove this module  -- it may be used in a future version of
# seisflows

# also, please leave dummy arguments in place for the time being


from seisflows.tools.tools import Struct



### isotropic maps

def phi_beta_forward(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = Struct()

    vp = input.vp
    vs = input.vs
    rho = input.rho

    kappa = rho*(vp**2.-(4./3.)*vs**2.)

    output.bulk_c = (kappa/rho)**0.5
    output.bulk_beta = vs
    output.rho = rho

    return output


def phi_beta_inverse(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = Struct()

    phi = input.bulk_c
    vs = input.bulk_beta
    rho = input.rho

    kappa = rho*phi**2.
    mu = rho*vs**2.

    output.vp = ((kappa+(4./3.)*mu)/rho)**0.5
    output.vs = vs
    output.rho = rho

    return output


def kappa_mu_forward(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = Struct()

    vp = input.vp
    vs = input.vs
    rho = input.rho

    output.kappa = rho*(vp**2.-(4./3.)*vs**2.)
    output.mu = rho*vs**2.
    output.rho = rho

    return output


def kappa_mu_inverse(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = Struct()

    kappa = input.kappa
    mu = input.mu
    rho = input.rho

    output.vp = ((kappa+(4./3.)*mu)/rho)**0.5
    output.vs = (mu/rho)**0.5
    output.rho = rho

    return output


def lambda_mu_forward(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = Struct()

    vp = input.vp
    vs = input.vs
    rho = input.rho

    output.lame1 = rho*(vp**2. - 2.*vs**2.)
    output.lame2 = rho*vs**2.
    output.rho = rho    

    return output


def lambda_mu_inverse(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = Struct()

    lame1 = input.lame1
    lame2 = input.lame2
    rho = input.rho

    output.vp = ((lame1 + 2.*lame2)/rho)**0.5
    output.vs = (lame2/rho)**0.5
    output.rho = rho

    return output


def vs_forward(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = input
    return output


def vs_inverse(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = input
    return output


def vp_vs_forward(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = input
    return output


def vp_vs_inverse(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = input
    return output




### transverse isotropy

def thomsen_voigt_2d(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = Struct()

    vp = input.vp
    vs = input.vs
    rho = input.rho
    epslion = input.epsilon
    delta = input.delta

    output.c11 = rho * vp**2. * (1. + 2.*epsilon)
    output.c12 = rho * vp**2. * (1. + 2.*epsilon) - 2.*rho * vs**2. * (1. + 2.*epsilon)
    output.c13 = rho * vp**2. * (1. + 2.*delta) - 2.*rho * vs**2.
    output.c33 = rho * vp**2.
    output.c55 = rho * vs**2.

    return output



### rotated transverse istoropy

def tti_voight_2d(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = Struct()

    vp = input.vp
    vs = input.vs
    rho = input.rho
    epsilon = input.epsilon
    delta = input.delta
    theta = input.theta

    c11 = rho * vp**2. * (1. + 2.*epsilon)
    c12 = rho * vp**2. * (1. + 2.*epsilon) - 2.*rho * vs**2. * (1. + 2.*epsilon)
    c13 = rho * vp**2. * (1. + 2.*delta) - 2.*rho * vs**2.
    c33 = rho * vp**2.
    c55 = rho * vs**2.

    sint = sin(PI/180. * theta)
    cost = cos(PI/180. * theta)
    sin2t = sin(2.*PI/180. * theta)
    cos2t = cos(2.*PI/180. * theta)

    output.c11 = c11*cost**4 + c33*sint**4 + 2*c13*sint**2*cost**2 + c55*sin2t**2
    output.c12 = 1.e-6
    output.c13 = c13*(cost**4+sint**4) + (c11+c33)*sint**2*cost**2 - c55*sin2t**2
    output.c15 = ((c11-c13)*cost**2 + (c13-c33)*sint**2)*sint*cost - c55*sin2t*cos2t
    output.c23 = 1.e-6
    output.c25 = 0.0
    output.c33 = c11*sint**4 + c33*cost**4 + 2*c13*sint**2*cost**2 + c55*sin2t**2
    output.c35 = ((c11-c13)*sint**2 + (c13-c33)*cost**2)*sint*cost + c55*sin2t*cos2t
    output.c55 = (c11 - 2.*c13+c33)*sint**2*cost**2 + c55*cos2t**2

    return output



### anisotropic maps

def chentromp_voigt_2d(dummy, keys, vals):
    pass


def voigt_chentromp_2d(dummy, keys, vals):
    pass

def voigt_voigt_2d(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = input
    return output



### density maps

def rho_gardner(dummy, keys, vals):
    input = Struct(zip(keys, vals))

    vp = input.vp
    vs = input.vs

    rho = 310.*vp**0.25

    rho_water = 1050

    thresh = 1.0e-3

    rho[vs < thresh] = rho_water

    return rho


def phi_beta_gardner_forward(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = Struct()

    vp = input.vp
    vs = input.vs
    rho = input.rho

    kappa = rho*(vp**2.-(4./3.)*vs**2.)

    output.rho = rho
    output.bulk_c = (kappa/rho)**0.5
    output.bulk_beta = vs

    return output


def phi_beta_gardner_inverse(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = Struct()

    phi = input.phi
    beta = input.beta
    alpha = (phi**2. + (4./3.)*beta**2.)**0.5

    output.rho = 310.*alpha**0.25
    output.vp = alpha
    output.vs = beta

    return output


def lambda_mu_forward(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = Struct()

    vp = input.vp
    vs = input.vs
    rho = input.rho

    output.lame1 = rho*(vp**2. - 2.*vs**2.)
    output.lame2 = rho*vs**2.
    output.rho = rho

    return output




def kernel_map_alpha_beta(self, model, kernels):
    output = Struct()

    vp = model.vp
    vs = model.vs

    rho_alpha = (1./4.)*310.*vp**(-3./4.)
    rho_beta = 0.

    output.vp = [kernels.vp + rho_alpha*kernels.rho]
    output.vs = [kernels.vs]

    return output



def kernel_map_lambda_mu(self, model, kernels):
    output = Struct()

    vp = model.vp
    vs = model.vs
    rho = model.rho

    rho_kappa = (1./9.)*310.**(8./9.)*(rho*vp**2.)**(-8./9.)
    rho_mu = (4./27.)*310.**(8./9.)*(rho*vp**2.)**(-8./9.)

    output.kappa = [kernels.kappa + rho_kappa*kernels.rho]
    output.mu = [kernels.mu + rho_mu*kernels.rho]
    #output.rho = [kernels.rho]

    return output



def kernel_map_kappa_mu(self, model, kernels):
    output = Struct()

    vp = model.vp
    vs = model.vs
    rho = model.rho

    rho_kappa = (1./9.)*310.**(8./9.)*(rho*vp**2.)**(-8./9.)
    rho_mu = (4./27.)*310.**(8./9.)*(rho*vp**2.)**(-8./9.)

    output.kappa = [kernels.kappa + rho_kappa*kernels.rho]
    output.mu = [kernels.mu + rho_mu*kernels.mu]
    output.rho = [kernels.rho]

    return output



def kernel_map_phi_beta(self, model, kernels):
    output = Struct()

    vp = model.vp
    vs = model.vs
    rho = model.rho

    rho_bulk_c = (1./4.)*310.*(vp/rho)**0.5*(vp**2.)**(-7./8.)
    rho_bulk_beta = (1./3.)*310.*(vs/rho)**0.5*(vp**2.)**(-7./8.)

    output.kappa = [kernels.bulk_c + rho_bulk_c*kernels.rho]
    output.mu = [kernels.bulk_beta + rho_bulk_beta*kernels.rho]
    output.rho = [kernels.rho]

    return output

