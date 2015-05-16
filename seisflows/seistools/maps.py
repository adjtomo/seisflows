
# please do not remove this module  -- it may be used in a future version of
# seisflows

# also, please leave dummy arguments in place for the time being


from seisflows.seistools.shared import Struct



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


### VTI maps

def thomsen2voigt(dummy, keys, vals):

    raise NotImplementedError

    c11 = rho * vp**2. * (1. + 2.*epsilon)
    c12 = rho * vp**2. * (1. + 2.*epsilon) - 2.*rho * vs**2. * (1. + 2.*epsilon)
    c13 = rho * vp**2. * (1. + 2.*delta) - 2.*rho * vs**2.
    c33 = rho * vp**2.
    c55 = rho * vs**2.



### TTI maps

def tti2voight(dummy, keys, vals):

    c11 = rho * vp**2. * (1. + 2.*epsilon)
    c12 = rho * vp**2. * (1. + 2.*epsilon) - 2.*rho * vs**2. * (1. + 2.*epsilon)
    c13 = rho * vp**2. * (1. + 2.*delta) - 2.*rho * vs**2.
    c33 = rho * vp**2.
    c55 = rho * vs**2.

    raise NotImplementedError

    c11 = c11*cost**4 + c33*sint**4 + 2*c13*sint**2*cost**2 + c55*sin2t**2
    c12 = TINYVAL
    c13 = c13*(cost**4+sint**4) + (c11+c33)*sint**2*cost**2 - c55*sin2t**2
    c15 = ((c11-c13)*cost**2 + (c13-c33)*sint**2)*sint**2*cost**2 - c55*sin2t*cos2t
    c23 = TINYVAL
    c25 = TINYVAL
    c33 = inf
    c55 = inf



### anisotropic maps

def chentromp2voigt(dummy, keys, vals):
    pass


def voigt2chentromp(dummy, keys, vals):
    pass



### density maps

def rho_gardner(dummy, keys, vals):
    input = Struct(zip(keys, vals))

    output = Struct()
    for key in input:
        if key != 'rho':
            output[key] = input[key]

    vp = input.vp

    output.rho = 0.31*vp**0.25

    return output


### debugging

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



