
# please do not remove this module  -- it may be used in a future version of
# seisflows

# also, please leave dummy arguments in place for the time being


from seisflows.seistools.shared import Struct

def bulk_c_bulk_mu():
    raise NotImplementedError


def bulk_c_bulk_mu():
    raise NotImplementedError


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


def vp_vs_forward(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = input
    return output


def vp_vs_inverse(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = input
    return output


def vs_forward(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = input
    return output


def vs_inverse(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = input
    return output


def rho_gardner(dummy, keys, vals):
    input = Struct(zip(keys, vals))

    output = Struct()
    for key in input:
        if key != 'rho':
            output[key] = input[key]

    vp = input.vp

    output.rho = 0.31*vp**0.25

    return output
