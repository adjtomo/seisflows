
# please do not remove this module  -- it may be used in a future version of
# seisflows

# also, please leave dummy arguments in place for the time being


from seisflows.seistools.shared import Struct

def forward_bulk_c_bulk_mu():
    raise NotImplementedError


def inverse_bulk_c_bulk_mu():
    raise NotImplementedError


def forward_kappa_mu(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = Struct()

    vp = input.vp
    vs = input.vs
    rho = input.rho

    output.kappa = rho*(vp**2.-(4./3.)*vs**2.)
    output.mu = rho*vs**2.
    output.rho = rho

    return output


def inverse_kappa_mu(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = Struct()

    kappa = input.kappa
    mu = input.mu
    rho = input.rho

    output.vp = ((kappa+(4./3.)*mu)/rho)**0.5
    output.vs = (mu/rho)**0.5
    output.rho = rho

    return output


def forward_vp_vs(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = input
    return output


def inverse_vp_vs(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = input
    return output


def forward_vs(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = input
    return output


def inverse_vs(dummy, keys, vals):
    input = Struct(zip(keys, vals))
    output = input
    return output




