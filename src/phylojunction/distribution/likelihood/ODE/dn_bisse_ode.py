import typing as ty
import numpy as np
import math
from scipy.integrate import solve_ivp

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


def bisse_e_eqns(es, qs, mus, lambdas):
    """Return Maddison et al. d/dt E(t) equations.

    Args:
        es (numpy.ndarray): Values of Es for states 0 and 1.
        qs (numpy.ndarray): Two-dimensional array of values of
            transition rate parameters for states 0 and 1.
        mus (numpy.ndarray): Values of extinction rate parameters for
            states 0 and 1.
        lambdas (numpy.ndarray): Two-dimensional array of values of
            birth-rate parameters for states 0 and 1.

    Returns:
        (tuple): Values (floats) of E for states 0 and 1
    """

    e0 = mus[0] - (mus[0] + qs[0,1] + lambdas[0]) * es[0] \
         + qs[0,1] * es[1] + lambdas[0] * es[0] ** 2
    e1 = mus[1] - (mus[1] + qs[1,0] + lambdas[1]) * es[1] \
         + qs[1,0] * es[0] + lambdas[1] * es[1] ** 2

    return e0, e1

def bisse_d_eqns(t, ds_es, qs, mus, lambdas):
    """Return Maddison et al. d/dt D(t) equations.

    Args:
        t (float): Time (not used).
        ds_es (numpy.ndarray): Values of Ds for states 0 and 1,
            followed by values of Es, for states 0 and 1.
        qs (numpy.ndarray): Two-dimensional array of values of
            transition rate parameters for states 0 and 1.
        mus (numpy.ndarray): Values of extinction rate parameters for
            states 0 and 1.
        lambdas (numpy.ndarray): Values of birth-rate parameters for
            states 0 and 1.

    Returns
        (numpy.ndarray): Values (floats) of D for all states
            and of E for all states.
    """

    # doing d/dt D(t)
    es = ds_es[2:] # grab last two elements
    e0, e1 = bisse_e_eqns(es, qs, mus, lambdas)

    # doing d/dt E(t)
    ds = ds_es[:2]  # grab first two elements
    d0 = -(lambdas[0] + mus[0] + qs[0,1]) * ds[0] + \
         qs[0,1] * ds[1] + \
         2 * lambdas[0] * e0 * ds[0]

    d1 = -(lambdas[1] + mus[1] + qs[1,0]) * ds[1] + \
         qs[1,0] * ds[0] + \
         2 * lambdas[1] * e1 * ds[1]

    # NOTE: don't try to update 'y' (dn_es),
    # by assigning d0 to y[0], for example.
    # solve_ivp() seems to use dn_ds to store
    # intermediate values or something, and then
    # doing the above overwrites them

    return d0, d1, e0, e1

def solve_bisse_ds_es(ds_es, t_start, t_end, qs, mus, lambdas,
                      verbose: ty.Optional[bool]=False):
    pars = (qs, mus, lambdas)
    t_range = [t_start, t_end]

    ds_es = \
        solve_ivp(bisse_d_eqns,
                  t_range,
                  ds_es,
                  method='RK45',
                  args=pars,
                  rtol=1e-8,
                  atol=1e-8,
                  t_eval=[t_end])  # OdeResult object!

    ds_es_arr = ds_es.y[:,0]

    if verbose:
        print("From", t_start, "to", t_end)
        print("d0, d1, e0, e1")
        print(ds_es_arr)

    return ds_es_arr


if __name__ == "__main__":
    n_states = 2
    # d0, d1, e0, e1
    ds_es = np.array([1.0, 1.0, 0.0, 0.0])

    qs = np.matrix([[-1.1, 0.1], [0.1, -0.95]])
    mus = np.array([.25, .35])
    lambdas = np.array([.75, .5])

    t = 0.0
    dt = 0.00008

    solve_bisse_ds_es(ds_es, t, t + dt, qs, mus, lambdas, verbose=True)