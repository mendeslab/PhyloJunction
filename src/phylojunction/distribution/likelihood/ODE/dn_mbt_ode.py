import typing as ty
import numpy as np
from scipy.integrate import solve_ivp

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


def mbt_e_eqn(es, qs, mus, b):
    de_dt = mus + np.matmul(qs, es) + np.matmul(b, np.kron(es, es))

    return np.ravel(de_dt)

def mbt_bisse_d_eqn(t, ds_es, ds_es_buffer, qs, mus, b):
    # the simplest logic is to assign directly to ds_es_buffer
    # to prevent allocating more memory, but this makes the
    # numbers be slightly different from those of the pure BiSSE
    # ODEs (in the 5th decimal place!)
    es = ds_es[2:]  # grab last two elements
    ds_es_buffer[2:] = mbt_e_eqn(es, qs, mus, b)
    # another option
    # ds = mbt_e_eqn(es, qs, mus, b)

    ds = ds_es[:2]  # grab first two elements
    ds_es_buffer[:2] = np.ravel(
        np.matmul(qs, ds) + np.matmul(b, np.kron(es, ds))
    ) \
                       + np.matmul(b, np.kron(ds, es))
    # another option
    # es = np.ravel(np.matmul(qs, ds) + np.matmul(b, np.kron(es, ds))) \
    #                    + np.matmul(b, np.kron(ds, es))

    # another option
    # ds_es = np.concatenate([ds, es])

    # yet another option (brings it closer, in the 5th decimal place,
    # to the pure BiSSE numbers)
    # ... precision seems to change because of the simple act of
    # assigning values to arrays
    # for i in range(2):
    #     ds_es_buffer[i] = ds[i]
    #     ds_es_buffer[2+i] = es[i]

    # no new memory allocation, pointer swap!!
    return ds_es_buffer

def solve_mbt_ds_es(ds_es, t_start, t_end, ds_es_buffer, qs, mus, b,
                    verbose: ty.Optional[bool]=False):
    pars = (ds_es_buffer, qs, mus, b)
    t_range = [t_start, t_end]

    ds_es = \
        solve_ivp(mbt_bisse_d_eqn,
                  t_range,
                  ds_es,
                  method='RK45',
                  args=pars,
                  rtol=1e-10,
                  atol=1e-10,
                  t_eval=[0, .00002, .00004, .00006, t_end])
         #, t_eval=[t_start + t_end])

    print(ds_es.t)
    print(ds_es.y)
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
    ds_es_buffer = np.zeros(2 * n_states)

    qs = np.matrix([[-1.1, 0.1], [0.1, -0.95]])
    mus = np.array([.25, .35])
    b = np.zeros((2,4))
    b[0][0] = .75
    b[1][3] = .5

    t = 0.0
    dt = 0.00008

    solve_mbt_ds_es(ds_es, t, t + dt, ds_es_buffer, qs, mus, b, verbose=True)