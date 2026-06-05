import unittest
import numpy as np

# pj imports
import phylojunction.distribution.likelihood.ODE.dn_bisse_ode as pjbode
import phylojunction.distribution.likelihood.ODE.dn_mbt_ode as pjmbtode

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestMBTBiSSEODEs(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        cls.n_states = 2

        # parameters
        cls.qs = np.matrix([[-1.1, 0.1], [0.1, -0.95]])
        cls.mus = np.array([.25, .35])
        cls.lambdas = np.array([.75, .5])  # for BiSSE
        cls.b = np.zeros((2, 4))  # for MBT
        cls.b[0][0] = .75
        cls.b[1][3] = .5

        cls.t = 0.0
        cls.dt = 0.00008

    def test_mbt_bisse(self) -> None:
        # d0, d1, e0, e1
        ds_es = np.array([1.0, 1.0, 0.0, 0.0])
        res = \
            pjbode.solve_bisse_ds_es(ds_es,
                                 self.t,
                                 self.t + self.dt,
                                 self.qs,
                                 self.mus,
                                 self.lambdas)

        np.testing.assert_array_almost_equal(
            res,
            np.array([0.9999, 0.9999, 0.00002, 0.000028]),
            decimal=4)

    def test_bisse(self) -> None:
        # d0, d1, e0, e1
        ds_es = np.array([1.0, 1.0, 0.0, 0.0])
        ds_es_buffer = np.zeros(2* self.n_states)
        res =\
            pjmbtode.solve_mbt_ds_es(ds_es,
                                     self.t,
                                     self.t + self.dt,
                                     ds_es_buffer,
                                     self.qs,
                                     self.mus,
                                     self.b)

        np.testing.assert_array_almost_equal(
            res,
            np.array([0.9999, 0.9999, 0.00002, 0.000028]),
            decimal=4)
