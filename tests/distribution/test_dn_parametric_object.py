import sys
sys.path.extend(["../", "../phylojunction"])
import unittest
import math
import numpy as np
from unittest.mock import patch
from phylojunction.distribution.dn_iid import DnIID
from statistics import mean, stdev

# pj imports
import phylojunction.distribution.dn_parametric as dnpar

class TestParametricDns(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dummy_parent_node_tracker = dict()

    # Existing moment checks do not protect nested shape or bulk sampling. This
    # contract matters while IID vectors are explicit values rather than model replicates.
    def test_iid_draws(self):
        distribution = dnpar.DnNormal(2.0, 0.5)
        with patch.object(distribution, 'draw_normal', return_value=np.arange(8.0)) as draw:
            self.assertEqual(DnIID(8, distribution).draw(), list(range(8)))
            draw.assert_called_once_with(8, 2.0, 0.5)
        nested = DnIID(3, DnIID(4, distribution)).draw()
        self.assertEqual([len(row) for row in nested], [4, 4, 4])
        self.assertEqual(DnIID(3, DnIID(0, distribution)).draw(), [[], [], []])
        self.assertEqual(DnIID(0, distribution).draw(), [])
        for n in (-1, 1.5, True):
            with self.assertRaises(ValueError):
                DnIID(n, distribution)
        beta = dnpar.DnBeta(2, 3)
        with patch('phylojunction.distribution.dn_parametric.beta.rvs',
                   return_value=np.array([0.25, 0.75])) as draw:
            self.assertEqual(beta.draw_many(2), [0.25, 0.75])
            draw.assert_called_once_with(2, 3, size=2)
        for shapes in ((0, 1), (1, -1), (float('nan'), 1), (1, float('inf'))):
            with self.assertRaises(ValueError):
                dnpar.DnBeta(*shapes)

    # Cauchy has no moments: quantiles protect parameterization without unstable mean tests.
    # Language checks also protect the object, vectorization, and IID integration paths.
    def test_cauchy(self):
        from phylojunction.interface.cmdbox.cmd_parse import compile_model
        for location, scale in [(0, 0), (0, -1), (float('nan'), 1), (0, float('inf'))]:
            with self.assertRaises(ValueError):
                dnpar.DnCauchy(location, scale)
        state = np.random.get_state()
        try:
            np.random.seed(42)
            values = dnpar.DnCauchy(2, 3).draw_many(20000)
            np.testing.assert_allclose(np.quantile(values, [.25, .5, .75]), [-1, 2, 5], atol=.15)
            self.assertIsInstance(dnpar.DnCauchy(0, 1).draw(), float)
            self.assertEqual(DnIID(0, dnpar.DnCauchy(0, 1)).draw(), [])
            a = compile_model('d := cauchy(2,3)\nx ~ iid(4,d)').execute(2,42)
            b = compile_model('x ~ iid(4,cauchy(location=2,scale=3))').execute(2,42)
            self.assertEqual(a.values('x'), b.values('x'))
            a = compile_model('x ~ cauchy([0,1],1)').execute(1,42)
            b = compile_model('x ~ [cauchy(mu,1) for mu in [0,1]]').execute(1,42)
            self.assertEqual(a.values('x'), b.values('x'))
        finally:
            np.random.set_state(state)

    def test_unif(self):
        """Test uniform draws with non-vectorized input."""

        unif_rv = dnpar.DnUnif(-1.0, 1.0).draw_many((100000) * (1))
        self.assertAlmostEqual(0.0, mean(unif_rv), delta=1e-2)
        self.assertLessEqual(-1.0, min(unif_rv))
        self.assertGreater(1.0, max(unif_rv))

    def test_unif_vectorized(self):
        """Test uniform draws with vectorized input"""
        repl_size = 2
        mins = [0.0, 1.0, 2.0, 3.0, 4.0]
        maxs = [1.0, 2.0, 3.0, 4.0, 5.0]
        tups = tuple((i, maxs[idx]) for idx, i in enumerate(mins))
        # unif_rv = dnpar.DnUnif([5, repl_size, mins, maxs]).generate()
        unif_rv = [value for lo, hi in zip(mins, maxs) for value in dnpar.DnUnif(lo, hi).draw_many(repl_size)]

        for idx, tup in enumerate(tups):
            self.assertTrue(tup[0] <= unif_rv[idx * repl_size] < tup[1])
            self.assertTrue(tup[0] <= unif_rv[idx * repl_size + 1] < tup[1])

    def test_exp(self):
        """Test exponential draws with non-vectorized input"""
        # exp_rv1 = dnpar.DnExponential([100000, 1, 0.5, True]).generate() # rate parameterization (True is default)
        # exp_rv2 = dnpar.DnExponential([100000, 1, 0.5, False]).generate()
        exp_rv1 = dnpar.DnExponential(0.5, True).draw_many((100000) * (1)) # rate parameterization (True is default)
        exp_rv2 = dnpar.DnExponential(0.5, False).draw_many((100000) * (1))
        self.assertAlmostEqual(2.0, mean(exp_rv1), delta=0.05)
        self.assertAlmostEqual(0.5, mean(exp_rv2), delta=0.05)

    def test_gamma(self):
        """Test gamma draws with non-vectorized input"""
        # gamma_rv1 = dnpar.DnGamma([100000, 1, [0.5], [0.5], False]).generate() # rate parameterization (False is default)
        # gamma_rv2 = dnpar.DnGamma([100000, 1, [0.5], [0.5], True]).generate() # rate parameterization (False is default)
        gamma_rv1 = dnpar.DnGamma(0.5, 0.5, False).draw_many((100000) * (1)) # rate parameterization (False is default)
        gamma_rv2 = dnpar.DnGamma(0.5, 0.5, True).draw_many((100000) * (1)) # rate parameterization (False is default)
        self.assertAlmostEqual(0.25, mean(gamma_rv1), delta=0.05)
        self.assertAlmostEqual(1.0, mean(gamma_rv2), delta=0.05)

    def test_normal(self):
        """Test normal draws with non-vectorized input"""
        # normal_rv = dnpar.DnNormal([100000, 1, [0.5], [0.1]]).generate()
        normal_rv = dnpar.DnNormal(0.5, 0.1).draw_many((100000) * (1))
        self.assertAlmostEqual(0.5, mean(normal_rv), delta=0.05)

    def test_lognormal(self):
        """Test log-normal draws with non-vectorized input."""

        # log-space = True, same as R's default
        lognormal_rv1 = dnpar.DnLogNormal(-3.25, 0.25, True).draw_many((200000) * (1))
        # real space
        lognormal_rv2 = dnpar.DnLogNormal(math.exp(-3.25), 0.25, False).draw_many((200000) * (1))
        self.assertAlmostEqual(0.04, mean(lognormal_rv1), delta=0.01)
        self.assertAlmostEqual(0.04, mean(lognormal_rv2), delta=0.01)
        self.assertAlmostEqual(0.01, stdev(lognormal_rv1), delta=0.01)
        self.assertAlmostEqual(0.01, stdev(lognormal_rv2), delta=0.01)


if __name__ == '__main__':
    # Assuming you opened the PhyloJunction/ (repo root) folder
    # on VSCode and that you want to run this as a standalone script,
    # i.e., "Run Without Debugging", you will need to configure your
    # launch.json to have:
    # 
    # "env": {"PYTHONPATH": "${workspaceRoot}/src/phylojunction/"}
    # 
    # and your settings.json to have:
    #   
    # "python.analysis.extraPaths": [ "${workspaceFolder}/src/phylojunction/" ]
    # 
    # If you want to run this as a standalone from PhyloJunction/
    # on the terminal, remember to add "src/phylojunction" to
    # PYTHONPATH (system variable), or to set it if it does not
    # exist -- don't forget to export it!
    # 
    # Then you can do:
    # $ python3 tests/distribution/test_dn_parametric_object.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_dn_parametric_object
    #
    # or 
    #
    # $ python3 -m unittest tests.distribution.test_dn_parametric_object.TestParametricDns.test_unif
    
    unittest.main()