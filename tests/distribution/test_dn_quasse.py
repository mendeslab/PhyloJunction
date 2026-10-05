"""Deterministic checks of QuaSSE's numerical and observation contracts."""

import unittest
from unittest.mock import patch
import numpy as np

from phylojunction.calculation.continuous_sse import (
    LogisticRate, ConstantRate, GaussianRate, StepRate, LinearRate,
    SkewGaussianRate, QuadraticRate, AbsoluteRate,
)
from phylojunction.data.trait import ContinuousTrait
from phylojunction.distribution.dn_quasse import DnQuaSSE
from phylojunction.interface.cmdbox.cmd_parse import cmdline2dag
from phylojunction.pgm.pgm import DirectedAcyclicGraph
from phylojunction.utility import exception_classes as ec


class TestQuaSSE(unittest.TestCase):
    # Direct checks cover rate shapes and caps that stochastic tree tests cannot isolate.
    # Retire if the rate families move to another library with equivalent mathematical tests.
    def test_rate_families(self):
        cases = [(ConstantRate, (2,), [2, 2, 2]),
                 (LogisticRate, (1, 3, 0, 0), [2, 2, 2]),
                 (GaussianRate, (2, 0, 0, 1), [2 * (1 - np.exp(-.5)), 0, 2 * (1 - np.exp(-.5))]),
                 (StepRate, (1, 3, 0), [1, 3, 3]),
                 (LinearRate, (0, 2), [0, 0, 2]),
                 (SkewGaussianRate, (1, 2, 0, 1, 0), [1 + 2 * np.exp(-.5), 3, 1 + 2 * np.exp(-.5)]),
                 (QuadraticRate, (1, 2, 0), [3, 1, 3]),
                 (AbsoluteRate, (1, 2, 0), [3, 1, 3])]
        for cls, args, expected in cases:
            with self.subTest(rate=cls.__name__):
                rate = cls(*args)
                np.testing.assert_allclose(rate([-1, 0, 1]), expected)
                self.assertAlmostEqual(float(rate(0)), expected[1])
                capped = cls(*args, cap=.5)
                np.testing.assert_allclose(capped([-1, 0, 1]), np.minimum(expected, .5))
                self.assertLessEqual(capped.upper_bound, .5)
                self.assertIn('cap=0.5', str(capped))
                self.assertTrue(cls(*args, cap=0).is_identically_zero)
                self.assertEqual(cls(*args, cap=0).bound_on_interval(-1, 1), 0)
                self.assertAlmostEqual(capped.bound_on_interval(-1, 1), .5)
                for cap in (-1, np.inf, np.nan):
                    with self.assertRaises(ValueError):
                        cls(*args, cap=cap)
        x = np.array([-2., 0., 2.])
        np.testing.assert_allclose(SkewGaussianRate(1, 2, 0, 1, 3)(x),
                                   SkewGaussianRate(1, 2, 0, 1, -3)(-x))
        for rate in (LinearRate(0, 1), QuadraticRate(0, 1, 0), AbsoluteRate(0, 1, 0)):
            self.assertIsNone(rate.upper_bound)
            self.assertFalse(rate.is_identically_zero)
            self.assertEqual(float(rate(0)), 0)
            self.assertGreater(float(rate(100)), float(rate(10)))
        self.assertTrue(LinearRate(-1, 0).is_identically_zero)
        np.testing.assert_allclose(QuadraticRate(0, 1e-300, 0, cap=1e200)(1e200), 1e100)
        self.assertEqual(float(QuadraticRate(1, 1, 0, cap=3)(1e300)), 3)
        for cls, args in [(ConstantRate, (-1,)), (GaussianRate, (0, 1, 0, 0)),
                          (StepRate, (0, -1, 0)), (LinearRate, (0, np.inf)),
                          (SkewGaussianRate, (0, -1, 0, 1, 0)),
                          (QuadraticRate, (0, -1, 0)), (AbsoluteRate, (0, 1, np.nan))]:
            with self.assertRaises(ValueError):
                cls(*args)

    # Exact maxima supplement grid checks: tree tests can miss a narrow interior peak.
    # Retire this group if rate bounds are supplied by another verified implementation.
    def test_interval_bounds(self):
        cases = [(ConstantRate(2), -1, 1, 2),
                 (LogisticRate(3, 1, 0, 2), -1, 1, float(LogisticRate(3, 1, 0, 2)(-1))),
                 (GaussianRate(1, 5, .1, .01), -1, 1, 5),
                 (GaussianRate(5, 1, 0, 1), -1, 2, float(GaussianRate(5, 1, 0, 1)(2))),
                 (StepRate(3, 1, 0), 0, 1, 1), (StepRate(1, 3, 0), -1, 0, 3),
                 (LinearRate(1, -2), -1, 1, 3), (QuadraticRate(1, 2, 0), -1, 2, 9),
                 (AbsoluteRate(1, 2, 0), -1, 2, 5), (QuadraticRate(2, 0, 0), -1, 2, 2),
                 (SkewGaussianRate(1, 2, 0, 1, 3), -1, 2, 5),
                 (QuadraticRate(3, 2, 0, cap=1), -1, 2, 1)]
        for rate, left, right, expected in cases:
            with self.subTest(rate=str(rate)):
                bound = rate.bound_on_interval(left, right)
                self.assertAlmostEqual(bound, expected)
                self.assertTrue(np.all(rate(np.linspace(left, right, 101)) <= bound))
                self.assertGreaterEqual(rate.bound_on_interval(left, left), float(rate(left)))
        for left, right in [(1, 0), (np.nan, 1), (0, np.inf)]:
            with self.assertRaises(ValueError):
                ConstantRate(1).bound_on_interval(left, right)
        with patch.object(QuadraticRate, '_interval_bound', return_value=float('inf')):
            with self.assertRaises(ValueError):
                QuadraticRate(0, 1, 0).bound_on_interval(0, 1)

    # Discrete tests cannot detect cancellation in positive logistic tails or resulting simulation failures.
    # Retire this group if another shared rate implementation tests this same contract.
    def test_logistic(self):
        rate = LogisticRate(1, 3, 0, 2)
        np.testing.assert_allclose(rate([-1e300, 0, 1e300]), [1, 2, 3])
        np.testing.assert_allclose(LogisticRate(1, 3, 0, -2)([-1e300, 0, 1e300]), [3, 2, 1])
        np.testing.assert_allclose(LogisticRate(2, 2, 0, 2)([-1, 0, 1]), 2)
        np.testing.assert_allclose(LogisticRate(1, 3, 1e308, 0)([-1e308]), 2)
        zero = LogisticRate(0, 0, 0, 0)
        for slope in (-1, 1):
            for y0, y1, x in ((1, 0, 40 * slope), (0, 1, -40 * slope)):
                tail = LogisticRate(y0, y1, 0, slope)
                np.testing.assert_allclose(tail(x), 1 / (1 + np.exp(40)), rtol=1e-14, atol=0)
                with patch('numpy.random.random', return_value=.5):
                    tree = DnQuaSSE(tail, zero, 'age', 1, start_trait=x).generate()[0]
                self.assertEqual(tree.n_extant_terminal_nodes, 1)
        for values in [(-1, 2, 0, 0), (0, 1, np.inf, 0)]:
            with self.assertRaises(ValueError):
                LogisticRate(*values)

    # Supplied random outcomes detect event-order, inheritance, and final-step regressions.
    # These checks are needed until the diversitree method is removed or replaced.
    def test_final_step_birth_and_death(self):
        zero, one = LogisticRate(0, 0, 0, 0), LogisticRate(1, 1, 0, 0)
        dn = DnQuaSSE(one, zero, 'age', .25, k=2, drift=2, diffusion=4)
        with patch('numpy.random.random', side_effect=[.2, 0]), \
                patch('numpy.random.choice', return_value=0), \
                patch('numpy.random.normal', return_value=np.array([1., -1.])) as normal:
            tree = dn.simulate()
        normal.assert_called_once_with(.5, 1., 2)
        self.assertEqual(tree.root_node.edge_length, 0)
        self.assertEqual(tree.root_node.trait, 0)
        self.assertEqual([nd.trait for nd in tree.tree.leaf_node_iter()], [1., -1.])
        self.assertEqual([nd.edge_length for nd in tree.tree.leaf_node_iter()], [.25, .25])
        # .3 would cause an event at the original 1/k probability, but not at R*delta.
        with patch('numpy.random.random', return_value=.3):
            self.assertIsNone(dn.simulate().root_node)
        death = DnQuaSSE(zero, one, 'age', 1, cond_surv=False, k=1)
        with patch('numpy.random.random', return_value=0), patch('numpy.random.choice', return_value=0):
            tree = death.generate()[0]
        self.assertTrue(tree.tree_died)
        self.assertEqual(tree.node_ages_dict['brosc'], 1.)
        self.assertEqual(tree.brosc_node.edge_length, 0.)
        self.assertEqual(tree.n_extinct_terminal_nodes, 1)

    # Size stopping must omit both the terminating birth and its following trait update.
    # Age-only tests cannot detect this; retain while the N+1 stopping convention exists.
    def test_size_stopping(self):
        one, zero = LogisticRate(1, 1, 0, 0), LogisticRate(0, 0, 0, 0)
        dn = DnQuaSSE(one, zero, 'size', 2, k=1, drift=1, cond_obs_both_sides=True)
        with patch('numpy.random.random', return_value=0), patch('numpy.random.choice', return_value=0):
            tree = dn.generate()[0]
        self.assertEqual(tree.seed_age, 1)
        self.assertEqual(tree.n_extant_terminal_nodes, 2)
        self.assertEqual([nd.trait for nd in tree.tree.leaf_node_iter()], [1, 1])
        with patch('numpy.random.random', return_value=0), patch('numpy.random.choice', return_value=0):
            single = DnQuaSSE(one, zero, 'size', 1, k=1).generate()[0]
        self.assertEqual(single.seed_age, 0)
        self.assertEqual(single.n_extant_terminal_nodes, 1)

    # DAG integration protects vector broadcasting and retry-independent parameter ownership.
    # Discrete rate objects follow a different path; retire if shared object-vector tests replace it.
    def test_dag_and_zero_rates(self):
        dag = DirectedAcyclicGraph()
        for line in ['x <- [0,2]',
                     'b := quasse_logistic(y0=0,y1=0,midpoint=x,slope=1)',
                     't ~ quasse(n=2,nr=2,birth_rate=b,death_rate=b,stop="age",'
                     'stop_value=[0,2],start_trait=x,drift=1,sampling_prob=[1,0])']:
            cmdline2dag(dag, line)
        trees = dag.name_node_dict['t'].value
        self.assertEqual([t.seed_age for t in trees], [0, 0, 2, 2])
        self.assertEqual([t.brosc_node.trait for t in trees], [0, 0, 4, 4])
        self.assertEqual([t.n_extant_sampled_terminal_nodes for t in trees], [1, 1, 0, 0])
        self.assertEqual(trees[0].extract_reconstructed_tree().seed_node.child_nodes()[0].trait, 0)
        self.assertEqual(str(trees[-1].extract_reconstructed_tree()), ';')
        self.assertIn(dag.name_node_dict['b'], dag.name_node_dict['t'].parent_nd_list)
        constructors = ["constant(rate=2", "gaussian(baseline=1,center_rate=3,center=0,width=1",
                        "step(left=1,right=2,threshold=0", "linear(intercept=1,slope=2",
                        "skew_gaussian(baseline=1,amplitude=2,location=0,width=1,skew=3",
                        "quadratic(baseline=1,strength=2,center=0", "absolute(baseline=1,strength=2,center=0"]
        for constructor in constructors:
            cmdline2dag(dag, f"r := quasse_{constructor},cap=x)")
            rates = dag.name_node_dict['r'].value
            self.assertEqual(len(rates), 2)
            self.assertTrue(rates[0].is_identically_zero)
            self.assertIn(dag.name_node_dict['x'], dag.name_node_dict['r'].parent_nd_list)
        cmdline2dag(dag, 'u ~ quasse(n=2,birth_rate=r,death_rate=b,stop="age",stop_value=0,dt_max=.1,max_nodes=2)')
        self.assertEqual(len(dag.name_node_dict['u'].value), 2)
        for spec in ('quasse_constant(rate=[1,2],cap=[1,2,3])', 'quasse_constant(rate=1,cap=-1)',
                     'quasse_gaussian(baseline=1,center_rate=2,center=0)', 'quasse_constant(rate=1,typo=2)'):
            with self.assertRaises(ec.ParseDetFnInitFailError):
                cmdline2dag(dag, f'invalid := {spec}')

    # Frozen-rate steps must cross zero-rate regions and scale event chance with dt_max.
    # Random tree tests cannot isolate either; retire when this integration method is replaced.
    def test_dt_max(self):
        zero = ConstantRate(0)
        dn = DnQuaSSE(StepRate(0, 1, .5), zero, 'age', 1, drift=1, dt_max=.5, k=1)
        with patch('numpy.random.random', side_effect=[0., 0.]), \
                patch('numpy.random.choice', return_value=0):
            tree = dn.simulate()
        self.assertEqual(tree.n_extant_terminal_nodes, 2)
        self.assertEqual(tree.root_node.trait, .5)
        with patch('numpy.random.random', return_value=.3):
            tree = DnQuaSSE(ConstantRate(1), zero, 'age', 1, k=1, dt_max=.25).simulate()
        self.assertEqual(tree.n_extant_terminal_nodes, 1)
        with patch('numpy.random.normal', return_value=np.array([2.])) as normal:
            tree = DnQuaSSE(zero, zero, 'age', 1, diffusion=1, dt_max=.01).simulate()
        normal.assert_called_once_with(0., 1., 1)
        self.assertEqual(tree.brosc_node.trait, 2.)

    # Historical nodes consume memory even after death; max_alive cannot protect that bound.
    # Retire if simulation no longer retains dead nodes or uses a shared allocation limiter.
    def test_max_nodes(self):
        one, zero = ConstantRate(1), ConstantRate(0)
        with patch.object(DnQuaSSE, '_new_node') as allocate:
            with self.assertRaisesRegex(ec.GenerateFailError, 'max_nodes'):
                DnQuaSSE(one, zero, 'age', 1, max_nodes=1).generate()
            allocate.assert_not_called()
        dn = DnQuaSSE(one, one, 'age', 2, k=1, max_nodes=4)
        # Birth, death, birth: only one lineage remains, but the third event needs nodes 5/6.
        with patch('numpy.random.random', side_effect=[0., 0., 0., .9, 0., 0.]), \
                patch('numpy.random.choice', return_value=0), \
                patch.object(dn, '_new_node', wraps=dn._new_node) as allocate:
            with self.assertRaisesRegex(ec.GenerateFailError, 'max_nodes'):
                dn.generate()
            self.assertEqual(allocate.call_count, 4)
        with patch('numpy.random.random', return_value=0.), patch('numpy.random.choice', return_value=0):
            tree = DnQuaSSE(one, zero, 'size', 2, k=1, max_nodes=4).generate()[0]
        self.assertEqual(len(list(tree.tree)), 4)

    # Count and root-side conditions must inspect observations, without changing parameters.
    # These tests supplement numerical checks and become obsolete with a shared acceptance layer.
    def test_sampling_conditions_and_retries(self):
        one, zero = LogisticRate(1, 1, 0, 0), LogisticRate(0, 0, 0, 0)
        dn = DnQuaSSE(one, zero, 'age', .5, k=1, sampling_prob=.5,
                     min_rec_taxa=1, max_rec_taxa=1, cond_spn=True)
        with patch('numpy.random.random', side_effect=[0., 0., .1, .9]), \
                patch('numpy.random.choice', return_value=0):
            tree = dn.simulate()
        self.assertEqual(tree.n_extant_sampled_terminal_nodes, 1)
        self.assertTrue(dn._is_tree_accepted(tree))
        dn.cond_obs_both_sides = True
        self.assertFalse(dn._is_tree_accepted(tree))
        dn.cond_obs_both_sides = False
        rec = tree.extract_reconstructed_tree()
        self.assertEqual(len(list(rec.leaf_node_iter())), 1)
        self.assertEqual(tree.root_node.label, 'root')
        dn.n_sim, dn.n_repl = 1, 2
        with patch.object(dn, 'simulate', return_value=tree) as simulate, \
                patch.object(dn, '_is_tree_accepted', side_effect=[False, True, True]):
            self.assertEqual(len(dn.generate()), 2)
        self.assertEqual([call.args[0] for call in simulate.call_args_list], [0, 0, 0])
        dn.max_n_attempts = 2
        with patch.object(dn, 'simulate', return_value=tree), \
                patch.object(dn, '_is_tree_accepted', return_value=False):
            with self.assertRaises(ec.MaxNFailedAttemptsLimit):
                dn.generate()
        # Resource errors propagate instead of silently filtering the model's realizations.
        with patch.object(dn, 'simulate', side_effect=ec.GenerateFailError('DnQuaSSE', 'limit')):
            with self.assertRaises(ec.GenerateFailError):
                dn.generate()

    # Invalid controls and numerical stalls must not become implicit conditioning.
    # Retain until another layer enforces the same continuous-simulation preconditions.
    def test_failures(self):
        one, zero = LogisticRate(1, 1, 0, 0), LogisticRate(0, 0, 0, 0)
        for kwargs in [{'method': 'unknown'}, {'k': 1.5}, {'diffusion': -1},
                       {'sampling_prob': 2}, {'min_rec_taxa': 3, 'max_rec_taxa': 2},
                       {'drift': [0, 1]}, {'dt_max': 0}, {'dt_max': np.inf},
                       {'max_nodes': 0}, {'max_nodes': 2.5}]:
            with self.assertRaises(ValueError):
                DnQuaSSE(one, zero, 'age', 1, **kwargs)
        for kwargs in [{'sampling_prob': .5}, {'min_rec_taxa': 0}]:
            with self.assertRaises(ValueError):
                DnQuaSSE(one, zero, 'size', 2, **kwargs)
        with self.assertRaises(ec.GenerateFailError):
            DnQuaSSE(zero, zero, 'size', 2).generate()
        underflow = LogisticRate(0, 1, 10000, 1)
        with self.assertRaises(ec.GenerateFailError):
            DnQuaSSE(underflow, zero, 'age', 1).generate()
        with patch('numpy.random.random', return_value=.99):
            with self.assertRaises(ec.GenerateFailError):
                DnQuaSSE(one, zero, 'age', 1, max_steps=1).generate()
        with self.assertRaises(ec.RunTimeLimit):
            DnQuaSSE(one, zero, 'age', 1).simulate(deadline=0)

    # Preserve colors and layout when redrawing or switching DAG values, including CLI/GUI axes.
    # Discrete tests lack colorbars; retire if shared continuous plotting covers these transitions.
    def test_tip_coloring(self):
        from matplotlib import pyplot as plt
        from phylojunction.interface.pjcli.cli_plotting import start_fig_and_axes
        one, zero = LogisticRate(1, 1, 0, 0), LogisticRate(0, 0, 0, 0)
        dn = DnQuaSSE(one, zero, 'size', 2, k=1)
        with patch('numpy.random.random', return_value=0), patch('numpy.random.choice', return_value=0):
            tree = dn.generate()[0]
        tips = list(tree.tree.leaf_node_iter())
        tips[0].trait, tips[1].trait = -2, 3
        # Use only a custom attribute to catch hard-coded names in both plotting views;
        # the simulator tests cover the default. Retire if general rename tests replace this.
        tree.trait = ContinuousTrait(name="body_size")
        for node in tree.tree:
            node.body_size = node.trait
            node.annotations.drop(name="trait")
            node.annotations.add_bound_attribute("body_size")
            del node.trait
        dag = DirectedAcyclicGraph()
        for line in ['numbers <- [1,2,3]', 'rate := quasse_logistic(y0=1,y1=2,midpoint=0,slope=1)']:
            cmdline2dag(dag, line)
        empty = DnQuaSSE(zero, zero, 'age', 0, sampling_prob=0).generate()[0]
        for factory in (plt.subplots, start_fig_and_axes):
            with self.subTest(axes=factory.__name__):
                fig, ax = factory()
                original = ax.get_position().bounds
                colored = None
                try:
                    for reconstructed in (False, True, False, True):
                        tree.plot_node(ax, draw_reconstructed=reconstructed)
                        fig.canvas.draw()
                        if colored is None:
                            colored = ax.get_position().bounds
                        np.testing.assert_allclose(ax.get_position().bounds, colored)
                        self.assertEqual(len(fig.axes), 2)
                        self.assertEqual(ax._pj_trait_colorbar.mappable.get_clim(), (-2, 3))
                        np.testing.assert_array_equal(ax.collections[-1].get_array(), [-2, 3])
                    empty.plot_node(ax, draw_reconstructed=True)
                    self.assertEqual(len(fig.axes), 1)
                    self.assertIsNone(ax._pj_trait_colorbar)
                    np.testing.assert_allclose(ax.get_position().bounds, original)
                    for name in ('numbers', 'rate'):
                        tree.plot_node(ax)
                        dag.name_node_dict[name].plot_node(ax, sample_idx=None)
                        fig.canvas.draw()
                        self.assertEqual(len(fig.axes), 1)
                        self.assertIsNone(ax._pj_trait_colorbar)
                        self.assertIsNone(ax._pj_trait_position)
                        np.testing.assert_allclose(ax.get_position().bounds, original)
                        self.assertFalse(ax.collections)
                        if name == 'numbers':
                            self.assertTrue(ax.patches)
                        else:
                            self.assertFalse(ax.patches or ax.lines or ax.texts)
                        # Blank/histogram redraws with no colorbar also remain valid.
                        dag.name_node_dict[name].plot_node(ax, sample_idx=None)
                    tree.plot_node(ax)
                    self.assertEqual(len(fig.axes), 2)
                    np.testing.assert_allclose(ax.get_position().bounds, colored)
                finally:
                    plt.close(fig)

    # Continuous export must retain all node traits without generating discrete state matrices.
    # Existing discrete export tests cover the other format; retire if a shared trait writer replaces both.
    def test_export(self):
        import io
        import tempfile
        from pathlib import Path
        import pandas as pd
        from phylojunction.readwrite.pj_write import prep_data_df, prep_data_filepaths_dfs
        from phylojunction.inference.revbayes.rb_inference import dag_obj_to_rev_inference_spec
        dag = DirectedAcyclicGraph()
        for line in ['b := quasse_logistic(y0=0,y1=0,midpoint=0,slope=1)',
                     't ~ quasse(n=2,nr=2,birth_rate=b,death_rate=b,stop="age",'
                     'stop_value=[0,2000],start_trait=[1,3],sampling_prob=[1,0])']:
            cmdline2dag(dag, line)
        # Export must read the description, not assume the simulator's default attribute.
        # Existing expected values cover output equivalence; no separate export test is needed.
        for tree in dag.name_node_dict['t'].value:
            tree.trait = ContinuousTrait(name="body_size")
            for node in tree.tree:
                node.body_size = node.trait
                node.annotations.drop(name="trait")
                node.annotations.add_bound_attribute("body_size")
                del node.trait
        filenames, contents = prep_data_filepaths_dfs(*prep_data_df(dag, write_nex_states=True))
        output = dict(zip(filenames, contents))
        self.assertIn('t_annotated_complete.tsv', output)
        self.assertIn('t_reconstructed.tsv', output)
        self.assertFalse(any(name.endswith('.nex') or '_anc_states' in name for name in filenames))
        traits = pd.read_csv(io.StringIO(output['t_traits.tsv']), sep='\t')
        self.assertEqual(list(traits.columns), ['sample', 'replicate', 'node', 'trait', 'alive', 'sampled'])
        self.assertEqual(len(traits), 8)
        np.testing.assert_array_equal(traits['trait'], [1, 1, 1, 1, 3, 3, 3, 3])
        self.assertEqual(output['t_stats.csv'].iloc[0]['origin_age'], 0)
        for tree in dag.name_node_dict['t'].value:
            self.assertFalse(tree.state_count_dict)
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(NotImplementedError):
                dag_obj_to_rev_inference_spec(dag, directory)
            self.assertEqual(list(Path(directory).iterdir()), [])
