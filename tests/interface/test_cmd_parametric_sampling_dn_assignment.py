from phylojunction.utility.exception_classes import ScriptEvaluationError
import unittest
from unittest.mock import patch
import random
import numpy as np
from phylojunction.distribution.dn_parametric import DnNormal, DnCategorical
import re
import math
from statistics import mean, stdev

# pj imports #
import phylojunction.pgm.pgm as pgm
import phylojunction.interface.cmdbox.cmd_parse as cmdp
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"

class TestParametricSamplingDnAssignment(unittest.TestCase):
    
    # Exercise math calls through the script binder, including mapped inputs and source errors.
    # Sampler tests miss these cases; retire if shared expression tests cover the same calls.
    def test_log_exp_expressions(self):
        result = cmdp.compile_model('x ~ normal(0, sd=log(2)/2)\ny := exp(x=[0,x])').execute(2, 42)
        for x, y in zip(result.values('x'), result.values('y')):
            self.assertEqual(y, [1., math.exp(x)])
        for expression in ('log(0)', 'log(-1)', 'exp(true)', 'log("2")', 'exp(1000)'):
            with self.subTest(expression=expression), self.assertRaises(ec.ScriptEvaluationError):
                cmdp.compile_model('y := ' + expression).execute()

    def test_sampling_unif_assignment(self):
        """
        Test if uniform sampling distribution assignments are correctly
        evaluated and result in the right probabilistic graphical model
        """
        
        dag_obj = pgm.DirectedAcyclicGraph().execute()
        
        cmd_line = "u ~ iid(100000,unif(min=-1.0,max=1.0))"

        
        dag_obj = dag_obj.append(cmd_line)
        a_node_dag = dag_obj.value("u")
  
        self.assertTrue(isinstance(a_node_dag, list))
        self.assertEqual(len(a_node_dag), 100000)
        self.assertAlmostEqual(0.0, mean(a_node_dag), delta=1e-2)
        self.assertLessEqual(-1.0, min(a_node_dag))
        self.assertGreater(1.0, max(a_node_dag))
        self.assertEqual(1, dag_obj.model.n_nodes)

    def test_sampling_unif_vectorized_assignment(self):
        """Select parameter vectors explicitly, then draw IID values at each setting."""
        for i in range(1, 6):
            source = (f'i <- {i}\nmins <- [0,1,2,3,4]\n'
                      'u ~ iid(2,unif(min=mins[i],max=[1,2,3,4,5][i]))')
            results = cmdp.compile_model(source).execute()
            values = results.value('u')
            self.assertIsInstance(values, list)
            self.assertEqual(len(values), 2)
            self.assertTrue(all(i - 1 <= value < i for value in values))
            self.assertEqual(results.model.n_nodes, 3)

    def test_sampling_exp_assignment(self):
        """
        Test if exponential sampling distribution assignments are
        correctly evaluated and result in the right probabilistic
        graphical model
        """

        ################################################
        # Exponential, rate parameterization (default) #
        ################################################
        dag_obj = pgm.DirectedAcyclicGraph().execute()

        cmd_line1 = "e1 ~ iid(100000,exponential(rate=0.5))"
        
        
        
        dag_obj = dag_obj.append(cmd_line1)
        
        a_node_dag1 = dag_obj.value("e1")

        self.assertTrue(isinstance(a_node_dag1, list))
        self.assertEqual(len(a_node_dag1), 100000)
        self.assertAlmostEqual(2.0, mean(a_node_dag1), delta=0.05)
        self.assertEqual(1, dag_obj.model.n_nodes)

        #######################################
        # Exponential, scale parameterization #
        #######################################
        cmd_line2 = ('e2 ~ iid(100000,exponential(rate=0.5, '
                     'rate_parameterization=\"false\"))')

        
        dag_obj = dag_obj.append(cmd_line2)
        a_node_dag2 = dag_obj.value("e2")

        self.assertTrue(isinstance(a_node_dag2, list))
        self.assertEqual(len(a_node_dag2), 100000)
        self.assertAlmostEqual(0.5, mean(a_node_dag2), delta=0.05)
        self.assertEqual(2, dag_obj.model.n_nodes)

    def test_sampling_gamma_assignment(self):
        """
        Test if gamma sampling distribution assignments are correctly
        evaluated and result in the right probabilistic graphical model
        """

        ###########################################
        # Gamma, scale parameterization (default) #
        ###########################################
        dag_obj = pgm.DirectedAcyclicGraph().execute()

        cmd_line1 = "g1 ~ iid(100000,gamma(shape=0.5,scale=0.5))"
        
        
        
        dag_obj = dag_obj.append(cmd_line1)
        
        a_node_dag1 = dag_obj.value("g1")
        
        self.assertTrue(isinstance(a_node_dag1, list))
        self.assertEqual(len(a_node_dag1), 100000)
        self.assertAlmostEqual(0.25, mean(a_node_dag1), delta=0.05)
        self.assertEqual(1, dag_obj.model.n_nodes)

        ################################
        # Gamma, rate parameterization #
        ################################
        cmd_line2 = "g2 ~ iid(100000,gamma(shape=0.5,scale=0.5, " \
            + "rate_parameterization=\"true\"))"

        

        dag_obj = dag_obj.append(cmd_line2)
        
        a_node_dag2 = dag_obj.value("g2")

        self.assertTrue(isinstance(a_node_dag2, list))
        self.assertEqual(len(a_node_dag2), 100000)
        self.assertAlmostEqual(1.0, mean(a_node_dag2), delta=0.05)
        self.assertEqual(2, dag_obj.model.n_nodes)

    def test_sampling_normal_assignment(self):
        """
        Test if normal sampling distribution assignments are correctly
        evaluated and result in the right probabilistic graphical
        model
        """
        
        dag_obj = pgm.DirectedAcyclicGraph().execute()
        
        cmd_line1 = "n ~ iid(100000,normal(mean=0.5,sd=0.1))"

        
        
        dag_obj = dag_obj.append(cmd_line1)
        
        a_node_dag = dag_obj.value("n")
        
        self.assertTrue(isinstance(a_node_dag, list))
        self.assertEqual(len(a_node_dag), 100000)
        self.assertAlmostEqual(0.5, mean(a_node_dag), delta=0.1)
        self.assertEqual(1, dag_obj.model.n_nodes)

    def test_sampling_ln_assignment(self):
        """
        Test if log-normal sampling distribution assignments are
        correctly evaluated and result in the right probabilistic
        graphical model
        """
        
        ###################################
        # Log-normal, log space (default) #
        ###################################
        dag_obj = pgm.DirectedAcyclicGraph().execute()
        
        cmd_line1 = "ln1 ~ iid(100000,lognormal(meanlog=-3.25,sdlog=0.25))"

        
        
        dag_obj = dag_obj.append(cmd_line1)
        
        a_node_dag = dag_obj.value("ln1")
        
        self.assertTrue(isinstance(a_node_dag, list))
        self.assertEqual(len(a_node_dag), 100000)
        self.assertAlmostEqual(0.04, mean(a_node_dag), delta=0.1)
        self.assertAlmostEqual(0.01, stdev(a_node_dag), delta=0.1)
        self.assertEqual(1, dag_obj.model.n_nodes)
        
        ##########################
        # Log-normal, real space #
        ##########################

        exp_mean = str(math.exp(-3.25))
        cmd_line2 = "ln2 ~ iid(100000,lognormal(meanlog=" + exp_mean \
            + ", sdlog=0.25, log_space=\"false\"))"

        
        
        dag_obj = dag_obj.append(cmd_line2)
        
        a_node_dag = dag_obj.value("ln2")

        self.assertTrue(isinstance(a_node_dag, list))
        self.assertEqual(len(a_node_dag), 100000)
        self.assertAlmostEqual(0.04, mean(a_node_dag), delta=0.1)
        self.assertAlmostEqual(0.01, stdev(a_node_dag), delta=0.1)
        self.assertEqual(2, dag_obj.model.n_nodes)

    def test_unif_misspec(self):
        """Test uniform sampling distribution assignment throws exception
if mandatory parameters are missing"""
        for source in ['u ~ unif(max=1.0)', 'u ~ unif(rate=1.0)', 'u ~ unif()']:
            with self.subTest(source=source), self.assertRaises(ScriptEvaluationError) as caught:
                cmdp.compile_model(source).execute()
            error = caught.exception
            self.assertEqual((error.source, error.variable, error.execution), (source, source.split()[0], 1))
            self.assertIn('missing a required argument', str(error))
            self.assertIsInstance(error.__cause__, TypeError)

    def test_exp_misspec(self):
        """Test exponential sampling distribution assignment throws
exception if mandatory parameters are missing"""
        for source in ['e ~ exponential()', 'e ~ exponential(min=-1)', 'e ~ exponential()']:
            with self.subTest(source=source), self.assertRaises(ScriptEvaluationError) as caught:
                cmdp.compile_model(source).execute()
            error = caught.exception
            self.assertEqual((error.source, error.variable, error.execution), (source, source.split()[0], 1))
            self.assertIn('missing a required argument', str(error))
            self.assertIsInstance(error.__cause__, TypeError)

    def test_gamma_misspec(self):
        """Test gamma sampling distribution assignment throws exception
if mandatory parameters are missing"""
        for source in ['g ~ gamma(shape=0.5)', 'g ~ gamma(rate=0.5)', 'g ~ gamma()']:
            with self.subTest(source=source), self.assertRaises(ScriptEvaluationError) as caught:
                cmdp.compile_model(source).execute()
            error = caught.exception
            self.assertEqual((error.source, error.variable, error.execution), (source, source.split()[0], 1))
            self.assertIn('missing a required argument', str(error))
            self.assertIsInstance(error.__cause__, TypeError)

    def test_normal_misspec(self):
        """Test normal sampling distribution assignment throws exception
if mandatory parameters are missing"""
        for source in ['n ~ normal(sd=0.1)', 'n ~ normal(rate=1.0)', 'n ~ normal()']:
            with self.subTest(source=source), self.assertRaises(ScriptEvaluationError) as caught:
                cmdp.compile_model(source).execute()
            error = caught.exception
            self.assertEqual((error.source, error.variable, error.execution), (source, source.split()[0], 1))
            self.assertIn('missing a required argument', str(error))
            self.assertIsInstance(error.__cause__, TypeError)


    def test_ln_misspec(self):
        """Test log-normal sampling distribution assignment throws
exception if mandatory parameters are missing"""
        for source in ['ln ~ lognormal(sdlog=0.25)', 'ln ~ lognormal(rate=1.0)', 'ln ~ lognormal()']:
            with self.subTest(source=source), self.assertRaises(ScriptEvaluationError) as caught:
                cmdp.compile_model(source).execute()
            error = caught.exception
            self.assertEqual((error.source, error.variable, error.execution), (source, source.split()[0], 1))
            self.assertIn('missing a required argument', str(error))
            self.assertIsInstance(error.__cause__, TypeError)

    # Whole-model rejection must redraw parents and recompute descendants; distribution-level
    # retry tests cannot cover this. Retire if equivalent execution coverage replaces this group.
    def test_conditions_restart_model(self):
        model = cmdp.compile_model('x ~ normal(0, 1)\ny := 2*x\n'
                                   'condition(y > 0)\ncondition(y < 6)\nz := y + 1')
        completed = []
        with patch.object(DnNormal, 'draw', side_effect=[-1., 4., 2., -2., 1.]) as draw:
            result = model.execute(runs=2, max_attempts=3, on_replicate_complete=completed.append)
        self.assertEqual(result.values('x'), [2., 1.])
        self.assertEqual(result.values('y'), [4., 2.])
        self.assertEqual(result.values('z'), [5., 3.])
        self.assertEqual(draw.call_count, 5)
        self.assertEqual(completed, [1, 2])

    # Execution errors must not become rejection, and attempt limits must count draws exactly.
    # Sampler tests miss this boundary; retire if a shared acceptance layer tests it equivalently.
    def test_condition_limits_and_errors(self):
        source = 'x ~ normal(0, 1)\ncondition(x > 0)'
        with patch.object(DnNormal, 'draw', return_value=-1.) as draw:
            with self.assertRaises(ScriptEvaluationError) as caught:
                cmdp.compile_model(source).execute(max_attempts=3)
        self.assertEqual(draw.call_count, 3)
        error = caught.exception
        self.assertEqual((error.source, error.line_number, error.execution), ('condition(x > 0)', 2, 1))
        self.assertIn('after 3 attempts', str(error))
        for expression in ('1', '[true, false]'):
            with self.subTest(expression=expression), patch.object(DnNormal, 'draw', return_value=1.) as draw:
                with self.assertRaises(ScriptEvaluationError) as caught:
                    cmdp.compile_model('x ~ normal(0, 1)\ncondition(' + expression + ')').execute()
                self.assertEqual(draw.call_count, 1)
                self.assertIsInstance(caught.exception.__cause__, TypeError)
                self.assertIn('Boolean scalar', str(caught.exception))
        failure = RuntimeError('sampling failed')
        with patch.object(DnNormal, 'draw', side_effect=failure) as draw:
            with self.assertRaises(ScriptEvaluationError) as caught:
                cmdp.compile_model(source).execute()
        self.assertEqual(draw.call_count, 1)
        self.assertIs(caught.exception.__cause__, failure)

    # Failed condition append must leave earlier results and definitions untouched, even after
    # preceding assignments. Ordinary append tests miss this; retire if append supports conditions.
    def test_condition_append_is_transactional(self):
        result = cmdp.compile_model('x := 1').execute()
        with self.assertRaisesRegex(ScriptEvaluationError, 'Cannot append conditions'):
            result.append('y := 2\ncondition(x > 0)')
        self.assertEqual(result.value('x'), 1)
        self.assertEqual(set(result.model.name_node_dict), {'x'})
        self.assertEqual(result.model.n_nodes, 1)
        self.assertEqual(len(result.model.statements), 1)

    # Real draws detect reseeding on rejection and failure to restore RNG state on errors, which
    # mocked draws cannot. Retire if execution stops using these global random-number generators.
    def test_condition_random_stream(self):
        saved_python, saved_numpy = random.getstate(), np.random.get_state()
        try:
            for expression, distribution in (
                    ('normal(0, 1)', DnNormal(0, 1)),
                    ('categorical(cats=[1,-1], probs=[1,1])', DnCategorical([1, -1], [1, 1]))):
                with self.subTest(expression=expression):
                    random.seed(1)
                    np.random.seed(1)
                    for attempt in range(10):
                        expected = distribution.draw()
                        if expected < 0:
                            break
                    else:
                        self.fail('Reference loop did not accept within ten attempts')
                    self.assertGreater(attempt, 0)
                    next_numpy, next_python = np.random.random(), random.random()
                    result = cmdp.compile_model('x ~ ' + expression + '\ncondition(x < 0)').execute(
                        random_seed=1, max_attempts=10)
                    self.assertEqual(result.value('x'), expected)
                    self.assertEqual(np.random.random(), next_numpy)
                    self.assertEqual(random.random(), next_python)

            before_python, before_numpy = random.getstate(), np.random.get_state()
            model = cmdp.compile_model('x ~ normal(0, 1)\ncondition(false)')
            # Count real draws to ensure the error follows sampling, rather than failing before it.
            with patch.object(DnNormal, 'draw', autospec=True, side_effect=DnNormal.draw) as draw:
                with self.assertRaises(ScriptEvaluationError) as caught:
                    model.execute(random_seed=2, max_attempts=3)
            self.assertEqual(draw.call_count, 3)
            self.assertEqual((caught.exception.source, caught.exception.execution), ('condition(false)', 1))
            self.assertIn('after 3 attempts', str(caught.exception))
            self.assertEqual(random.getstate(), before_python)
            after_numpy = np.random.get_state()
            self.assertEqual(after_numpy[0], before_numpy[0])
            np.testing.assert_array_equal(after_numpy[1], before_numpy[1])
            self.assertEqual(after_numpy[2:], before_numpy[2:])
        finally:
            random.setstate(saved_python)
            np.random.set_state(saved_numpy)


if __name__ == "__main__":
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
    # on the terminal, remember to add "src/" to
    # PYTHONPATH (system variable), or to set it if it does not
    # exist -- don't forget to export it!
    # 
    # Then you can do:
    # $ python3 tests/interface/test_cmd_parametric_sampling_dn_assignment.py
    # 
    # or
    #
    # $ python3 -m tests.interface.test_cmd_parametric_sampling_dn_assignment
    #
    # or 
    #
    # $ python3 -m unittest tests.interface.test_cmd_parametric_sampling_dn_assignment.TestParametricSamplingDnAssignment.test_sampling_unif_assignment
    
    unittest.main()