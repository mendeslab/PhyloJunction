from phylojunction.utility.exception_classes import ScriptEvaluationError
import unittest
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