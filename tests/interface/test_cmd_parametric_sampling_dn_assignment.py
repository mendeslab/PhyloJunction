import unittest
import re
import math
from statistics import mean, stdev

# pj imports #
import phylojunction.pgm.pgm as pgm
import phylojunction.interface.cmdbox.cmd_parse as cmdp
import phylojunction.interface.cmdbox.cmd_parse_utils as cmdu
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"

class TestParametricSamplingDnAssignment(unittest.TestCase):
    
    def test_sampling_unif_assignment(self):
        """
        Test if uniform sampling distribution assignments are correctly
        evaluated and result in the right probabilistic graphical model
        """
        
        dag_obj = pgm.DirectedAcyclicGraph()
        
        cmd_line = "u ~ unif(n=100000, min=-1.0, max=1.0)"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line)
        cmdp.parse_samp_dn_assignment(dag_obj, stoch_node_name, stoch_node_spec, cmd_line)
        a_node_dag = dag_obj.get_node_dag_by_name("u")
  
        self.assertTrue(isinstance(a_node_dag.value, list))
        self.assertEqual(len(a_node_dag.value), 100000)
        self.assertAlmostEqual(0.0, mean(a_node_dag.value), delta=1e-2)
        self.assertLessEqual(-1.0, min(a_node_dag.value))
        self.assertGreater(1.0, max(a_node_dag.value))
        self.assertEqual(1, dag_obj.n_nodes)

    def test_sampling_unif_vectorized_assignment(self):
        """
        Test if uniform (with vectorized inputs) sampling distribution
        assignments are correctly evaluated and result in the right
        probabilistic graphical model
        """

        mins = [0.0, 1.0, 2.0, 3.0, 4.0]
        maxs = [1.0, 2.0, 3.0, 4.0, 5.0]
        tups = tuple((i, maxs[idx])for idx, i in enumerate(mins))

        dag_obj = pgm.DirectedAcyclicGraph()
        
        cmd_line1 = "mins <- [0.0, 1.0, 2.0, 3.0, 4.0]"
        
        stoch_node_name, _, stoch_node_spec = re.split(cmdu.assign_regex, cmd_line1)
        cmdp.parse_variable_assignment(dag_obj, stoch_node_name, stoch_node_spec, cmd_line1)
        
        cmd_line2 = "u ~ unif(n=5, nr=2, min=mins, max=[1.0, 2.0, 3.0, 4.0, 5.0])"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line2)
        cmdp.parse_samp_dn_assignment(dag_obj, stoch_node_name, stoch_node_spec, cmd_line2)
        a_node_dag = dag_obj.get_node_dag_by_name("u")
  
        self.assertTrue(isinstance(a_node_dag.value, list))
        for idx, tup in enumerate(tups):
            self.assertTrue(tup[0] <= a_node_dag.value[idx * 2] < tup[1])
            self.assertTrue(tup[0] <= a_node_dag.value[idx * 2 + 1] < tup[1])
        self.assertEqual(2, dag_obj.n_nodes)

    def test_sampling_exp_assignment(self):
        """
        Test if exponential sampling distribution assignments are
        correctly evaluated and result in the right probabilistic
        graphical model
        """

        ################################################
        # Exponential, rate parameterization (default) #
        ################################################
        dag_obj = pgm.DirectedAcyclicGraph()

        cmd_line1 = "e1 ~ exponential(n=100000, nr=1, rate=0.5)"
        
        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.sampled_as_regex, cmd_line1)
        
        cmdp.parse_samp_dn_assignment(dag_obj,
                                      stoch_node_name,
                                      stoch_node_spec,
                                      cmd_line1)
        
        a_node_dag1 = dag_obj.get_node_dag_by_name("e1")

        self.assertTrue(isinstance(a_node_dag1.value, list))
        self.assertEqual(len(a_node_dag1.value), 100000)
        self.assertAlmostEqual(2.0, mean(a_node_dag1.value), delta=0.05)
        self.assertEqual(1, dag_obj.n_nodes)

        #######################################
        # Exponential, scale parameterization #
        #######################################
        cmd_line2 = ('e2 ~ exponential(n=100000, nr=1, rate=0.5, '
                     'rate_parameterization=\"false\")')

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.sampled_as_regex, cmd_line2)
        cmdp.parse_samp_dn_assignment(dag_obj,
                                      stoch_node_name,
                                      stoch_node_spec,
                                      cmd_line2)
        a_node_dag2 = dag_obj.get_node_dag_by_name("e2")

        self.assertTrue(isinstance(a_node_dag2.value, list))
        self.assertEqual(len(a_node_dag2.value), 100000)
        self.assertAlmostEqual(0.5, mean(a_node_dag2.value), delta=0.05)
        self.assertEqual(2, dag_obj.n_nodes)

    def test_sampling_gamma_assignment(self):
        """
        Test if gamma sampling distribution assignments are correctly
        evaluated and result in the right probabilistic graphical model
        """

        ###########################################
        # Gamma, scale parameterization (default) #
        ###########################################
        dag_obj = pgm.DirectedAcyclicGraph()

        cmd_line1 = "g1 ~ gamma(n=100000, nr=1, shape=0.5, scale=0.5)"
        
        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.sampled_as_regex, cmd_line1)
        
        cmdp.parse_samp_dn_assignment(dag_obj,
                                      stoch_node_name,
                                      stoch_node_spec,
                                      cmd_line1)
        
        a_node_dag1 = dag_obj.get_node_dag_by_name("g1")
        
        self.assertTrue(isinstance(a_node_dag1.value, list))
        self.assertEqual(len(a_node_dag1.value), 100000)
        self.assertAlmostEqual(0.25, mean(a_node_dag1.value), delta=0.05)
        self.assertEqual(1, dag_obj.n_nodes)

        ################################
        # Gamma, rate parameterization #
        ################################
        cmd_line2 = "g2 ~ gamma(n=100000, nr=1, shape=0.5, scale=0.5, " \
            + "rate_parameterization=\"true\")"

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.sampled_as_regex, cmd_line2)

        cmdp.parse_samp_dn_assignment(dag_obj,
                                      stoch_node_name,
                                      stoch_node_spec,
                                      cmd_line2)
        
        a_node_dag2 = dag_obj.get_node_dag_by_name("g2")

        self.assertTrue(isinstance(a_node_dag2.value, list))
        self.assertEqual(len(a_node_dag2.value), 100000)
        self.assertAlmostEqual(1.0, mean(a_node_dag2.value), delta=0.05)
        self.assertEqual(2, dag_obj.n_nodes)

    def test_sampling_normal_assignment(self):
        """
        Test if normal sampling distribution assignments are correctly
        evaluated and result in the right probabilistic graphical
        model
        """
        
        dag_obj = pgm.DirectedAcyclicGraph()
        
        cmd_line1 = "n ~ normal(n=100000, mean=0.5, sd=0.1)"

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.sampled_as_regex, cmd_line1)
        
        cmdp.parse_samp_dn_assignment(dag_obj,
                                      stoch_node_name,
                                      stoch_node_spec,
                                      cmd_line1)
        
        a_node_dag = dag_obj.get_node_dag_by_name("n")
        
        self.assertTrue(isinstance(a_node_dag.value, list))
        self.assertEqual(len(a_node_dag.value), 100000)
        self.assertAlmostEqual(0.5, mean(a_node_dag.value), delta=0.1)
        self.assertEqual(1, dag_obj.n_nodes)

    def test_sampling_ln_assignment(self):
        """
        Test if log-normal sampling distribution assignments are
        correctly evaluated and result in the right probabilistic
        graphical model
        """
        
        ###################################
        # Log-normal, log space (default) #
        ###################################
        dag_obj = pgm.DirectedAcyclicGraph()
        
        cmd_line1 = "ln1 ~ lognormal(n=100000, meanlog=-3.25, sdlog=0.25)"

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.sampled_as_regex, cmd_line1)
        
        cmdp.parse_samp_dn_assignment(dag_obj,
                                      stoch_node_name,
                                      stoch_node_spec,
                                      cmd_line1)
        
        a_node_dag = dag_obj.get_node_dag_by_name("ln1")
        
        self.assertTrue(isinstance(a_node_dag.value, list))
        self.assertEqual(len(a_node_dag.value), 100000)
        self.assertAlmostEqual(0.04, mean(a_node_dag.value), delta=0.1)
        self.assertAlmostEqual(0.01, stdev(a_node_dag.value), delta=0.1)
        self.assertEqual(1, dag_obj.n_nodes)
        
        ##########################
        # Log-normal, real space #
        ##########################

        exp_mean = str(math.exp(-3.25))
        cmd_line2 = "ln2 ~ lognormal(n=100000, meanlog=" + exp_mean \
            + ", sdlog=0.25, log_space=\"false\")"

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.sampled_as_regex, cmd_line2)
        
        cmdp.parse_samp_dn_assignment(dag_obj,
                                      stoch_node_name,
                                      stoch_node_spec,
                                      cmd_line2)
        
        a_node_dag = dag_obj.get_node_dag_by_name("ln2")

        self.assertTrue(isinstance(a_node_dag.value, list))
        self.assertEqual(len(a_node_dag.value), 100000)
        self.assertAlmostEqual(0.04, mean(a_node_dag.value), delta=0.1)
        self.assertAlmostEqual(0.01, stdev(a_node_dag.value), delta=0.1)
        self.assertEqual(2, dag_obj.n_nodes)

    def test_unif_misspec(self):
        """
        Test uniform sampling distribution assignment throws exception
        if mandatory parameters are missing
        """

        dag_obj = pgm.DirectedAcyclicGraph()

        cmd_line1 = "u ~ unif(n=1, nr=1, max=1.0)"

        with self.assertRaises(ec.ScriptSyntaxError) as exc_outer:
            cmdp.cmdline2dag(dag_obj, cmd_line1)

        self.assertEqual(
            str(exc_outer.exception),
            ("ERROR: u ~ unif(n=1, nr=1, max=1.0)\n\nThe line above "
             "had a syntax problem and could not be tokenized. Parsing "
             "the specification of \'unif\' failed. Parameter \'min\' is "
             "missing."))

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.sampled_as_regex, cmd_line1)

        with self.assertRaises(ec.ParseDnInitFailError) as exc_inner:
            cmdp.parse_samp_dn_assignment(
                dag_obj,
                stoch_node_name,
                stoch_node_spec,
                cmd_line1)

        self.assertEqual(str(exc_inner.exception),
                         ("Parsing the specification of \'unif\' failed. "
                          "Parameter \'min\' is missing."))

        cmd_line2 = "u ~ unif(n=1, rate=1.0)"

        with self.assertRaises(ec.ScriptSyntaxError) as exc_outer:
            cmdp.cmdline2dag(dag_obj, cmd_line2)

        self.assertEqual(
            str(exc_outer.exception),
            ("ERROR: u ~ unif(n=1, rate=1.0)\n\nThe line above had a "
             "syntax problem and could not be tokenized. Parsing the "
             "specification of \'unif\' failed. \'rate\' is not a valid "
             "parameter."))

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.sampled_as_regex, cmd_line2)

        with self.assertRaises(ec.ParseDnInitFailError) as exc:
            cmdp.parse_samp_dn_assignment(
                dag_obj,
                stoch_node_name,
                stoch_node_spec,
                cmd_line2)        

        self.assertEqual(str(exc.exception),
            ("Parsing the specification of \'unif\' failed. "
             "\'rate\' is not a valid parameter."))

        cmd_line3 = "u ~ unif()"

        with self.assertRaises(ec.ScriptSyntaxError) as exc_outer2:
            cmdp.cmdline2dag(dag_obj, cmd_line3)

        self.assertEqual(
            str(exc_outer2.exception),
            ("ERROR: u ~ unif()\n\nThe line above had a syntax problem "
             "and could not be tokenized. Something went wrong during "
             "sampling distribution specification. Could not find either "
             "the name of a distribution (e.g., \'normal\') or its "
             "specification (e.g., \'(mean=0.0, sd=1.0)\'), or both."))

    def test_exp_misspec(self):
        """
        Test exponential sampling distribution assignment throws
        exception if mandatory parameters are missing
        """

        dag_obj = pgm.DirectedAcyclicGraph()

        cmd_line1 = "e ~ exponential(n=1, nr=1)"

        with self.assertRaises(ec.ScriptSyntaxError) as exc_outer:
            cmdp.cmdline2dag(dag_obj, cmd_line1)

        self.assertEqual(
            str(exc_outer.exception),
            ("ERROR: e ~ exponential(n=1, nr=1)\n\nThe line above had a "
             "syntax problem and could not be tokenized. Parsing the "
             "specification of \'exponential\' failed. Parameter \'rate\' is "
             "missing."))

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.sampled_as_regex, cmd_line1)

        with self.assertRaises(ec.ParseDnInitFailError) as exc_inner:
            cmdp.parse_samp_dn_assignment(
                dag_obj,
                stoch_node_name,
                stoch_node_spec,
                cmd_line1)        

        self.assertEqual(str(exc_inner.exception),
            "Parsing the specification of \'exponential\' failed. " \
            + "Parameter \'rate\' is missing.")

        cmd_line2 = "e ~ exponential(n=1, min=-1)"

        with self.assertRaises(ec.ScriptSyntaxError) as exc_outer:
            cmdp.cmdline2dag(dag_obj, cmd_line2)

        self.assertEqual(
            str(exc_outer.exception),
            ("ERROR: e ~ exponential(n=1, min=-1)\n\nThe line above had a "
             "syntax problem and could not be tokenized. Parsing the "
             "specification of \'exponential\' failed. \'min\' is not a valid "
             "parameter."))

        cmd_line3 = "e ~ exponential()"

        with self.assertRaises(ec.ScriptSyntaxError) as exc_outer2:
            cmdp.cmdline2dag(dag_obj, cmd_line3)

        self.assertEqual(
            str(exc_outer2.exception),
            ("ERROR: e ~ exponential()\n\nThe line above had a syntax problem "
             "and could not be tokenized. Something went wrong during "
             "sampling distribution specification. Could not find either "
             "the name of a distribution (e.g., \'normal\') or its "
             "specification (e.g., \'(mean=0.0, sd=1.0)\'), or both."))

    def test_gamma_misspec(self):
        """
        Test gamma sampling distribution assignment throws exception
        if mandatory parameters are missing
        """

        dag_obj = pgm.DirectedAcyclicGraph()

        cmd_line1 = "g ~ gamma(n=1, nr=1, shape=0.5)"

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.sampled_as_regex, cmd_line1)

        with self.assertRaises(ec.ParseDnInitFailError) as exc_inner:
            cmdp.parse_samp_dn_assignment(
                dag_obj, stoch_node_name, stoch_node_spec, cmd_line1)      

        self.assertEqual(str(exc_inner.exception),
            "Parsing the specification of \'gamma\' failed. "
            "Parameter \'scale\' is missing.")

        cmd_line2 = "g ~ gamma(n=1, nr=1, rate=0.5)"

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.sampled_as_regex, cmd_line2)

        with self.assertRaises(ec.ScriptSyntaxError) as exc_outer:
            cmdp.cmdline2dag(dag_obj, cmd_line2)

        self.assertEqual(
            str(exc_outer.exception),
            ("ERROR: g ~ gamma(n=1, nr=1, rate=0.5)\n\nThe line above had a "
             "syntax problem and could not be tokenized. Parsing the "
             "specification of \'gamma\' failed. \'rate\' is not a valid "
             "parameter."))

        cmd_line3 = "g ~ gamma()"

        with self.assertRaises(ec.ScriptSyntaxError) as exc_outer2:
            cmdp.cmdline2dag(dag_obj, cmd_line3)

        self.assertEqual(
            str(exc_outer2.exception),
            ("ERROR: g ~ gamma()\n\nThe line above had a syntax problem "
             "and could not be tokenized. Something went wrong during "
             "sampling distribution specification. Could not find either "
             "the name of a distribution (e.g., \'normal\') or its "
             "specification (e.g., \'(mean=0.0, sd=1.0)\'), or both."))

    def test_normal_misspec(self):
            """
            Test normal sampling distribution assignment throws exception
            if mandatory parameters are missing
            """

            dag_obj = pgm.DirectedAcyclicGraph()

            cmd_line1 = "n ~ normal(n=1, nr=1, sd=0.1)"

            stoch_node_name, _, stoch_node_spec = \
                re.split(cmdu.sampled_as_regex, cmd_line1)

            with self.assertRaises(ec.ParseDnInitFailError) as exc_inner:
                cmdp.parse_samp_dn_assignment(
                    dag_obj, stoch_node_name, stoch_node_spec, cmd_line1)        

            self.assertEqual(str(exc_inner.exception),
                "Parsing the specification of \'normal\' failed. " \
                + "Parameter \'mean\' is missing.")

            cmd_line2 = "n ~ normal(n=1, rate=1.0)"

            with self.assertRaises(ec.ScriptSyntaxError) as exc_outer:
                cmdp.cmdline2dag(dag_obj, cmd_line2)

            self.assertEqual(
                str(exc_outer.exception),
                ("ERROR: n ~ normal(n=1, rate=1.0)\n\nThe line above had a "
                "syntax problem and could not be tokenized. Parsing the "
                "specification of \'normal\' failed. \'rate\' is not a valid "
                "parameter."))

            cmd_line3 = "n ~ normal()"

            with self.assertRaises(ec.ScriptSyntaxError) as exc_outer2:
                cmdp.cmdline2dag(dag_obj, cmd_line3)

            self.assertEqual(
                str(exc_outer2.exception),
                ("ERROR: n ~ normal()\n\nThe line above had a syntax problem "
                "and could not be tokenized. Something went wrong during "
                "sampling distribution specification. Could not find either "
                "the name of a distribution (e.g., \'normal\') or its "
                "specification (e.g., \'(mean=0.0, sd=1.0)\'), or both."))


    def test_ln_misspec(self):
        """
        Test log-normal sampling distribution assignment throws
        exception if mandatory parameters are missing
        """

        dag_obj = pgm.DirectedAcyclicGraph()

        cmd_line1 = "ln ~ lognormal(n=1, nr=1, sdlog=0.25)"

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.sampled_as_regex, cmd_line1)

        with self.assertRaises(ec.ParseDnInitFailError) as exc_inner:
            cmdp.parse_samp_dn_assignment(
                dag_obj, stoch_node_name, stoch_node_spec, cmd_line1)        

        self.assertEqual(str(exc_inner.exception), 
            "Parsing the specification of \'lognormal\' failed. " \
            "Parameter \'meanlog\' is missing.")

        cmd_line2 = "ln ~ lognormal(n=1, rate=1.0)"

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.sampled_as_regex, cmd_line2)

        with self.assertRaises(ec.ScriptSyntaxError) as exc_outer:
            cmdp.cmdline2dag(dag_obj, cmd_line2)        

        self.assertEqual(
                str(exc_outer.exception),
                ("ERROR: ln ~ lognormal(n=1, rate=1.0)\n\nThe line above had a "
                "syntax problem and could not be tokenized. Parsing the "
                "specification of \'lognormal\' failed. \'rate\' is not a valid "
                "parameter."))

        cmd_line3 = "ln ~ lognormal()"

        with self.assertRaises(ec.ScriptSyntaxError) as exc_outer2:
            cmdp.cmdline2dag(dag_obj, cmd_line3)

        self.assertEqual(
            str(exc_outer2.exception),
            ("ERROR: ln ~ lognormal()\n\nThe line above had a syntax problem "
            "and could not be tokenized. Something went wrong during "
            "sampling distribution specification. Could not find either "
            "the name of a distribution (e.g., \'normal\') or its "
            "specification (e.g., \'(mean=0.0, sd=1.0)\'), or both."))


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