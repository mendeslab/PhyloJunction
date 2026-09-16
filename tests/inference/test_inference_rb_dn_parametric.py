import unittest
import re
import math
from statistics import mean, stdev

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.interface.cmdbox.cmd_parse as cmdp
import phylojunction.interface.cmdbox.cmd_parse_utils as cmdu
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"

class TestInferenceRevBayesParametricDn(unittest.TestCase):
    # Literal and named scale inputs must export the same rate, retaining dependencies.
    # Moment tests never exercise this translation; remove if rate/scale export is replaced.
    def test_scale_rate_export(self):
        for distribution in ["exponential", "gamma"]:
            for as_rate in [False, True]:
                for named in [False, True]:
                    with self.subTest(distribution=distribution, as_rate=as_rate, named=named):
                        prefix = "s <- [2.0, 4.0]\n" if named else ""
                        argument = "s" if named else "[2.0, 4.0]"
                        params = "rate=" if distribution == "exponential" else "shape=[3.0, 3.0], scale="
                        model = prefix + (f'x ~ {distribution}(n=2, {params}{argument}, '
                                          f'rate_parameterization="{str(as_rate).lower()}")')
                        dag = cmdp.script2dag(model, in_pj_file=False, random_seed=123)
                        rates = (["s" if as_rate else "1.0 / s"] * 2 if named else
                                 (["2.0", "4.0"] if as_rate else ["0.5", "0.25"]))
                        template = ("dnExponential(lambda={})" if distribution == "exponential"
                                    else "dnGamma(3.0, rate={})")
                        self.assertEqual(dag.get_node_dag_by_name("x").sampling_dn.get_rev_inference_spec_info(),
                                         [template.format(rate) for rate in rates])

    def test_pj2rb_uniform(self):
        """
        Test if PJ -> .Rev conversion functionality is correct for
        uniform distribution, by comparing against expected strings
        and expected sampled values (from running RevBayes)
        """
        
        dag_obj = pgm.DirectedAcyclicGraph()
        
        n_repl = 100000
        cmd_line = "u ~ unif(n=1, nr=" + str(n_repl) + ", min=-1.0, max=1.0)"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line)
        cmdp.parse_samp_dn_assignment(dag_obj, stoch_node_name, stoch_node_spec, cmd_line)
        a_node_dag = dag_obj.get_node_dag_by_name("u")

        rb_sample_mean = -0.002172273 # mean(u) in RB
        rb_sample_sd = 0.5781801 # stdev(u) in RB
        rb_sample_min = -0.9999855 # min(u) in RB
        rb_sample_max = 0.9999579 # max(u) in RB

        pj_unif_values = a_node_dag.value # 100000 floats in list

        self.assertAlmostEqual(rb_sample_mean, mean(pj_unif_values), delta=1e-2)
        self.assertAlmostEqual(rb_sample_sd, stdev(pj_unif_values), delta=1e-2)
        self.assertAlmostEqual(rb_sample_min, min(pj_unif_values), delta=1e-2)
        self.assertAlmostEqual(rb_sample_max, max(pj_unif_values), delta=1e-2)
        

    def test_pj2rb_exponential(self):
        """
        Test if PJ -> .Rev conversion functionality is correct for
        exponential distribution, by comparing against expected strings
        and expected sampled values (from running RevBayes)
        """
        
        dag_obj = pgm.DirectedAcyclicGraph()
        
        n_repl = 100000
        cmd_line1 = "e1 ~ exponential(n=1, nr=" + str(n_repl) \
            + ", rate=0.5, rate_parameterization=\"true\")" # default is true

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.sampled_as_regex, cmd_line1)
        cmdp.parse_samp_dn_assignment(dag_obj,
                                      stoch_node_name,
                                      stoch_node_spec,
                                      cmd_line1)
        a_node_dag1 = dag_obj.get_node_dag_by_name("e1")

        # rev_str_list = a_node_dag1.sampling_dn.get_rev_inference_spec_info()
        # str_to_run_in_rb_for_unittest1 = "for (i in 1:" + str(n_repl) + ") {\n" + \
        #     "    e1[i] ~ " + rev_str_list[0] + "\n" + \
        #     "}"
        # print(str_to_run_in_rb_for_unittest1)

        rb_sample_mean1 = 2.013223 # mean(e1) in RB
        rb_sample_sd1 = 2.017317 # stdev(e1) in RB

        pj_exponential_values1 = a_node_dag1.value # 100000 floats in list

        cmd_line2 = "e2 ~ exponential(n=1, nr=" + str(n_repl) \
            + ", rate=0.5, rate_parameterization=\"false\")"

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.sampled_as_regex, cmd_line2)
        cmdp.parse_samp_dn_assignment(dag_obj,
                                      stoch_node_name,
                                      stoch_node_spec,
                                      cmd_line2)
        a_node_dag2 = dag_obj.get_node_dag_by_name("e2")

        # rev_str_list = a_node_dag2.sampling_dn.get_rev_inference_spec_info()
        # str_to_run_in_rb_for_unittest2 = "for (i in 1:" + str(n_repl) + ") {\n" + \
        #     "    e2[i] ~ " + rev_str_list[0] + "\n" + \
        #     "}"
        # print(str_to_run_in_rb_for_unittest2)

        rb_sample_mean2 = 0.50126 # mean(e2) in RB
        rb_sample_sd2 = 0.5042986 # stdev(e2) in RB

        pj_exponential_values2 = a_node_dag2.value # 100000 floats in list

        self.assertAlmostEqual(rb_sample_mean1,
                               mean(pj_exponential_values1),
                               delta=1e-1)
        self.assertAlmostEqual(rb_sample_sd1,
                               stdev(pj_exponential_values1),
                               delta=1e-1)
        self.assertAlmostEqual(rb_sample_mean2,
                               mean(pj_exponential_values2),
                               delta=1e-1)
        self.assertAlmostEqual(rb_sample_sd2,
                               stdev(pj_exponential_values2),
                               delta=1e-1)


    def test_pj2rb_gamma(self):
        """
        Test if PJ -> .Rev conversion functionality is correct for
        gamma distribution, by comparing against expected strings
        and expected sampled values (from running RevBayes)
        """
        pass


    def test_pj2rb_normal(self):
        """
        Test if PJ -> .Rev conversion functionality is correct for
        normal distribution, by comparing against expected strings
        and expected sampled values (from running RevBayes)
        """
        
        dag_obj = pgm.DirectedAcyclicGraph()
        
        n_repl = 100000
        cmd_line = "n ~ normal(n=1, nr=" + str(n_repl) + ", mean=0.5, sd=0.1)"

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.sampled_as_regex, cmd_line)
        cmdp.parse_samp_dn_assignment(dag_obj,
                                      stoch_node_name,
                                      stoch_node_spec,
                                      cmd_line)
        a_node_dag = dag_obj.get_node_dag_by_name("n")

        rb_sample_mean = 0.4998006 # mean(n) in RB
        rb_sample_sd = 0.09983804 # stdev(n) in RB

        pj_normal_values = a_node_dag.value # 1000 floats in list

        self.assertAlmostEqual(rb_sample_mean, mean(pj_normal_values), delta=1e-2)
        self.assertAlmostEqual(rb_sample_sd, stdev(pj_normal_values), delta=1e-2)


    # Export must preserve both parameterization and parent dependencies; sample-moment
    # tests do not inspect Rev code. Remove if this distribution exporter is replaced.
    def test_pj2rb_lognormal(self):
        for log_space, means in [(True, [0.0, 1.0]), (False, [1.0, math.e])]:
            for named in [False, True]:
                with self.subTest(log_space=log_space, named=named):
                    prefix = "m <- " + str(means) + "\ns <- [0.5, 1.0]\n" if named else ""
                    mean_arg = "m" if named else str(means)
                    sd_arg = "s" if named else "[0.5, 1.0]"
                    model = prefix + ('x ~ lognormal(n=2, meanlog=' + mean_arg + ', sdlog=' + sd_arg
                                      + ', log_space="' + str(log_space).lower() + '")')
                    dag = cmdp.script2dag(model, in_pj_file=False, random_seed=123)
                    dn = dag.get_node_dag_by_name("x").sampling_dn
                    means_out = (["m" if log_space else "ln(m)"] * 2 if named else ["0.0", "1.0"])
                    sds_out = ["s", "s"] if named else ["0.5", "1.0"]
                    self.assertEqual(dn.get_rev_inference_spec_info(),
                                     [f"dnLognormal(mean={m}, sd={s})" for m, s in zip(means_out, sds_out)])

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
    # $ python3 tests/inference/test_inference_rb_dn_parametric.py
    # 
    # or
    #
    # $ python3 -m tests.inference.test_inference_rb_dn_parametric
    #
    # or 
    #
    # $ python3 -m unittest tests.inference.test_inference_rb_dn_parametric.TestInferenceRevBayesParametricDn.test_pj2rb_uniform

    unittest.main()