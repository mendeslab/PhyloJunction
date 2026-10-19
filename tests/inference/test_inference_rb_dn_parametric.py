import unittest
import re
import math
from statistics import mean, stdev

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.interface.cmdbox.cmd_parse as cmdp
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"

class TestInferenceRevBayesParametricDn(unittest.TestCase):
    # Literal and named scale inputs must export the same rate, retaining dependencies.
    # Moment tests never exercise this translation; remove if rate/scale export is replaced.
    def test_scale_rate_export(self):
        from phylojunction.inference.revbayes.execution import execution_to_rev
        for distribution in ('exponential', 'gamma'):
            for as_rate in (False, True):
                for named in (False, True):
                    for value in (2.0, 4.0):
                        prefix = f's <- {value}\n' if named else ''
                        argument = 's' if named else str(value)
                        parameters = 'rate=' if distribution == 'exponential' else 'shape=3.0,scale='
                        source = prefix + f'x ~ {distribution}({parameters}{argument},rate_parameterization={as_rate})'
                        result = cmdp.compile_model(source).execute(2,123)
                        models, _, _ = execution_to_rev(result,'results')
                        rate = argument if as_rate else f'1.0/({argument})'
                        expected = (f'dnExponential(lambda={rate})' if distribution == 'exponential'
                                    else f'dnGamma(shape=3.0, rate={rate})')
                        for model in models:
                            self.assertIn('x ~ ' + expected, model)

    # Legacy string checks do not exercise distribution aliases, nested IID, or clamping
    # per replicate. These contracts remain required by the shared-graph exporter.
    def test_execution_export(self):
        from phylojunction.inference.revbayes.execution import execution_to_rev
        source = ('a ~ normal(0,1)\nd := normal(a,1)\nx ~ iid(2,iid(2,d),clamp=true)'
                  '\na <- 5\nb ~ normal(a,1)\np ~ beta(2,3)')
        results = cmdp.compile_model(source).execute(2, 42)
        models, logging, _ = execution_to_rev(results, 'results', 10)
        self.assertEqual(len(models), 2)
        for i, model in enumerate(models):
            self.assertIn('a__definition1 ~ dnNormal(mean=0, sd=1)', model)
            self.assertIn('d = dnNormal(mean=a__definition1, sd=1)', model)
            self.assertRegex(model, r'x\[pjInternal\d+\]\[pjInternal\d+\] ~ d')
            self.assertIn(f'x[2][2].clamp({results.value("x", i)[1][1]})', model)
            self.assertIn('b ~ dnNormal(mean=a, sd=1)', model)
            self.assertIn('dnBeta(alpha=2, beta=3)', model)
            self.assertIn('model(a__definition1, x[1][1], x[1][2], x[2][1], x[2][2], b, p)', logging[i])
        for source in ('x ~ iid(0,normal(0,1))', 'b := quasse_constant(1)\nx ~ quasse(b,b,"age",1)'):
            with self.assertRaises(ValueError):
                execution_to_rev(cmdp.compile_model(source).execute(), 'results')

    # Real Rev execution catches invalid collection types and lost dependencies that
    # string checks miss. Keep this integration check while translating distribution values.
    def test_collection_export_in_rev(self):
        import shutil
        import subprocess
        import tempfile
        from pathlib import Path
        from phylojunction.inference.revbayes.execution import execution_to_rev
        binary = shutil.which('rb')
        if binary is None:
            self.skipTest('RevBayes executable is not installed')
        source = ('pjInternal1 ~ normal(0,1)\nrows := [[pjInternal1],[pjInternal1,2]]\n'
                  'ds := [[normal(mu,1) for mu in row] for row in rows]\nxs ~ ds\n'
                  'ys ~ [[normal(pjInternal1,1)]]\n'
                  'mus ~ iid(2,normal(0,1))\n'
                  'zs ~ [normal(mu,sd) for (mu,sd) in zip_strict(mus,[1,2])]\n'
                  'empty ~ iid(0,normal(0,1),clamp=true)\n'
                  'nested ~ [normal([0,1],1)]\n'
                  'indexed ~ [normal(pair[1],pair[2]) for pair in zip_strict(mus,[1,2])]')
        result = cmdp.compile_model(source).execute(1,42)
        with tempfile.TemporaryDirectory() as directory:
            models, logging, paths = execution_to_rev(result,directory,5)
            self.assertIn('ys[1][1] ~ dnNormal(mean=pjInternal1, sd=1)', models[0])
            self.assertIn('pjInternalX', models[0])
            Path(paths[2]).mkdir()
            script = Path(directory) / 'model.rev'
            script.write_text(models[0] + '\nprint(model(pjInternal1,xs[1][1],xs[2][2]))\n' + logging[0])
            completed = subprocess.run([binary, '-q', str(script)], text=True,
                                       capture_output=True, timeout=30)
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertIn('_parents      = [ pjInternal1,', completed.stdout)
        for source in ('x ~ normal([],1)', 'x <- []'):
            with self.assertRaises(ValueError):
                execution_to_rev(cmdp.compile_model(source).execute(), 'unused')

    def test_pj2rb_uniform(self):
        """
        Test if PJ -> .Rev conversion functionality is correct for
        uniform distribution, by comparing against expected strings
        and expected sampled values (from running RevBayes)
        """
        
        dag_obj = pgm.DirectedAcyclicGraph().execute()
        
        n_repl = 100000
        cmd_line = "u ~ iid(" + str(n_repl) + ",unif(min=-1.0,max=1.0))"

        
        dag_obj = dag_obj.append(cmd_line)
        a_node_dag = dag_obj.value("u")

        rb_sample_mean = -0.002172273 # mean(u) in RB
        rb_sample_sd = 0.5781801 # stdev(u) in RB
        rb_sample_min = -0.9999855 # min(u) in RB
        rb_sample_max = 0.9999579 # max(u) in RB

        pj_unif_values = a_node_dag # 100000 floats in list

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
        
        dag_obj = pgm.DirectedAcyclicGraph().execute()
        
        n_repl = 100000
        cmd_line1 = "e1 ~ iid(" + str(n_repl) \
            + ",exponential(rate=0.5,rate_parameterization=\"true\"))" # default is true

        
        dag_obj = dag_obj.append(cmd_line1)
        a_node_dag1 = dag_obj.value("e1")

        # rev_str_list = a_node_dag1.sampling_dn.get_rev_inference_spec_info()
        # str_to_run_in_rb_for_unittest1 = "for (i in 1:" + str(n_repl) + ") {\n" + \
        #     "    e1[i] ~ " + rev_str_list[0] + "\n" + \
        #     "}"
        # print(str_to_run_in_rb_for_unittest1)

        rb_sample_mean1 = 2.013223 # mean(e1) in RB
        rb_sample_sd1 = 2.017317 # stdev(e1) in RB

        pj_exponential_values1 = a_node_dag1 # 100000 floats in list

        cmd_line2 = "e2 ~ iid(" + str(n_repl) \
            + ",exponential(rate=0.5,rate_parameterization=\"false\"))"

        
        dag_obj = dag_obj.append(cmd_line2)
        a_node_dag2 = dag_obj.value("e2")

        # rev_str_list = a_node_dag2.sampling_dn.get_rev_inference_spec_info()
        # str_to_run_in_rb_for_unittest2 = "for (i in 1:" + str(n_repl) + ") {\n" + \
        #     "    e2[i] ~ " + rev_str_list[0] + "\n" + \
        #     "}"
        # print(str_to_run_in_rb_for_unittest2)

        rb_sample_mean2 = 0.50126 # mean(e2) in RB
        rb_sample_sd2 = 0.5042986 # stdev(e2) in RB

        pj_exponential_values2 = a_node_dag2 # 100000 floats in list

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
        
        dag_obj = pgm.DirectedAcyclicGraph().execute()
        
        n_repl = 100000
        cmd_line = "n ~ iid(" + str(n_repl) + ",normal(mean=0.5,sd=0.1))"

        
        dag_obj = dag_obj.append(cmd_line)
        a_node_dag = dag_obj.value("n")

        rb_sample_mean = 0.4998006 # mean(n) in RB
        rb_sample_sd = 0.09983804 # stdev(n) in RB

        pj_normal_values = a_node_dag # 1000 floats in list

        self.assertAlmostEqual(rb_sample_mean, mean(pj_normal_values), delta=1e-2)
        self.assertAlmostEqual(rb_sample_sd, stdev(pj_normal_values), delta=1e-2)


    # Export must preserve both parameterization and parent dependencies; sample-moment
    # tests do not inspect Rev code. Remove if this distribution exporter is replaced.
    def test_pj2rb_lognormal(self):
        from phylojunction.inference.revbayes.execution import execution_to_rev
        for log_space, means in [(True, [0.0,1.0]), (False, [1.0,math.e])]:
            for named in (False,True):
                for mean_value, sd_value in zip(means,[0.5,1.0]):
                    prefix = f'm <- {mean_value}\ns <- {sd_value}\n' if named else ''
                    mean_arg, sd_arg = ('m','s') if named else (str(mean_value),str(sd_value))
                    source = prefix + f'x ~ lognormal({mean_arg},{sd_arg},log_space={log_space})'
                    result = cmdp.compile_model(source).execute(2,123)
                    models, _, _ = execution_to_rev(result,'results')
                    mean = mean_arg if log_space else f'ln({mean_arg})'
                    for model in models:
                        self.assertIn(f'x ~ dnLognormal(mean={mean}, sd={sd_arg})', model)

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