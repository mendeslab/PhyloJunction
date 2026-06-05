import unittest
import math
import statistics

# pj imports
import phylojunction.utility.helper_functions as pjh
import phylojunction.calculation.math_utils as pjmath
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.distribution.dn_discrete_sse as distsse

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestYuleTrees(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        # not state-dependent (just state 0, and no transition)
        rates_t0_s0 = [ 
            sseobj.DiscreteStateDependentRate(
                name="lambda", val=1.0,
                event=sseobj.MacroevolEvent.W_SPECIATION,
                states=[0,0,0]
            )
        ]
    
        matrix_atomic_rate_params = [rates_t0_s0]
        
        state_dep_par_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                matrix_atomic_rate_params,
                1)

        cls.event_handler = \
            sseobj.MacroevolEventHandler(state_dep_par_manager)

        cls.sse_stash = sseobj.SSEStash(cls.event_handler)


    def test_expected_root_height_yule(self):
        """
        Test if Yule trees simulated with PJ have the expected root age
        (against the theoretical expectation).
        """

        # setting up stopping conditions
        stop_condition = "size"
        stop_condition_value = [20] ## 20 species
        start_at_origin = True

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]
        # seeds_list = [i+1 for i in range(n_sim)]
    
        # simulations
        print("\n\nRunning TestYuleTrees.test_expected_root_height_yule")
        sim_batches = list()
        for i in range(n_batches):
            # print("Doing batch " + str(n_batches - i))
            sse_sim = distsse.DnSSE(
                self.sse_stash,
                n=n_sim,
                origin=start_at_origin,
                start_states_list=start_states_list,
                stop=stop_condition,
                stop_value=stop_condition_value,
                condition_on_survival=True,
                epsilon=1e-12,
                runtime_limit=3000)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)

        # parsing results
        exp_age = pjmath.exp_root_height_yule_ntaxa(1.0,
                                                    stop_condition_value[0])
        
        extant_count_set = set()
        batch_cis = list()
        in_ci_count = 0
        for i, batch in enumerate(sim_batches):
            extant_count_set = extant_count_set.union(
                set(ann_tr.n_extant_terminal_nodes for ann_tr in batch)
            )
            root_ages = [ann_tr.root_age for ann_tr in batch]
            stdevs = statistics.stdev(root_ages)
            sterr = stdevs / math.sqrt(n_sim)
            diff = 1.96 * sterr
            avg_age = statistics.mean(root_ages)
            batch_cis = (avg_age - diff, avg_age + diff)

            if exp_age >= batch_cis[0] and exp_age <= batch_cis[1]:
                in_ci_count += 1

        self.assertEqual(
            extant_count_set, 
            {20},
            ("All trees should have only 20 extant observable nodes. But "
             "found trees of size: ") +
            ", ".join(str(s) for s in extant_count_set) + "."
        )

        print("\n95% CI includes expectation " + str(in_ci_count) + " times.")
        
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        
        self.assertAlmostEqual(
            in_ci_count,
            exp_count,
            msg="Truth should be contained within 95%-CI " + str(exp_count) \
                + " (+/- " + str(a_delta) + ") out of 100 times.",
            delta=a_delta
        )


    def test_expected_size_yule(self):
        """
        Test if Yule trees simulated with PJ have the expected number
        of extant observable nodes (against the theoretical
        expectation).
        """

        # setting up stopping conditions
        stop_condition = "age"
        stop_condition_value = [2.0] ## 2 time units
        start_at_origin = True

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]
        # seeds_list = [i+1 for i in range(n_sim)]
    
        # simulations
        print("\n\nRunning TestYuleTrees.test_expected_size_yule")
        sim_batches = list()
        for i in range(n_batches):
            sse_sim = distsse.DnSSE(self.sse_stash,
                                    n=n_sim,
                                    origin=start_at_origin,
                                    start_states_list=start_states_list,
                                    stop=stop_condition,
                                    stop_value=stop_condition_value,
                                    epsilon=1e-8,
                                    runtime_limit=3000,
                                    debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)

        # parsing results
        exp_size = pjmath.exp_extant_count_bd(1.0, 0.0, 2.0)
        
        batch_cis = list()
        in_ci_count = 0
        for i, batch in enumerate(sim_batches):
            n_extant_terminal_nodes = [ann_tr.n_extant_terminal_nodes
                                       for ann_tr in batch]
            stdevs = statistics.stdev(n_extant_terminal_nodes)
            sterr = stdevs / math.sqrt(n_sim)
            diff = 1.96 * sterr
            avg_size = statistics.mean(n_extant_terminal_nodes)
            batch_cis = (avg_size - diff, avg_size + diff)

            if exp_size >= batch_cis[0] and exp_size <= batch_cis[1]:
                in_ci_count += 1

        print("\n95% CI includes expectation " + str(in_ci_count) + " times.")
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        
        self.assertAlmostEqual(
            in_ci_count, exp_count,
            msg="Truth should be contained within 95%-CI " \
            + str(exp_count) + " (+/- " + str(a_delta) \
            + ") out of 100 times.", delta=a_delta)


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
    # $ python3 tests/distribution/test_dn_discrete_sse_yule.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_dn_discrete_sse_yule
    #
    # or 
    #
    # $ python3 -m unittest tests.distribution.test_dn_discrete_sse_yule.TestYuleTrees.test_expected_root_height_yule

    unittest.main()