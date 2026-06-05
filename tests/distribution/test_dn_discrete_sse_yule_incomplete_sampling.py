import unittest
import statistics
import math

# pj imports
import phylojunction.utility.helper_functions as pjh
import phylojunction.calculation.math_utils as pjmath
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.distribution.dn_discrete_sse as distsse

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"

class TestYuleIncompleteSamplingTrees(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # not state-dependent (just state 0, and no transition)
        total_n_states = 1

        birth_rate = [
            sseobj.DiscreteStateDependentRate(name="lambda", val=1.0, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0])
        ]

        sampling_prob = [
            sseobj.DiscreteStateDependentProbability(name="rho", val=0.5, state=0)
        ]

        matrix_atomic_rate_params = [ birth_rate ]
        
        matrix_state_dep_probs = [ sampling_prob ]

        state_dep_rate_manager = sseobj.DiscreteStateDependentParameterManager(
            matrix_atomic_rate_params,
            total_n_states
        )

        state_dep_prob_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                matrix_state_dep_probs, total_n_states
            )

        event_handler = sseobj.MacroevolEventHandler(state_dep_rate_manager)
        
        prob_handler = \
            sseobj.DiscreteStateDependentProbabilityHandler(
                state_dep_prob_manager
            )

        cls.sse_stash_complete = sseobj.SSEStash(event_handler)
        cls.sse_stash_incomplete = sseobj.SSEStash(event_handler, prob_handler)


    def test_yule_terminal_counts_complete_sampling(self):
        """
        Test if when no sampling probability information is provided,
        if they are automatically assigned 1.0, and then the number of
        sampled terminal taxa is the same as the number of terminal taxa
        """

        # setting up stopping conditions
        stop_condition = "age"
        stop_condition_value = [2.25] ## 2.25 time units
        start_at_origin = True

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]
    
        # simulations
        sim_batches = list()
        for i in range(n_batches):
            sse_sim = distsse.DnSSE(
                self.sse_stash_complete,
                n=n_sim,
                stop=stop_condition,
                stop_value=stop_condition_value,
                origin=start_at_origin,
                start_states_list=start_states_list,
                epsilon=1e-8,
                runtime_limit=3000,
                debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)

        batch_avgs_extant_terminals = list()
        batch_avgs_extant_sampled_terminals = list()
        for i, batch in enumerate(sim_batches):
            n_extant_terminal_nodes = [ann_tr.n_extant_terminal_nodes for ann_tr in batch]
            avg_extant_terminal = statistics.mean(n_extant_terminal_nodes)
            batch_avgs_extant_terminals.append(avg_extant_terminal)

            n_extant_sampled_terminal_nodes = [ann_tr.n_extant_sampled_terminal_nodes for ann_tr in batch]
            avg_extant_sampled_terminal = statistics.mean(n_extant_sampled_terminal_nodes)
            batch_avgs_extant_sampled_terminals.append(avg_extant_sampled_terminal)

        # if no information on sampling probability is provided
        # then both counts must be the same
        self.assertEqual(
            statistics.mean(batch_avgs_extant_terminals),
            statistics.mean(batch_avgs_extant_sampled_terminals)
        )


    def test_yule_terminal_counts_incomplete_sampling(self):
        """
        Test if sampling probability of 0.5 leads to Yule trees
        with half the expected size when compared to a sampling
        probability of 1.0
        """

        # setting up stopping conditions
        stop_condition = "age"
        stop_condition_value = [2.25] ## 2.25 time units
        start_at_origin = True

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]
    
        # simulations
        sim_batches = list()
        for i in range(n_batches):
            sse_sim = distsse.DnSSE(
                self.sse_stash_incomplete,
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
        #
        # dividing by two because rho = 0.5
        exp_sampled_size = pjmath.exp_extant_count_bd(1.0, 0.0, 2.25) / 2.0 

        batch_avgs_extant_terminals = list()
        batch_avgs_extant_sampled_terminals = list()
        batch_cis = list()
        in_ci_count = 0
        for i, batch in enumerate(sim_batches):
            # n_extant_terminal_nodes = [ann_tr.n_extant_terminal_nodes for ann_tr in batch]
            # avg_extant_terminal = statistics.mean(n_extant_terminal_nodes)
            # batch_avgs_extant_terminals.append(avg_extant_terminal)

            n_extant_sampled_terminal_nodes = [ann_tr.n_extant_sampled_terminal_nodes for ann_tr in batch]
            avg_extant_sampled_terminal = statistics.mean(n_extant_sampled_terminal_nodes)
            batch_avgs_extant_sampled_terminals.append(avg_extant_sampled_terminal)

            stdevs = statistics.stdev(n_extant_sampled_terminal_nodes)
            sterr = stdevs / math.sqrt(n_sim)
            diff = 1.96 * sterr
            batch_cis = (avg_extant_sampled_terminal - diff, avg_extant_sampled_terminal + diff)

            if exp_sampled_size >= batch_cis[0] and exp_sampled_size <= batch_cis[1]:
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
    # $ python3 tests/distribution/test_dn_discrete_sse_yule_incomplete_sampling.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_dn_discrete_sse_yule_incomplete_sampling
    #
    # or 
    #
    # $ python3 -m unittest tests.distribution.test_dn_discrete_sse_yule_incomplete_sampling.TestYuleIncompleteSamplingTrees.test_yule_terminal_counts_incomplete_sampling

    unittest.main()