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


class TestBiSSEIncompleteSamplingTrees(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        total_n_states = 2

        rates_s0 = [
            sseobj.DiscreteStateDependentRate(
                name="lambda0",
                val=1.0,
                event=sseobj.MacroevolEvent.W_SPECIATION,
                states=[0,0,0]),
            sseobj.DiscreteStateDependentRate(
                name="mu0",
                val=0.5,
                event=sseobj.MacroevolEvent.EXTINCTION,
                states=[0]),
            sseobj.DiscreteStateDependentRate(
                name="q01",
                val=0.25,
                event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION,
                states=[0,1])
        ]

        rates_s1 = [
            sseobj.DiscreteStateDependentRate(
                name="lambda1",
                val=1.0,
                event=sseobj.MacroevolEvent.W_SPECIATION,
                states=[1,1,1]),
            sseobj.DiscreteStateDependentRate(
                name="mu1",
                val=0.5,
                event=sseobj.MacroevolEvent.EXTINCTION,
                states=[1]),
            sseobj.DiscreteStateDependentRate(
                name="q10",
                val=0.25,
                event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION,
                states=[1,0])
        ]

        rates = rates_s0 + rates_s1
        
        sampling_probs = [
            sseobj.DiscreteStateDependentProbability(
                name="rho",
                val=0.25,
                state=0),
            sseobj.DiscreteStateDependentProbability(
                name="rho", 
                val=0.75,
                state=1)
        ]

        matrix_atomic_rates = [rates]
        
        matrix_state_dep_probs = [sampling_probs]

        state_dep_rate_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                matrix_atomic_rates,
                total_n_states)

        state_dep_prob_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                matrix_state_dep_probs,
                total_n_states)

        event_handler = \
            sseobj.MacroevolEventHandler(state_dep_rate_manager)
        
        prob_handler = \
            sseobj.DiscreteStateDependentProbabilityHandler(
                state_dep_prob_manager
            )

        cls.sse_stash_complete = sseobj.SSEStash(event_handler)
        cls.sse_stash_incomplete = \
            sseobj.SSEStash(event_handler, prob_handler)


    def test_bisse_sampled_terminal_counts_incomplete_sampling(self):
        """
        Test if BiSSE trees simulated with different sampling
        probabilities for each state (using PJ's code made specifically
        for that) match BiSSE trees with perfect sampling but whose
        taxa are masked off "as if" they had been sampled differentially
        (this latter simulation is the analog of a ground truth
        provided by a different simulator in R, for example). Only
        PJ code is used in this test.
        """

        # setting up stopping conditions
        stop_condition = "age"
        stop_condition_value = [2.75]  # 2.75 time units
        start_at_origin = True

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        # start_states_list = [0 for i in range(n_sim)]
        start_states_list = [0, 1] * 50
    
        # we will simulate batches of trees with perfect sampling, but then
        # artificially "mask"-off some of the extant taxa in both states,
        # proportionally to the sampling probabilities we will use in the
        # batch of simulations below
        # 
        # these batches here could also have come from a different simulator
        # in R, for example; it is the same logic used in the other unit
        # tests -- they will serve as the expectation for the model
        # implementation (non-zero sampling probabilities) we are testing
        print(("\n\nRunning TestBiSSEIncompleteSamplingTrees.test_bisse_"
               "sampled_terminal_counts_incomplete_sampling "
                "(first 100 batches)"))
        sim_batches_complete = list()
        for i in range(n_batches):
            sse_sim = distsse.DnSSE(
                self.sse_stash_complete,
                n=n_sim,
                origin=start_at_origin,
                start_states_list=start_states_list,
                stop=stop_condition,
                stop_value=stop_condition_value,
                epsilon=1e-8,
                runtime_limit=3000,
                debug=False)

            trs = sse_sim.generate()

            sim_batches_complete.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)

        batch_avgs_extant_sampled_terminals_complete_n0 = list()
        batch_avgs_extant_sampled_terminals_complete_n1 = list()
        n0_ci_width_sampled_complete_list = list()
        n1_ci_width_sampled_complete_list = list()
        for i, batch in enumerate(sim_batches_complete):
            approx_n0s_sampled_complete = \
                [ann_tr.extant_terminal_state_count_dict[0] * 0.25 for ann_tr in batch]
            approx_n1s_sampled_complete = \
                [ann_tr.extant_terminal_state_count_dict[1] * 0.75 for ann_tr in batch]
            
            mean_n0_sampled_complete = \
                statistics.mean(approx_n0s_sampled_complete)
            mean_n1_sampled_complete = \
                statistics.mean(approx_n1s_sampled_complete)
            
            batch_avgs_extant_sampled_terminals_complete_n0.\
                append(mean_n0_sampled_complete)
            batch_avgs_extant_sampled_terminals_complete_n1.\
                append(mean_n1_sampled_complete)

            stdevs_n0_sampled_complete = \
                statistics.stdev(approx_n0s_sampled_complete)
            stdevs_n1_sampled_complete = \
                statistics.stdev(approx_n1s_sampled_complete)
            
            sterr_n0_sampled_complete = \
                stdevs_n0_sampled_complete / math.sqrt(n_sim)
            sterr_n1_sampled_complete = \
                stdevs_n1_sampled_complete / math.sqrt(n_sim)

            n0_ci_width_sampled_complete = 1.96 * sterr_n0_sampled_complete
            n0_ci_width_sampled_complete_list.append(
                n0_ci_width_sampled_complete)
            
            n1_ci_width_sampled_complete = 1.96 * sterr_n1_sampled_complete
            n1_ci_width_sampled_complete_list.append(
                n1_ci_width_sampled_complete)
        
        # these are the "expected" taxon counts (100 values, one per batch)
        # against which we will compare our simulations with non-zero
        # sampling probabilities below
        expected_avg_sampled_incomplete_n0 = \
            statistics.mean(batch_avgs_extant_sampled_terminals_complete_n0)
        expected_avg_sampled_incomplete_n1 = \
            statistics.mean(batch_avgs_extant_sampled_terminals_complete_n1)

        # now we simulate with the model we are actually testing
        print(("\n\nRunning TestBiSSEIncompleteSamplingTrees.test_bisse_"
               "sampled_terminal_counts_incomplete_sampling "
                "(second 100 batches)"))
        sim_batches_incomplete = list()
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

            sim_batches_incomplete.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)

        batch_avgs_extant_sampled_terminals_incomplete_n0 = list()
        batch_avgs_extant_sampled_terminals_incomplete_n1 = list()
        n0_ci_overlap_count = 0
        n1_ci_overlap_count = 0
        global_mean_n0 = 0.0
        global_mean_n1 = 0.0
        for i, batch in enumerate(sim_batches_incomplete):
            approx_n0s_sampled_incomplete = \
                [ann_tr.extant_terminal_sampled_state_count_dict[0] for ann_tr in batch]
            approx_n1s_sampled_incomplete = \
                [ann_tr.extant_terminal_sampled_state_count_dict[1] for ann_tr in batch]
            
            mean_n0_sampled_incomplete = \
                statistics.mean(approx_n0s_sampled_incomplete)
            mean_n1_sampled_incomplete = \
                statistics.mean(approx_n1s_sampled_incomplete)
            
            batch_avgs_extant_sampled_terminals_incomplete_n0.\
                append(mean_n0_sampled_incomplete)
            batch_avgs_extant_sampled_terminals_incomplete_n1.\
                append(mean_n1_sampled_incomplete)

            global_mean_n0 += mean_n0_sampled_incomplete
            global_mean_n1 += mean_n1_sampled_incomplete

            stdevs_n0 = statistics.stdev(approx_n0s_sampled_incomplete)
            stdevs_n1 = statistics.stdev(approx_n1s_sampled_incomplete)
            
            sterr_n0 = stdevs_n0 / math.sqrt(n_sim)
            sterr_n1 = stdevs_n1 / math.sqrt(n_sim)

            n0_ci_width = 1.96 * sterr_n0
            n1_ci_width = 1.96 * sterr_n1

            if abs(mean_n0_sampled_incomplete \
                   - expected_avg_sampled_incomplete_n0) \
                    <= (n0_ci_width + n0_ci_width_sampled_complete_list[i]):
                n0_ci_overlap_count += 1
            
            if abs(mean_n1_sampled_incomplete \
                   - expected_avg_sampled_incomplete_n1) \
                    <= (n1_ci_width + n1_ci_width_sampled_complete_list[i]):
                n1_ci_overlap_count += 1

            # if (expected_avg_sampled_incomplete_n0 < \
            #     (mean_n0_sampled_incomplete + n0_ci_width)) \
            #     and (expected_avg_sampled_incomplete_n0 > \
            #     (mean_n0_sampled_incomplete - n0_ci_width)):
            #     n0_ci_overlap_count += 1

            # if (expected_avg_sampled_incomplete_n1 < \
            #     (mean_n1_sampled_incomplete + n1_ci_width)) \
            #     and (expected_avg_sampled_incomplete_n1 > \
            #     (mean_n1_sampled_incomplete - n1_ci_width)):
            #     n1_ci_overlap_count += 1

        print("\n\nExpected global mean n0 taxon count = " \
              + str(global_mean_n0 / 100.0))
        print("Simulated (with incomplete sampling) mean n0 taxon count = " \
              + str(expected_avg_sampled_incomplete_n0))
        print("\nExpected global mean n1 taxon count = " \
              + str(global_mean_n1 / 100.0))
        print("Simulated (with incomplete sampling) mean n0 taxon count = " \
              + str(expected_avg_sampled_incomplete_n1))

        # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
        # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
        # 95% of the time

        print("\n95% CI for sampled taxa in state 0 includes expectation " \
            + str(n0_ci_overlap_count) + " times.")
        print("\n95% CI for sampled taxa in state 1 includes expectation " \
            + str(n1_ci_overlap_count) + " times.")
        
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)

        self.assertAlmostEqual(
            n0_ci_overlap_count,
            exp_count,
            msg="Truth should be contained within 95%-CI " + str(exp_count) \
                + " (+/- " + str(a_delta) \
            + ") out of 100 times.",
            delta=a_delta)
        
        self.assertAlmostEqual(
            n1_ci_overlap_count,
            exp_count,
            msg="Truth should be contained within 95%-CI " + str(exp_count) \
                + " (+/- " + str(a_delta) \
            + ") out of 100 times.",
            delta=a_delta)


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
    # $ python3 tests/distribution/test_dn_discrete_sse_bisse_incomplete_sampling.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_dn_discrete_sse_bisse_incomplete_sampling
    #
    # or 
    #
    # $ python3 -m unittest tests.distribution.test_dn_discrete_sse_bisse_incomplete_sampling.TestBiSSEIncompleteSamplingTrees.test_bisse_sampled_terminal_counts_incomplete_sampling

    unittest.main()