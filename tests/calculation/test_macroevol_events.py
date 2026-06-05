import unittest

# pj imports
import phylojunction.calculation.discrete_sse as sseobj

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"

class TestMacroEvolEvent(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # birth-death
        total_n_states = 1
        cls.bd_rates_t0_s0 = [
            sseobj.DiscreteStateDependentRate(
                name="lambda", val=1.0,
                event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]
            ),
            sseobj.DiscreteStateDependentRate(
                name="mu", val=0.5,
                event=sseobj.MacroevolEvent.EXTINCTION, states=[0]
            )
        ]
        
        # 1D: time slices
        # 2D: states
        # 3D: parameters of state, several parameters -> matrix
        bd_matrix_atomic_rate_params = [ cls.bd_rates_t0_s0 ]
        bd_discrete_state_dep_param_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                bd_matrix_atomic_rate_params, total_n_states
            )
        cls.bd_event_handler = \
            sseobj.MacroevolEventHandler(bd_discrete_state_dep_param_manager)
        cls.bd_state_representation_dict = \
            { 0: ["nd3", "nd4", "nd5", "nd6"] }

        # BiSSE
        total_n_states = 2
        cls.bisse_rates_t0_s0 = [
            sseobj.DiscreteStateDependentRate(
                name="lambda0", val=0.5,
                event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]
            ),
            sseobj.DiscreteStateDependentRate(
                name="mu0", val=0.25,
                event=sseobj.MacroevolEvent.EXTINCTION, states=[0]
            ),
            sseobj.DiscreteStateDependentRate(
                name="q01", val=0.5,
                event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION,
                states=[0,1]
            )
        ]
        cls.bisse_rates_t0_s1 = [
            sseobj.DiscreteStateDependentRate(
                name="lambda1", val=1.5,
                event=sseobj.MacroevolEvent.W_SPECIATION,
                states=[1,1,1]
            ),
            sseobj.DiscreteStateDependentRate(
                name="mu1", val=0.25,
                event=sseobj.MacroevolEvent.EXTINCTION,
                states=[1]
            ),
            sseobj.DiscreteStateDependentRate(
                name="q10", val=0.5,
                event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION,
                states=[1,0]
            )
        ]
        cls.bisse_probs_t0_s0 = [
            sseobj.DiscreteStateDependentProbability(
                name="rho0", val=0.0,
                state=0
            )
        ]
        cls.bisse_probs_t0_s1 = [
            sseobj.DiscreteStateDependentProbability(
                name="rho1", val=1.0,
                state=1
            )
        ]

        bisse_rates_t0 = cls.bisse_rates_t0_s0 + cls.bisse_rates_t0_s1
        bisse_probs_t0 = cls.bisse_probs_t0_s0 + cls.bisse_probs_t0_s1
        
        # 1D: time slices (i)
        # 2D: all rates from all states in i-th time slice
        bisse_matrix_state_dep_rates = [ bisse_rates_t0 ]
        bisse_state_dep_rates_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                bisse_matrix_state_dep_rates, total_n_states
            )
        cls.bisse_event_handler = \
            sseobj.MacroevolEventHandler(bisse_state_dep_rates_manager)
        
        bisse_matrix_state_dep_probs = [ bisse_probs_t0 ]
        bisse_state_dep_probs_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                bisse_matrix_state_dep_probs, total_n_states
            )
        
        cls.bisse_state_dep_prob_handler = \
            sseobj.DiscreteStateDependentProbabilityHandler(
                bisse_state_dep_probs_manager
            )

        cls.bisse_state_representation_dict = \
            { 0: ["nd3"], 1: ["nd4", "nd5", "nd6"] }


    def test_total_rate_single_epoch_bd(self):
        """
        Test if birth-death total rates are correctly computed for single epoch
        """

        a_time = 0.0

        # should not matter what time
        total_rate_t0, state_total_rates = \
            self.bd_event_handler.total_rate(
                a_time,
                self.bd_state_representation_dict)
        total_rate_t100, state_total_rates = \
            self.bd_event_handler.total_rate(
                a_time,
                self.bd_state_representation_dict)

        # 4 * 1.5 = 6.0
        self.assertEqual(
            total_rate_t0,
            6.0,
            "Total rate at time 0.0 should be 6.0.")
        self.assertEqual(
            total_rate_t100,
            6.0,
            "Total rate at time 100.0 should be 6.0.")

    def test_total_rate_single_epoch_bisse(self):
        """
        Test if global and state-conditioned total rates are correctly calculated
        """

        a_time = 0.0

        # overall rate for drawing time to next event
        total_rate_t0, state_total_rates = \
            self.bisse_event_handler.total_rate(
                a_time,
                self.bisse_state_representation_dict)
        total_rate_t100, state_total_rates = \
            self.bisse_event_handler.total_rate(
                a_time,
                self.bisse_state_representation_dict) # should not matter what time

        # 1 * (0.5 + 0.25 + 0.5) + 3 * (1.5 + 0.25 + 0.5) = 8.0
        self.assertEqual(total_rate_t0,
                         8.0,
                         "Total rate at time 0.0 should be 8.0.")
        self.assertEqual(total_rate_t100,
                         8.0,
                         "Total rate at time 100.0 should be 8.0.")

        # 0.5 + 0.25 + 0.5 = 1.25
        self.assertEqual(
            state_total_rates[0],
            1.25,
            "Total rate of state 0 at time 0.0 should be 1.25.")
        self.assertEqual(
            state_total_rates[0],
            1.25,
            "Total rate of state 0 at time 100.0 should be 1.25.")

        # 1.5 + 0.25 + 0.5 = 2.25
        self.assertEqual(
            state_total_rates[1],
            2.25,
            "Total rate of state 1 at time 0.0 should be 2.25.")
        self.assertEqual(
            state_total_rates[1],
            2.25,
            "Total rate of state 1 at time 100.0 should be 2.25.")

    def test_total_rate_two_epochs(self):
        """
        Test if total rate is correctly computed for two epochs
        """
        pass

    def test_event_sampling_bd(self):
        """
        Test if birth/death events are being randomly drawn
        proportionally to their weights.
        """
        
        a_time = 0.0
        n_events = n_samples = 10000
        event_outcomes = list()
        while n_samples > 0:
            i = n_events - n_samples
            n_samples -= 1
            total_rate, \
                state_total_rates = \
                self.bd_event_handler.total_rate(
                    a_time,
                    self.bd_state_representation_dict)

            event_outcomes.append(
                self.bd_event_handler \
                    .sample_event_sse_rate_param(
                    state_total_rates[0],
                    a_time,
                    [0])[0].event.value
            )

        # first [0] indexing gets lambda parameter, second [0] gets the first item inside the vectorized value
        expected_proportion_event_1 = \
            self.bd_rates_t0_s0[0].value[0] / (self.bd_rates_t0_s0[0].value[0] + self.bd_rates_t0_s0[1].value[0])
        obs_proportion_event_1 = \
            sum(1 for v in event_outcomes if v == 0) / n_events

        # let's see events if really drawn proportional to weights
        # print(expected_proportion_event_1)
        # print(obs_proportion_event_1) # pretty good!
        
        self.assertAlmostEqual(
            obs_proportion_event_1,
            expected_proportion_event_1,
            msg="Proportion of event is off by more than 0.01",
            delta=0.01)

    def test_event_sampling_bisse(self):
        """
        Test if state-dependent birth/death events and transition events are being
        randomly drawn proportionally to their weights and tree representation
        """
        a_time = 0.0

        n_events = n_samples = 10000
        seeds = [i for i in range(n_events)] # let's set the seed
        event_outcomes_s0, event_outcomes_s1 = list(), list()
        total_rate, state_total_rates = \
            self.bisse_event_handler.total_rate(
                a_time,
                self.bisse_state_representation_dict)
        while n_samples > 0:
            i = n_events - n_samples

            n_samples -= 1

            event_outcomes_s0.append(
                self.bisse_event_handler.\
                    sample_event_sse_rate_param(
                        state_total_rates[0],
                        a_time,
                        [0],
                        a_seed=seeds[i])[0].event.value
            )

            event_outcomes_s1.append(
                self.bisse_event_handler.\
                    sample_event_sse_rate_param(
                        state_total_rates[1],
                        a_time, [1],
                        a_seed=seeds[i])[0].event.value
            )

        # first [0] indexing gets lambda parameter,
        # second [0] gets the first item inside the vectorized value
        expected_proportion_birth_s0 = \
            self.bisse_rates_t0_s0[0].value[0] / state_total_rates[0]
        obs_proportion_birth_s0 = \
            sum(1 for v in event_outcomes_s0 if v == 0) / n_events
        # v = 0 is the value assigned to speciation (see tpsimulator_classes)

        expected_proportion_birth_s1 = \
            self.bisse_rates_t0_s1[0].value[0] / state_total_rates[1]
        obs_proportion_birth_s1 = \
            sum(1 for v in event_outcomes_s1 if v == 0) / n_events
        # v = 0 is the value assigned to speciation

        # let's see events if really drawn proportional to weights
        # print(expected_proportion_birth_s0)
        # print(obs_proportion_birth_s0) # pretty good!
        # print(expected_proportion_birth_s1)
        # print(obs_proportion_birth_s1) # pretty good!
        
        self.assertAlmostEqual(
            obs_proportion_birth_s0,
            expected_proportion_birth_s0,
            msg="Proportion of event is off by more than 0.01",
            delta=0.01)
        self.assertAlmostEqual(
            obs_proportion_birth_s1,
            expected_proportion_birth_s1,
            msg="Proportion of event is off by more than 0.01",
            delta=0.01)

    def test_sampling_probability_bisse(self):
        pass
        # print(self.bisse_state_dep_prob_handler.state_dep_prob_manager.state_dep_params_dict[0])
        # print(self.bisse_state_dep_prob_handler.randomly_decide_taxon_sampling_at_time_at_state(0.0, 0, 0))
        # print(self.bisse_state_dep_prob_handler.randomly_decide_taxon_sampling_at_time_at_state(0.0, 1, 0))


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
    # on the terminal, remember to add "src/" to
    # PYTHONPATH (system variable), or to set it if it does not
    # exist -- don't forget to export it!
    # 
    # Then you can do:
    # $ python3 tests/calculation/test_macroevol_events.py
    # 
    # or
    #
    # $ python3 -m tests.calculation.test_macroevol_events
    #
    # or 
    #
    # $ python3 -m unittest tests.calculation.test_macroevol_events.TestMacroEvolEvent.test_total_rate_single_epoch_bd
    
    unittest.main()