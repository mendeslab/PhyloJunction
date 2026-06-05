import unittest

# pj imports
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.utility.exception_classes as ec


__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"

class TestDiscreteSSE(unittest.TestCase):
    
    def test_state_dep_prob(self):
        """
        Test that state-dependent probability parameter is correctly
        initialized
        """

        state_dep_sampling_prob_A = \
            sseobj.DiscreteStateDependentProbability(
                name="rho0_A", val=0.5, state=0
            )

        state_dep_sampling_prob_B = \
            sseobj.DiscreteStateDependentProbability(
                name="rho0_B", val=[0.5, 0.5], state=1
            )

        state_dep_sampling_prob_C = \
            sseobj.DiscreteStateDependentProbability(
                name="rho0_C", val=["0.5", "0.5"], state=0
            )

        expected_str_A = "Discrete state-dependent probability\n" \
            + "   Name:   rho0_A\n" \
            + "   Value:  0.5\n" \
            + "   State:  0\n" \
            + "   Epoch:  1\n\n"

        expected_str_B = "Discrete state-dependent probability\n" \
            + "   Name:   rho0_B\n" \
            + "   Value:  0.5, 0.5\n" \
            + "   State:  1\n" \
            + "   Epoch:  1\n\n"

        expected_str_C = "Discrete state-dependent probability\n" \
            + "   Name:   rho0_C\n" \
            + "   Value:  0.5, 0.5\n" \
            + "   State:  0\n" \
            + "   Epoch:  1\n\n"

        self.assertEqual(str(state_dep_sampling_prob_A),
            expected_str_A)

        self.assertEqual(str(state_dep_sampling_prob_B),
            expected_str_B)

        self.assertEqual(str(state_dep_sampling_prob_C),
            expected_str_C)
            

    def test_state_dep_prob_init_error(self):
        """
        Test that state-dependent probability parameter throws
        exception if value(s) < 0.0 or value(s) > 1.0
        """
        
        with self.assertRaises(ec.NotBetweenZeroAndOneError) as exc1:
            sseobj.DiscreteStateDependentProbability(
                name="rho0_A", val=-0.1, state=0
            )

        with self.assertRaises(ec.NotBetweenZeroAndOneError) as exc2:
            sseobj.DiscreteStateDependentProbability(
                name="rho0_B", val=1.1, state=0
            )

        with self.assertRaises(ec.NotBetweenZeroAndOneError) as exc3:
            sseobj.DiscreteStateDependentProbability(
                name="rho0_C", val=[0.1, -0.1], state=0
            )

        with self.assertRaises(ec.NotBetweenZeroAndOneError) as exc4:
            sseobj.DiscreteStateDependentProbability(
                name="rho0_D", val=[1.0, 1.1], state=0
            )

        self.assertEqual(str(exc1.exception),
            "ERROR: 'rho0_A' must be contained in [0,1], but one or " \
                + "more of its values was negative."
        )

        self.assertEqual(str(exc2.exception),
            "ERROR: 'rho0_B' must be contained in [0,1], but one or " \
                + "more of its values was positive."
        )

        self.assertEqual(str(exc3.exception),
            "ERROR: 'rho0_C' must be contained in [0,1], but one or " \
                + "more of its values was negative."
        )

        self.assertEqual(str(exc4.exception),
            "ERROR: 'rho0_D' must be contained in [0,1], but one or " \
                + "more of its values was positive."
        )


    def test_state_dep_param_manager(self):
        """
        Test probabilities are correctly retrieved by parameter manager
        depending on which time-slice they are
        """
        
        total_n_states = 2

        probs_t0 = [
            sseobj.DiscreteStateDependentProbability(name="rho0_0", val=0.1, state=0),
            sseobj.DiscreteStateDependentProbability(name="rho0_1", val=0.2, state=1),
        ]

        probs_t1 = [
            sseobj.DiscreteStateDependentProbability(name="rho1_0", val=0.2, state=0),
            sseobj.DiscreteStateDependentProbability(name="rho1_1", val=0.4, state=1),
        ]

        # 1D: time slices (i)
        # 2D: all rates from all states in i-th time slice
        matrix_state_dep_probs = [ probs_t0, probs_t1 ]

        state_dep_par_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                matrix_state_dep_probs, total_n_states,
                seed_age_for_time_slicing=2.0,
                list_time_slice_age_ends=[1.0]
            )

        params_matrix_state0 = \
            state_dep_par_manager.state_dep_params_dict[0]
        params_matrix_state1 = \
            state_dep_par_manager.state_dep_params_dict[1]
        
        params_state0_at_time11 = \
            state_dep_par_manager.state_dep_params_at_time(
                1.1, params_matrix=params_matrix_state0)
        params_state1_at_time11 = \
            state_dep_par_manager.state_dep_params_at_time(
                1.1, params_matrix=params_matrix_state1)

        params_state0_at_time09 = \
            state_dep_par_manager.state_dep_params_at_time(
                0.9, params_matrix=params_matrix_state0)
        params_state1_at_time09 = \
            state_dep_par_manager.state_dep_params_at_time(
                0.9, params_matrix=params_matrix_state1)
        
        self.assertEqual(params_state0_at_time11[0].name, "rho1_0")
        self.assertEqual(params_state1_at_time11[0].name, "rho1_1")
        self.assertEqual(params_state0_at_time09[0].name, "rho0_0")
        self.assertEqual(params_state1_at_time09[0].name, "rho0_1")


    def test_state_dep_param_prob_handler(self):
        """
        Test probability handler correctly applies state-dependent
        probability depending on time-slice
        """

        total_n_states = 2

        probs_t0 = [
            sseobj.DiscreteStateDependentProbability(name="rho0_t0", val=0.0, state=0),
            sseobj.DiscreteStateDependentProbability(name="rho1_t0", val=1.0, state=1),
        ]

        probs_t1 = [
            sseobj.DiscreteStateDependentProbability(name="rho0_t1", val=1.0, state=0),
            sseobj.DiscreteStateDependentProbability(name="rho1_t1", val=0.0, state=1),
        ]

        # 1D: time slices (i)
        # 2D: all rates from all states in i-th time slice
        matrix_state_dep_probs = [ probs_t0, probs_t1 ]

        state_dep_par_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                matrix_state_dep_probs, total_n_states,
                seed_age_for_time_slicing=2.0,
                list_time_slice_age_ends=[1.0]
            )

        state_dep_prob_handler = \
            sseobj.DiscreteStateDependentProbabilityHandler(
                state_dep_par_manager
            )

        state0_sampled_t0 = state_dep_prob_handler. \
            randomly_decide_taxon_sampling_at_time_at_state(
                1.1, 0, 0
            )

        state1_sampled_t0 = state_dep_prob_handler. \
            randomly_decide_taxon_sampling_at_time_at_state(
                1.1, 1, 0
            )

        state0_sampled_t1 = state_dep_prob_handler. \
            randomly_decide_taxon_sampling_at_time_at_state(
                0.9, 0, 0
            )

        state1_sampled_t1 = state_dep_prob_handler. \
            randomly_decide_taxon_sampling_at_time_at_state(
                0.9, 1, 0
            )

        self.assertTrue(state0_sampled_t0)
        self.assertFalse(state1_sampled_t0)
        self.assertTrue(state1_sampled_t1)
        self.assertFalse(state0_sampled_t1)


    def test_state_dep_param_manager_init_error(self):
        """
        Test for different types of initialization errors related
        to forgetting state-dependent parameters of certain kinds
        in one or more epochs, or mixing state-dependent parameters
        of different types (e.g., rates vs. probabilities)
        """

        probs_t0 = [
            sseobj.DiscreteStateDependentProbability(name="rho0_t0", val=0.1, state=0),
            sseobj.DiscreteStateDependentProbability(name="rho1_t0", val=0.2, state=1),
        ]

        probs_t1 = [
            sseobj.DiscreteStateDependentProbability(name="rho0_t1", val=0.2, state=0),
        ]

        # issue: states 0 and 2 here, but in epoch 0
        # there were just states 0 and 1 (missing 2
        # in first epoch and 1 in second epoch)
        probs_t1_2 = [
            sseobj.DiscreteStateDependentProbability(name="rho0_t1", val=0.2, state=0),
            sseobj.DiscreteStateDependentProbability(name="rho2_t1", val=0.4, state=2)
        ]

        # issue: repeated state
        probs_t1_3 = [
            sseobj.DiscreteStateDependentProbability(name="rho0_t1", val=0.2, state=0),
            sseobj.DiscreteStateDependentProbability(name="rho0_t1", val=0.4, state=0)
        ]

        # issue: mixing probabilities and rates
        prob_rate_t1 = [
            sseobj.DiscreteStateDependentProbability(name="rho0_t1", val=0.2, state=0),
            sseobj.DiscreteStateDependentRate(name="lambda0_t1", val=1.0,
                event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0])
        ]

        rates_t0 = [
            sseobj.DiscreteStateDependentRate(name="lambda0_t0", val=1.0,
                event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
            sseobj.DiscreteStateDependentRate(name="lambda1_t0", val=1.0,
                event=sseobj.MacroevolEvent.W_SPECIATION, states=[1,1,1])
        ]

        # issue: repeated states
        rates_t1_1 = [
            sseobj.DiscreteStateDependentRate(name="lambda0_t1", val=1.0,
                event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
            sseobj.DiscreteStateDependentRate(name="lambda0_t1", val=1.0,
                event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0])
        ]

        # issue: no rates in epoch!
        rates_t1_2 = []

        # issue: all rates in this epoch are zero, this is not allowed
        rates_t1_3 = [
            sseobj.DiscreteStateDependentRate(name="lambda0_t1", val=0.0,
                event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
            sseobj.DiscreteStateDependentRate(name="lambda1_t1", val=0.0,
                event=sseobj.MacroevolEvent.W_SPECIATION, states=[1,1,1])
        ]

        # 1D: time slices (i)
        # 2D: all rates from all states in i-th time slice
        matrix_state_dep_probs_error1 = [ probs_t0, probs_t1 ]
        matrix_state_dep_probs_error2 = [ probs_t0, prob_rate_t1 ]
        matrix_state_dep_probs_error3 = [ probs_t0, probs_t1_3 ]
                    
        matrix_state_dep_rates_error1 = [ rates_t0, rates_t1_1 ]
        matrix_state_dep_rates_error2 = [ rates_t0, rates_t1_2 ]
        matrix_state_dep_rates_error3 = [ rates_t0, rates_t1_3 ]

        # first epoch there are two probabilities (for state 0 and 1)
        # but second epoch has only one probability (for state 0);
        # we need 'total_n_states' probabilities per time slice!
        total_n_states = 2
        with self.assertRaises(ec.ObjInitIncorrectDimensionError) as exc1:
            state_dep_par_manager = \
                sseobj.DiscreteStateDependentParameterManager(
                    matrix_state_dep_probs_error1,
                    total_n_states
                )

        self.assertEqual(str(exc1.exception),
            ("Could not initialize the object of 'sse_stash' ('DiscreteState"
             "DependentParameterManager' in the backend). Incorrect "
             "dimension of container 'flat_prob_mat' ('matrix_state_dep_"
             "params[1]' in the backend), which was of size 1. The expected"
             " dimension was 2."))

        # rate and prob both passed together
        # but only one type allowed
        total_n_states = 2
        with self.assertRaises(ec.ObjInitRequireSameParameterTypeError) \
            as exc2:
            state_dep_par_manager = \
                sseobj.DiscreteStateDependentParameterManager(
                    matrix_state_dep_probs_error2,
                    total_n_states,
                    seed_age_for_time_slicing=2.0,
                    list_time_slice_age_ends=[1.0]
                )

        self.assertEqual(str(exc2.exception),
            ("When specifying object DiscreteStateDependentParameterManager "
             "only one type of parameter is allowed. Found 2."))

        total_n_states = 2
        with self.assertRaises(
            ec.ObjInitRepeatedStateDependentParameterError) as exc3:
            sseobj.DiscreteStateDependentParameterManager(
                matrix_state_dep_probs_error3,
                total_n_states,
                seed_age_for_time_slicing=2.0,
                list_time_slice_age_ends=[1.0]
            )

        self.assertEqual(str(exc3.exception),
            "State-dependent parameter defined by states 0 were repeated " \
            + "in epoch 1."
        )

        with self.assertRaises(
            ec.ObjInitRepeatedStateDependentParameterError) as exc4:
            sseobj.DiscreteStateDependentParameterManager(
                matrix_state_dep_rates_error1, total_n_states,
                seed_age_for_time_slicing=2.0,
                list_time_slice_age_ends=[1.0]
            )

        self.assertEqual(str(exc4.exception),
            "State-dependent parameter defined by states " \
            + "(event: W_SPECIATION) (0, 0, 0) were repeated in epoch 1."
        )

        total_n_states = 2
        with self.assertRaises(ec.ObjInitIncorrectDimensionError) as exc5:
            sseobj.DiscreteStateDependentParameterManager(
                matrix_state_dep_rates_error2,
                total_n_states,
                seed_age_for_time_slicing=2.0,
                list_time_slice_age_ends=[1.0]
            )

        self.assertEqual(str(exc5.exception),
            ("Could not initialize the object of 'sse_stash' ('DiscreteState"
             "DependentParameterManager' in the backend). Incorrect dimension"
             " of container 'flat_rate_mat' ('matrix_state_dep_params[2]' in"
             " the backend), which was of size 0. The expected dimension was"
             " at least 1.")
        )

        total_n_states = 2
        with self.assertRaises(
            ec.ObjInitRequireNonZeroStateDependentParameterError) as exc6:
            sseobj.DiscreteStateDependentParameterManager(
                matrix_state_dep_rates_error3,
                total_n_states,
                seed_age_for_time_slicing=2.0,
                list_time_slice_age_ends=[1.0]
            )

        self.assertEqual(str(exc6.exception),
            "When specifying object DiscreteStateDependentParameterManager, "
            "one of its dimensions (1) had zero-valued state-dependent "
            "parameters. At least one non-zero parameter must be provided."
        )


    def test_state2pattern_converter(self):
        """
        Test if conversion between integers representing compound
        states and bit patterns is correct
        """
        
        n_characters = 3 # regions A, B, C
        n_states_per_char = 2 # presence/absence (bit)
        # compound-states integer coding -> bit pattern
        # Null: 0 -> 000
        # A:    1 -> 100
        # B:    2 -> 010
        # C:    3 -> 001
        # AB:   4 -> 110
        # AC:   5 -> 101
        # BC:   6 -> 011
        # ABC:  7 -> 111
        
        svc = sseobj.StateIntoPatternConverter(
            n_characters, n_states_per_char)
        
        self.assertEqual(list(svc.int2set_dict.values()),
            ["000", "100", "010", "001", "110", "101", "011", "111"])
        self.assertEqual(list(svc.set2int_dict.values()),
        ["0", "1", "2", "3", "4", "5", "6", "7"])


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
    # $ python3 tests/calculation/test_discrete_sse.py
    # 
    # or
    #
    # $ python3 -m tests.calculation.test_discrete_sse
    #
    # or 
    #
    # $ python3 -m unittest tests.calculation.test_discrete_sse.TestDiscreteSSE.test_state_dep_param_manager_init_error

    unittest.main()