import unittest

# pj imports
import phylojunction.utility.exception_classes as ec
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.distribution.dn_discrete_sse as distsse

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestSSEStopConditionsFBD(unittest.TestCase):

    @classmethod
    def setUp(cls):
        # not state-dependent (just state 0, and no transition)
        rates_t0_s0 = [
            sseobj.DiscreteStateDependentRate(
                name="lambda",
                val=1.0,
                event=sseobj.MacroevolEvent.W_SPECIATION,
                states=[0,0,0]),
            sseobj.DiscreteStateDependentRate(
                name="mu",
                val=0.5,
                event=sseobj.MacroevolEvent.EXTINCTION,
                states=[0]),
            sseobj.DiscreteStateDependentRate(
                name="psi",
                val=0.8,
                event=sseobj.MacroevolEvent.ANCESTOR_SAMPLING,
                states=[0])
        ]
    
        matrix_atomic_rate_params = [rates_t0_s0]
        
        state_dep_par_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                matrix_atomic_rate_params,
                1)

        event_handler = sseobj.MacroevolEventHandler(state_dep_par_manager)

        cls.sse_stash = sseobj.SSEStash(event_handler)

    def test_tree_size_stop_condition_origin_fbd(self):
        """
        Test if fossilized birth-death trees have correct number of
        extant taxa, as specified by stop condition (input).
        Trees start from origin.
        """

        # setting up stopping conditions
        stop_condition = "size"
        stop_condition_value = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        start_at_origin = True

        # simulation initialization
        n_sim = 10
        start_states_list = [0 for i in range(n_sim)]
    
        # simulations
        dnsse = distsse.DnSSE(
            self.sse_stash,
            n=n_sim,
            origin=start_at_origin,
            start_states_list=start_states_list,
            stop=stop_condition,
            stop_value=stop_condition_value,
            condition_on_survival=True,
            epsilon=1e-12,
            debug=False)

        print(("\n\nRunning TestSSEStopConditionsFBD.test_tree_size_stop_"
               "condition_origin_fbd"))
        trs = dnsse.generate()

        tr_sizes = [ann_tr.n_extant_terminal_nodes for ann_tr in trs]

        self.assertEqual(tr_sizes, stop_condition_value)

    def test_tree_size_stop_condition_root_fbd(self):
        """
        Test if fossilized birth-death trees have correct number of
        extant taxa, as specified by stop condition (input).
        Trees start from root.
        """
        
        # setting up stopping conditions
        stop_condition = "size"
        stop_condition_value = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        start_at_origin = False

        # simulation initialization
        n_sim = 10
        start_states_list = [0 for i in range(n_sim)]
    
        # simulations
        dnsse = distsse.DnSSE(self.sse_stash,
            n=n_sim,
            origin=start_at_origin,
            start_states_list=start_states_list,
            stop=stop_condition,
            stop_value=stop_condition_value,
            condition_on_survival=True,
            epsilon=1e-12)
        
        print(("\n\nRunning TestSSEStopConditionsFBD.test_tree_size_stop_"
               "condition_root_fbd"))
        trs = dnsse.generate()

        tr_sizes = [ann_tr.n_extant_terminal_nodes for ann_tr in trs]

        self.assertEqual(tr_sizes, stop_condition_value)

    def test_tree_height_stop_condition_origin_fbd(self):
        """
        Test if fossilized birth-death trees have correct tree height,
        as specified by stop condition (input).
        Trees start from origin.
        """
        
        # setting up stopping conditions
        stop_condition = "age"
        stop_condition_value = \
            [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        start_at_origin = True

        # simulation initialization
        n_sim = 10
        start_states_list = [0 for i in range(n_sim)]
        # seeds_list = [i+1 for i in range(n_sim)]
    
        # simulations
        dnsse = distsse.DnSSE(
            self.sse_stash,
            n=n_sim,
            origin=start_at_origin,
            start_states_list=start_states_list,
            stop=stop_condition,
            stop_value=stop_condition_value,
            condition_on_survival=True,
            epsilon=1e-12)

        print(("\n\nRunning TestSSEStopConditionsFBD.test_tree_height_"
               "stop_condition_origin_fbd"))
        trs = dnsse.generate()

        tr_sizes = [ann_tr.origin_age for ann_tr in trs]

        for idx, tr_size in enumerate(tr_sizes):
            self.assertAlmostEqual(tr_size, stop_condition_value[idx], delta=1e-12)

    def test_tree_height_stop_condition_root_fbd(self):
        """
        Test if fossilized birth-death trees have correct tree height,
        as specified by stop condition (input).
        Trees start from root.
        """

        # setting up stopping conditions
        stop_condition = "age"
        stop_condition_value = \
            [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        start_at_origin = False

        # simulation initialization
        n_sim = 10
        start_states_list = [0 for i in range(n_sim)]
        # seeds_list = [i+1 for i in range(n_sim)]
    
        # simulations
        dnsse = distsse.DnSSE(
            self.sse_stash,
            n=n_sim,
            origin=start_at_origin,
            start_states_list=start_states_list,
            stop=stop_condition,
            stop_value=stop_condition_value,
            condition_on_survival=True,
            epsilon=1e-12)
        
        print(("\n\nRunning TestSSEStopConditionsFBD.test_tree_height_"
               "stop_condition_root_fbd"))
        trs = dnsse.generate()

        tr_sizes = [ann_tr.root_age for ann_tr in trs]

        for idx, tr_size in enumerate(tr_sizes):
            self.assertAlmostEqual(tr_size, stop_condition_value[idx], delta=1e-12)


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
    # on the terminal, remember to add "src/phylojunction" to
    # PYTHONPATH (system variable), or to set it if it does not
    # exist -- don't forget to export it!
    # 
    # Then you can do:
    # $ python3 tests/distribution/test_dn_discrete_sse_stop_conditions_fbd.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_dn_discrete_sse_stop_conditions_fbd
    #
    # or 
    #
    # $ python3 -m unittest tests.distribution.test_dn_discrete_sse_stop_conditions_fbd.TestSSEStopConditionsFBD.test_tree_size_stop_condition_origin_fbd

    unittest.main()