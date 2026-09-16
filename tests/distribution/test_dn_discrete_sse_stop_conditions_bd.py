import unittest
from types import SimpleNamespace

# pj imports
import phylojunction.utility.exception_classes as ec
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.distribution.dn_discrete_sse as distsse

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"

class TestSSEStopConditionsBD(unittest.TestCase):

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
                val=0.9,
                event=sseobj.MacroevolEvent.EXTINCTION,
                states=[0])
        ]
    
        matrix_atomic_rate_params = [rates_t0_s0]
        
        state_dep_par_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                matrix_atomic_rate_params,
                1)

        event_handler = sseobj.MacroevolEventHandler(state_dep_par_manager)
        
        cls.sse_stash = sseobj.SSEStash(event_handler)


    # Aborted trajectories must not enter accepted output under any conditioning mode.
    # Existing size/age tests do not force an abort; remove if aborts no longer return trees.
    def test_reject_aborted_trees(self):
        for origin in [False, True]:
            for survival in [False, True]:
                for stop, value in [("age", 10.0), ("size", 3)]:
                    with self.subTest(origin=origin, survival=survival, stop=stop):
                        sim = distsse.DnSSE(
                            self.sse_stash, origin=origin, start_states_list=[0],
                            stop=stop, stop_value=[value], condition_on_survival=survival)
                        self.assertFalse(sim._is_tr_ok(SimpleNamespace(tree_invalid=True), value))
        sim = distsse.DnSSE(
            self.sse_stash, origin=True, start_states_list=[0], stop="age", stop_value=[10.0],
            abort_at_alive_count=3, max_n_failed_attempts=0, rng_seed=123)
        self.sse_stash.get_meh().sse_rate_manager.matrix_state_dep_params[0][1].value = [0.0]
        with self.assertRaises(ec.MaxNFailedAttemptsLimit):
            sim.generate()

    # Unsupported observation models must fail before size simulation; existing tests
    # use default complete sampling. Remove when size stopping implements incomplete sampling.
    def test_size_requires_complete_sampling(self):
        probabilities = self.sse_stash.get_prob_handler().state_dep_prob_manager.matrix_state_dep_params
        for values in [[0.0], [0.5], [1.0], [1.0, 0.5]]:
            for stop in ["age", "size"]:
                with self.subTest(values=values, stop=stop):
                    probabilities[0][0].value = list(values)
                    kwargs = dict(n=len(values), origin=True, start_states_list=[0],
                                  stop=stop, stop_value=[3.0])
                    if stop == "size" and any(value != 1.0 for value in values):
                        with self.assertRaisesRegex(ec.ObjInitInvalidArgError, "complete sampling"):
                            distsse.DnSSE(self.sse_stash, **kwargs)
                    else:
                        distsse.DnSSE(self.sse_stash, **kwargs)

    def test_tree_size_stop_condition_origin(self):
        """
        Test if birth-death trees have correct number of
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
        # seeds_list = [i+1 for i in range(n_sim)]
    
        # simulations
        dnsse = distsse.DnSSE(
            self.sse_stash,
            n=n_sim,
            origin=start_at_origin,
            stop=stop_condition,
            stop_value=stop_condition_value,
            start_states_list=start_states_list,
            condition_on_survival=True,
            epsilon=1e-12)
        
        print(("\n\nRunning TestSSEStopConditionsBD.test_tree_size_stop_"
               "condition_origin"))
        trs = dnsse.generate()

        tr_sizes = [ann_tr.n_extant_terminal_nodes for ann_tr in trs]

        self.assertEqual(tr_sizes, stop_condition_value)


    def test_tree_size_stop_condition_root(self):
        """
        Test if birth-death trees have correct number of
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
        dnsse = distsse.DnSSE(
            self.sse_stash,
            n=n_sim,
            origin=start_at_origin,
            start_states_list=start_states_list,
            stop=stop_condition,
            stop_value=stop_condition_value,
            condition_on_survival=True,
            epsilon=1e-12)

        print(("\n\nRunning TestSSEStopConditionsBD.test_tree_size_stop_"
               "condition_root"))
        trs = dnsse.generate()

        tr_sizes = [ann_tr.n_extant_terminal_nodes for ann_tr in trs]

        self.assertEqual(tr_sizes, stop_condition_value)


    def test_tree_height_stop_condition_origin(self):
        """
        Test if birth-death trees have correct tree height, as
        specified by stop condition (input).
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
        
        print(("\n\nRunning TestSSEStopConditionsBD.test_tree_height_stop_"
               "condition_origin"))
        trs = dnsse.generate()

        tr_sizes = [ann_tr.origin_age for ann_tr in trs]

        for idx, tr_size in enumerate(tr_sizes):
            self.assertAlmostEqual(
                tr_size,
                stop_condition_value[idx],
                delta=1e-12)


    def test_tree_height_stop_condition_root(self):
        """
        Test if birth-death trees have correct tree height, as
        specified by stop condition (input).
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
        
        print(("\n\nRunning TestSSEStopConditionsBD.test_tree_height_stop_"
               "condition_root"))
        trs = dnsse.generate()

        tr_sizes = [ann_tr.root_age for ann_tr in trs]

        # for ann_tr in trs:
        #     print(ann_tr.tree.as_string(schema="newick", suppress_internal_taxon_labels=True))

        for idx, tr_size in enumerate(tr_sizes):
            self.assertAlmostEqual(tr_size, stop_condition_value[idx], delta=1e-12)


    def test_invalid_stop_condition_value(self):
        """
        Test if invalid stop condition values are correctly caught and
        the appropriate exceptions are raised:
            (1) Negative number of terminal nodes,
            (2) Number of terminal nodes being a float that cannot be
                converted to integer,
            (3) Negative tree height
        """
        # simulation initialization
        n_sim = 2
        start_states_list = [0 for i in range(n_sim)]

        stop_condition = "size"
        stop_condition_value = [1, -1]
        start_at_origin = True

        print(("\n\nRunning TestSSEStopConditionsBD.test_invalid_stop_"
               "condition_value"))
        
        # (1) Negative number of terminal nodes
        with self.assertRaises(ec.ObjInitInvalidArgError) as exc:
            distsse.DnSSE(
                self.sse_stash,
                n=n_sim,
                origin=start_at_origin,
                start_states_list=start_states_list,
                stop=stop_condition,
                stop_value=stop_condition_value,
                condition_on_survival=True,
                epsilon=1e-12)

        self.assertEqual(str(exc.exception),
            ("Could not initialize the object of 'DnSSE'. "
             "The argument provided to 'stop' was invalid. "
             "Stop condition value (number of terminal nodes) "
             "cannot be negative."))

        # (2) Number of terminal nodes being a float that cannot be
        # converted to integer
        stop_condition_value = [0.5, 2]
        with self.assertRaises(ec.ObjInitInvalidArgError) as exc:
            distsse.DnSSE(
                self.sse_stash,
                n=n_sim,
                origin=start_at_origin,
                start_states_list=start_states_list,
                stop=stop_condition,
                stop_value=stop_condition_value,
                condition_on_survival=True,
                epsilon=1e-12)       
        
        self.assertEqual(str(exc.exception),
            ("Could not initialize the object of 'DnSSE'. "
             "The argument provided to \'stop\' was invalid. "
             "Stop condition value (number of terminal "
             "nodes) could not be converted to integer."))

        # (3) Negative tree height
        stop_condition = "age"
        stop_condition_value = [0.0, -0.1]
        with self.assertRaises(ec.ObjInitInvalidArgError) as exc:
            distsse.DnSSE(
                self.sse_stash,
                n=n_sim,
                origin=start_at_origin, 
                start_states_list=start_states_list,
                stop=stop_condition,
                stop_value=stop_condition_value,
                condition_on_survival=True,
                epsilon=1e-12)

        self.assertEqual(str(exc.exception),
            ("Could not initialize the object of \'DnSSE\'. "
             "The argument provided to 'stop' was invalid. "
             "Stop condition value (tree height) cannot be negative."))


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
    # $ python3 tests/distribution/test_dn_discrete_sse_stop_conditions_bd.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_dn_discrete_sse_stop_conditions_bd
    #
    # or 
    #
    # $ python3 -m unittest tests.distribution.test_dn_discrete_sse_stop_conditions_bd.TestSSEStopConditionsBD.test_tree_size_stop_condition_origin

    unittest.main()