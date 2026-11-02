import unittest
import copy
from unittest.mock import patch

# pj imports
import phylojunction.utility.helper_functions as pjh
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.distribution.dn_discrete_sse as distsse

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestDnSSEObject(unittest.TestCase):

    # Explicit IID draws share a starting state but do not mutate caller-owned rates.
    # Statistical tree tests do not isolate argument reuse or state selection.
    def test_start_state_expansion(self):
        rate = sseobj.DiscreteStateDependentRate(
            name='lambda', val=1.0, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0])
        stash = sseobj.SSEStash(sseobj.MacroevolEventHandler(
            sseobj.DiscreteStateDependentParameterManager([[rate]], 2)))
        for state in (0, 1):
            sim = distsse.DnSSE(stash, start_state=state, stop='age', stop_value=1.0)
            with patch.object(sim, 'simulate', return_value=object()) as simulate, \
                    patch.object(sim, '_is_tr_ok', return_value=True):
                self.assertEqual(len(sim.draw_many(2)), 2)
                self.assertEqual([call.args[0] for call in simulate.call_args_list], [state, state])
            self.assertEqual(rate.value, [1.0])

    def test_dnsse_vectorization(self):
        """Test explicit scalar parameter settings across separate tree draws."""

        # NOTE
        # first 5 trees must be alive (l = 1.0, mu = 0.0)
        # second 5 trees must be dead (l = 0.0, mu = 1.0)
        total_n_states = 2

        l = [1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        lrate = \
            sseobj.DiscreteStateDependentRate(
                name="lambda0",
                val=l,
                event=sseobj.MacroevolEvent.W_SPECIATION,
                states=[0,0,0])
        
        mu = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        murate = \
            sseobj.DiscreteStateDependentRate(
                name="mu0",
                val=mu,
                event=sseobj.MacroevolEvent.EXTINCTION,
                states=[0])
        
        rates_t0 = [lrate, murate]

        matrix_atomic_rate_params = [rates_t0]
        
        state_dep_par_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                matrix_atomic_rate_params, \
                total_n_states)

        meh = sseobj.MacroevolEventHandler(state_dep_par_manager)

        sse_stash = sseobj.SSEStash(meh)

        ########
        # Tree #
        ########
        n_sim = 10
        stop_condition = "age"
        stop_condition_value = [4.0] # 4.0 time units
        start_at_origin = True
        start_states_list = [0 for i in range(n_sim)]

        trs = []
        for i in range(n_sim):
            stash = copy.deepcopy(sse_stash)
            for epoch in stash.meh.sse_rate_manager.matrix_state_dep_params:
                for parameter in epoch:
                    parameter.value = [parameter.value[i]]
            sim = distsse.DnSSE(stash, origin=start_at_origin, start_state=start_states_list[i],
                               stop=stop_condition, stop_value=stop_condition_value[0],
                               condition_on_speciation=False, condition_on_survival=False,
                               epsilon=1e-12, runtime_limit=3600, debug=False)
            # Force an event before the horizon so pure-death cases exercise extinction
            # deterministically, instead of relying on five exponential draws all being short.
            with patch('phylojunction.distribution.dn_parametric.DnExponential.draw_exp', return_value=[1.0]):
                trs.append(sim.draw())

        trs1to5 = [ tr for tr in trs[:5] ]
        trs6to10 = [ tr for tr in trs[5:] ]

        for tr in trs1to5:
            # first 5 trees are Yule trees and cannot have died
            self.assertFalse(tr.tree_died)

        for tr in trs6to10:
            # second 5 trees must have died
            self.assertTrue(tr.tree_died)


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
    # $ python3 tests/distribution/test_dn_discrete_sse_object.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_dn_discrete_sse_object
    #
    # or 
    #
    # $ python3 -m unittest tests.distribution.test_dn_discrete_sse_object.TestDnSSEObject.test_dnsse_vectorization

    unittest.main()