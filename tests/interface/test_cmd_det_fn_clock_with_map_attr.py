import unittest
from pathlib import Path

# pj imports
import phylojunction.interface.cmdbox.cmd_parse as cmd

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"

TREES_MAPS_DIR = Path(__file__).parent.parent.parent / "examples" / \
    "trees_maps_files"
TR_FILE_PATH_STR = str(TREES_MAPS_DIR / "geosse_dummy_tree2.tre")
MAPS_FILE_PATH_STR = str(TREES_MAPS_DIR / "geosse_dummy_tree2_maps.tsv")
TIP_STATES_FILE_PATH_STR = \
    str(TREES_MAPS_DIR / "geosse_dummy_tree2_tip_states.tsv")


class TestGammaSpikeClockWithMapAttr(unittest.TestCase):
    """Test 'gamma_spike_clock' fed a 'map_attr(fun="smap", ...)' tree.

    'map_attr' lets an externally-computed stochastic map "paint"
    discrete-state changes (anagenetic and cladogenetic, mid-branch)
    onto a tree read in from Newick ('read_tree(..., node_name_attr=
    ...)'). These tests confirm that painted tree can then be handed
    straight to 'gamma_spike_clock' (mirroring
    'test_cmd_n_rec_branches_and_n_nr_from.py's own
    'br_specific_rel_rate' pipeline, just fed a mapped/read tree
    instead of a 'discrete_sse'-simulated one), and that the state
    painting ('at_dict'/'clado_at_dict'/'node_attr_dict') survives
    into the clock model's output tree untouched.

    This combination used to crash (an 'AttributeError' deep inside
    'extract_reconstructed_tree') for any tree with no separate
    origin node (i.e. the seed node is itself the root) -- see
    'tests.readwrite.test_pj_read_annotated_newick.
    TestReadAnnotatedNewick.test_node_names_attribute_root_label_with_no_origin'
    for the root cause and its own, narrower regression test; this
    file instead exercises the full, real-world pipeline this bug
    actually blocked.
    """

    preamble = (
        'trs <- read_tree(n=1, file_path="' + TR_FILE_PATH_STR + '", '
        'node_name_attr="index")\n'
        'mapped_trs := map_attr(tree=trs, fun="smap", '
        'maps_file_path="' + MAPS_FILE_PATH_STR + '", '
        'tip_attr_file_path="' + TIP_STATES_FILE_PATH_STR + '", '
        'attr_name="state", n_regions=2, geo="true")\n'
        'n_br := n_rec_branches(tree=mapped_trs)\n'
        'br_specific_rel_rate ~ lognormal(n_nr_from=n_br, meanlog=-0.125, '
        'sdlog=0.5, log_space="true")\n'
    )

    def test_gamma_spike_clock_runs_end_to_end(self):
        """'gamma_spike_clock' does not crash on a smap-painted tree.

        Also confirms the reconstructed tree is correctly rooted at
        the tree's true root ("nd9", the Newick string's own
        '[&index=9]' seed node) rather than some other, wrong internal
        node -- the direct, visible consequence of the root-label bug
        this combination used to trigger.
        """

        script = self.preamble + (
            'phylo := gamma_spike_clock(tree=mapped_trs, '
            'attr_name="mutation_rate", clock_rate=0.001, spike_mean=0.01, '
            'spike_shape=2.0, br_specific_rel_rate=br_specific_rel_rate)\n')

        dag_obj = cmd.script2dag(script, in_pj_file=False, random_seed=7)
        phylo_node = dag_obj.get_node_dag_by_name("phylo")
        ann_tr = phylo_node.value[0]

        self.assertEqual(ann_tr.tree_reconstructed.seed_node.label, "nd9")
        self.assertEqual(ann_tr.n_rec_branches(), 8)

    def test_br_specific_rel_rate_matches_branch_raw_rate(self):
        """'br_specific_rel_rate' values land exactly as each branch's raw rate.

        Same check as
        'test_cmd_n_rec_branches_and_n_nr_from.py's
        'test_gamma_spike_clock_br_specific_rel_rate_end_to_end', here
        run on a tree read from Newick and mapped, rather than
        'discrete_sse'-simulated -- confirming the clock model itself
        behaves identically regardless of where the tree came from.
        """

        script = self.preamble + (
            'phylo := gamma_spike_clock(tree=mapped_trs, '
            'attr_name="mutation_rate", clock_rate=0.001, spike_mean=0.01, '
            'spike_shape=2.0, br_specific_rel_rate=br_specific_rel_rate, '
            'strict_clock_when_no_stubs="false")\n')

        dag_obj = cmd.script2dag(script, in_pj_file=False, random_seed=7)
        rate_node = dag_obj.get_node_dag_by_name("br_specific_rel_rate")
        phylo_node = dag_obj.get_node_dag_by_name("phylo")

        for i, ann_tr in enumerate(phylo_node.value):
            expected = list(rate_node.value[i])
            idx = 0

            for nd in ann_tr.tree_reconstructed.preorder_node_iter():
                if nd.branch_raw_rate is not None:
                    self.assertEqual(nd.branch_raw_rate, expected[idx])
                    idx += 1

            self.assertEqual(idx, len(expected))

    def test_state_painting_survives_into_clock_output(self):
        """Mid-branch/cladogenetic state changes survive 'gamma_spike_clock'.

        'GammaSpikeClockModel' deep-copies its input tree before
        adding its own 'mutation_rate' attribute (see
        'calculation.clock_model.ClockModel.generate'), so
        'mapped_trs's own state-painting members
        ('at_dict'/'clado_at_dict'/'node_attr_dict', populated by
        'map_attr') must come along for the ride, unchanged in
        content (though not as the very same objects, since they were
        deep-copied) -- demonstrating the two layers (discrete-state
        painting and the relaxed clock) coexist on the same tree
        without one clobbering the other.
        """

        script = self.preamble + (
            'phylo := gamma_spike_clock(tree=mapped_trs, '
            'attr_name="mutation_rate", clock_rate=0.001, spike_mean=0.01, '
            'spike_shape=2.0, br_specific_rel_rate=br_specific_rel_rate)\n')

        dag_obj = cmd.script2dag(script, in_pj_file=False, random_seed=7)
        mapped_ann_tr = dag_obj.get_node_dag_by_name("mapped_trs").value[0]
        ann_tr = dag_obj.get_node_dag_by_name("phylo").value[0]

        # sanity check on the fixture itself: if this ever becomes
        # empty, the rest of this test exercises nothing
        self.assertGreater(len(mapped_ann_tr.at_dict), 0)

        self.assertEqual(len(ann_tr.at_dict), len(mapped_ann_tr.at_dict))
        self.assertEqual(
            len(ann_tr.clado_at_dict), len(mapped_ann_tr.clado_at_dict))

        for nd_name, at_list in mapped_ann_tr.at_dict.items():
            cloned_at_list = ann_tr.at_dict[nd_name]
            self.assertEqual(len(cloned_at_list), len(at_list))

            for at, cloned_at in zip(at_list, cloned_at_list):
                self.assertIsNot(at, cloned_at)
                self.assertEqual(at.attr_label, cloned_at.attr_label)
                self.assertEqual(at.from_state, cloned_at.from_state)
                self.assertEqual(at.to_state, cloned_at.to_state)
                self.assertEqual(at.global_time, cloned_at.global_time)

        # discrete states ('state', from the stochastic map) and the
        # continuous clock attribute ('mutation_rate') both populated,
        # side by side, on the very same reconstructed-tree node
        checked_a_real_branch = False
        for rec_nd in ann_tr.tree_reconstructed.preorder_node_iter():
            if rec_nd.branch_raw_rate is not None:
                self.assertIn("state", ann_tr.node_attr_dict[rec_nd.label])
                self.assertIsNotNone(rec_nd.mutation_rate)
                checked_a_real_branch = True

        self.assertTrue(checked_a_real_branch)


if __name__ == "__main__":
    unittest.main()
