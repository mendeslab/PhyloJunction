import unittest
from unittest.mock import patch
import numpy as np
from dendropy import Tree, Node, Taxon

# pj imports
import phylojunction.readwrite.pj_read as pjr
import phylojunction.data.tree as pjtr
from phylojunction.data.trait import DiscreteTrait
import phylojunction.calculation.gamma_spike_clock_model as pjgsc


class TestGammaSpike(unittest.TestCase):

    # Valid-clock tests cannot detect partial mutation on invalid input. Keep this guard
    # while rates are populated in place and randomness comes from numpy.random.gamma.
    def test_invalid_rates_do_not_draw_or_annotate(self):
        tree = pjr.read_nwk_tree_str("(A:1,B:1);", in_file=False)
        rec = tree.extract_reconstructed_tree()
        valid = dict(clock_rate=0.1, spike_mean=0.01, spike_shape=2.0,
                     br_specific_rel_rates=[1.0, 1.0])
        invalid = [(name, value) for name in ("clock_rate", "spike_mean", "spike_shape")
                   for value in (-1.0, float("nan"), float("inf"), -float("inf"))]
        invalid += [("clock_rate", 0.0), ("spike_shape", 0.0)]
        invalid += [("br_specific_rel_rates", rates) for rates in
                    ([1.0], [1.0, -1.0], [1.0, float("nan")], [1.0, float("inf")])]
        for name, value in invalid:
            with self.subTest(name=name, value=value), patch("numpy.random.gamma") as gamma:
                with self.assertRaises(ValueError):
                    pjgsc.populate_gamma_spike_mol_rates_on_rec_tr_branches(
                        tree, **(valid | {name: value}))
                gamma.assert_not_called()
                self.assertTrue(all(not hasattr(nd, "mol_rate")
                                    for nd in rec.preorder_node_iter()))
        # Both zero components are legal, including on a branch leaving a speciation.
        pjgsc.populate_gamma_spike_mol_rates_on_rec_tr_branches(
            tree, **(valid | {"spike_mean": 0.0, "br_specific_rel_rates": [0.0, 0.0]}))
        self.assertTrue(all(nd.mol_rate == 0.0 for nd in rec.leaf_node_iter()))

    # The numerical tests do not exercise distribution objects or model repetition.
    # This protects value-based construction and independent copies until that interface changes.
    def test_distribution_expression_and_copies(self):
        from phylojunction.interface.cmdbox.cmd_parse import compile_model
        source = 'trees <- read_tree(string="(A:1,B:1);")\nrates ~ [iid(n_rec_branches(t), lognormal(0, 0.1)) for t in trees]\nd := gamma_spike_clock(trees[1], "rate", 0.1, 0, 2, rates[1])\nout ~ iid(2, d)\n'
        result = compile_model(source).execute(runs=2, random_seed=1)
        for replicate in range(2):
            original = result.value('trees', replicate)[0]
            first, second = result.value('out', replicate)
            self.assertIsNot(first, second)
            self.assertIsNone(original.clock_model_attr_name)
            for node, rate in zip(first.tree_reconstructed.leaf_node_iter(), result.value('rates', replicate)[0]):
                self.assertAlmostEqual(node.rate, 0.1 * rate)
            self.assertIn('rate=', first.rec_str())

    def test_gamma_spike(self):
        # We will construct the complete tree here

        origin_node = Node(
            taxon=Taxon(label="origin"),
            label="origin",
            edge_length=0.0
        )
        origin_node.state = 0
        origin_node.alive = False
        origin_node.sampled = False
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        root_node = Node(
            taxon=Taxon(label="root"),
            label="root",
            edge_length=1.0
        )
        root_node.state = 0
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        origin_node.add_child(root_node)

        # upper child of root
        nd3 = Node(
            taxon=Taxon(label="nd3"),
            label="nd3",
            edge_length=1.0
        )
        nd3.state = 0
        nd3.alive = False
        nd3.sampled = False
        nd3.is_sa = False
        nd3.is_sa_dummy_parent = False
        nd3.is_sa_lineage = False

        root_node.add_child(nd3)

        # upper child of nd3
        nd4 = Node(
            taxon=Taxon(label="nd4"),
            label="nd4",
            edge_length=1.0
        )
        nd4.state = 0
        nd4.alive = False
        nd4.sampled = False
        nd4.is_sa = False
        nd4.is_sa_dummy_parent = False
        nd4.is_sa_lineage = False

        nd3.add_child(nd4)

        # left extinct child of nd4
        extinct_sp1 = Node(
            taxon=Taxon(label="extinct_sp1"),
            label="extinct_sp1",
            edge_length=1.0
        )
        extinct_sp1.state = 0
        extinct_sp1.alive = False
        extinct_sp1.sampled = False
        extinct_sp1.is_sa = False
        extinct_sp1.is_sa_dummy_parent = False
        extinct_sp1.is_sa_lineage = False

        # right extinct child of nd4
        extinct_sp2 = Node(
            taxon=Taxon(label="extinct_sp2"),
            label="extinct_sp2",
            edge_length=1.0
        )
        extinct_sp2.state = 0
        extinct_sp2.alive = False
        extinct_sp2.sampled = False
        extinct_sp2.is_sa = False
        extinct_sp2.is_sa_dummy_parent = False
        extinct_sp2.is_sa_lineage = False

        nd4.add_child(extinct_sp1)
        nd4.add_child(extinct_sp2)

        # living sampled child of nd3
        nd10 = Node(
            taxon=Taxon(label="nd10"),
            label="nd10",
            edge_length=4.0
        )
        nd10.state = 0
        nd10.alive = True
        nd10.sampled = True
        nd10.is_sa = False
        nd10.is_sa_dummy_parent = False
        nd10.is_sa_lineage = False

        nd3.add_child(nd10)

        # lower child of root
        nd5 = Node(
            taxon=Taxon(label="nd5"),
            label="nd5",
            edge_length=1.0
        )
        nd5.state = 0
        nd5.alive = False
        nd5.sampled = False
        nd5.is_sa = False
        nd5.is_sa_dummy_parent = False
        nd5.is_sa_lineage = False

        root_node.add_child(nd5)

        # upper child of nd5
        nd11 = Node(
            taxon=Taxon(label="nd11"),
            label="nd11",
            edge_length=4.0
        )
        nd11.state = 0
        nd11.alive = True
        nd11.sampled = False
        nd11.is_sa = False
        nd11.is_sa_dummy_parent = False
        nd11.is_sa_lineage = False

        nd5.add_child(nd11)

        # lower child of nd5
        nd6 = Node(
            taxon=Taxon(label="nd6"),
            label="nd6",
            edge_length=1.0
        )
        nd6.state = 0
        nd6.alive = False
        nd6.sampled = False
        nd6.is_sa = False
        nd6.is_sa_dummy_parent = True
        nd6.is_sa_lineage = False

        nd5.add_child(nd6)

        # sampled ancestor child of nd6
        extinct_sp3 = Node(
            taxon=Taxon(label="extinct_sp3"),
            label="extinct_sp3",
            edge_length=0.0
        )
        extinct_sp3.state = 0
        extinct_sp3.alive = False
        extinct_sp3.sampled = True
        extinct_sp3.is_sa = True
        extinct_sp3.is_sa_dummy_parent = False
        extinct_sp3.is_sa_lineage = False

        # lower child of nd6
        nd7 = Node(
            taxon=Taxon(label="nd7"),
            label="nd7",
            edge_length=1.0
        )
        nd7.state = 0
        nd7.alive = False
        nd7.sampled = False
        nd7.is_sa = False
        nd7.is_sa_dummy_parent = False
        nd7.is_sa_lineage = True

        nd6.add_child(extinct_sp3)
        nd6.add_child(nd7)

        # upper living sampled child of nd7
        nd12 = Node(
            taxon=Taxon(label="nd12"),
            label="nd12",
            edge_length=2.0
        )
        nd12.state = 0
        nd12.alive = True
        nd12.sampled = True
        nd12.is_sa = False
        nd12.is_sa_dummy_parent = False
        nd12.is_sa_lineage = False

        # lower child of nd7
        nd8 = Node(
            taxon=Taxon(label="nd8"),
            label="nd8",
            edge_length=1.0
        )
        nd8.state = 0
        nd8.alive = False
        nd8.sampled = False
        nd8.is_sa = False
        nd8.is_sa_dummy_parent = False
        nd8.is_sa_lineage = False

        nd7.add_child(nd12)
        nd7.add_child(nd8)

        # upper extinct child of nd8
        extinct_sp4 = Node(
            taxon=Taxon(label="extinct_sp4"),
            label="extinct_sp4",
            edge_length=0.5
        )
        extinct_sp4.state = 0
        extinct_sp4.alive = False
        extinct_sp4.sampled = False
        extinct_sp4.is_sa = False
        extinct_sp4.is_sa_dummy_parent = False
        extinct_sp4.is_sa_lineage = False

        # lower living sampled child of nd8
        nd13 = Node(
            taxon=Taxon(label="nd13"),
            label="nd13",
            edge_length=1.0
        )
        nd13.state = 0
        nd13.alive = True
        nd13.sampled = True
        nd13.is_sa = False
        nd13.is_sa_dummy_parent = False
        nd13.is_sa_lineage = False

        nd8.add_child(extinct_sp4)
        nd8.add_child(nd13)

        tr_complete = Tree(seed_node=origin_node)

        tr_complete.taxon_namespace.add_taxon(origin_node.taxon)
        tr_complete.taxon_namespace.add_taxon(root_node.taxon)
        tr_complete.taxon_namespace.add_taxon(nd3.taxon)
        tr_complete.taxon_namespace.add_taxon(nd4.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp2.taxon)
        tr_complete.taxon_namespace.add_taxon(nd10.taxon)
        tr_complete.taxon_namespace.add_taxon(nd5.taxon)
        tr_complete.taxon_namespace.add_taxon(nd11.taxon)
        tr_complete.taxon_namespace.add_taxon(nd6.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp3.taxon)
        tr_complete.taxon_namespace.add_taxon(nd7.taxon)
        tr_complete.taxon_namespace.add_taxon(nd12.taxon)
        tr_complete.taxon_namespace.add_taxon(nd8.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp4.taxon)
        tr_complete.taxon_namespace.add_taxon(nd13.taxon)

        total_state_count = 1

        max_age = 6.0

        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            max_age=max_age,
            epsilon=1e-12
        )

        # debugging: checking complete tree is OK
        # print(ann_tr.tree.as_string(schema="newick"))

        tr_rec = ann_tr.extract_reconstructed_tree(
            plotting_overhead=False,
            require_obs_both_sides=False
        )

        ann_tr.populate_n_stubs_on_branches_of_rec_tr()

        # 'nd6' is the dummy parent of the sampled ancestor 'extinct_sp3'.
        # 'nd5' was unifurcation-suppressed when the reconstructed tree was
        # extracted (its other child, 'nd11', was alive but unsampled), so
        # the branch leading to the dummy node 'nd6' absorbs that stub.
        # The stub stays before the fossil, on nd6's incoming edge, not nd7's.
        expected_n_stub_by_label = {
            "root": 0,
            "nd10": 1,
            "nd6": 1,
            "extinct_sp3": 0,
            "nd7": 0,
            "nd12": 0,
            "nd13": 1,
        }

        for nd in tr_rec.preorder_node_iter():
            if nd.label in expected_n_stub_by_label:
                self.assertEqual(
                    nd.n_stub,
                    expected_n_stub_by_label[nd.label],
                    msg=f"n_stub mismatch on node '{nd.label}'"
                )

        # The continuation gets no parent-event spike: its parent marks sampling.
        rec_nd7 = tr_rec.find_node_with_label("nd7")
        self.assertTrue(rec_nd7.is_sa_lineage)

        clock_rate = 0.01
        spike_mean = 0.01
        spike_shape = 2.0
        seed = 123

        # the gradual, relaxed-clock component is always supplied
        # externally now (never drawn internally) -- a fixed, known
        # list stands in for whatever a 'lognormal(n_nr_from=...)' DAG
        # node would have produced, sized to match the tree's own
        # real-branch count exactly (checked below, in preorder, the
        # same order 'n_rec_branches()' counts them in)
        n_real_branches = ann_tr.n_rec_branches()
        br_specific_rel_rates = [
            0.8, 1.3, 0.95, 1.1, 1.05][:n_real_branches]

        # 'populate_gamma_spike_mol_rates_on_rec_tr_branches' takes no
        # 'seed' argument: it draws from NumPy's global legacy random
        # number generator, the same one PhyloJunction seeds once,
        # DAG-wide, from a '.pj' script's '-r' CLI flag or the GUI's
        # seed field (see 'pgm.DirectedAcyclicGraph.random_seed'). Here,
        # outside of that DAG machinery, we seed it ourselves. Only the
        # spike component is drawn from it now, since
        # 'br_specific_rel_rates' (above) is fixed, not RNG-driven.
        np.random.seed(seed)

        pjgsc.populate_gamma_spike_mol_rates_on_rec_tr_branches(
            ann_tr,
            clock_rate=clock_rate,
            spike_mean=spike_mean,
            spike_shape=spike_shape,
            br_specific_rel_rates=br_specific_rel_rates
        )

        # Independently replay the exact same (seeded) draw sequence the
        # function under test should have made, using the number-of-spikes
        # rule described in Douglas et al. (2025) and implemented in the
        # GammaSpikeModel BEAST2 package's 'BranchSpikePrior.getNSpikes()':
        # a branch gets one gamma-distributed spike per stub, plus one
        # more for its own (observed) bifurcation event -- UNLESS the
        # branch is the surviving sibling lineage of a sampled ancestor
        # ('is_sa_lineage'), in which case there is no such event and the
        # "+1" is dropped. Re-seeding with the same seed and replaying
        # through the same global generator, in the same preorder
        # traversal order, keeps the replay in lockstep with the call
        # above, so any deviation from this rule (e.g. always adding
        # "+1") will show up as a mismatched draw as soon as it happens,
        # not just a wrong final 'mol_rate'.
        np.random.seed(seed)

        n_checked_branches = 0
        branch_idx = 0
        for nd in tr_rec.preorder_node_iter():
            if nd.parent_node is None or \
                    nd.edge_length is None or \
                    nd.edge_length <= ann_tr.epsilon:
                self.assertIsNone(nd.branch_raw_rate)
                self.assertIsNone(nd.spike_effect)
                self.assertIsNone(nd.mol_rate)
                continue

            expected_raw_rate = br_specific_rel_rates[branch_idx]
            branch_idx += 1

            n_spikes = nd.n_stub if nd.is_sa_lineage else nd.n_stub + 1

            if n_spikes == 0:
                expected_spike_effect = 0.0
            else:
                expected_spike_effect = spike_mean * np.random.gamma(
                    shape=spike_shape * n_spikes,
                    scale=1.0 / spike_shape
                )

            self.assertEqual(
                nd.branch_raw_rate,
                expected_raw_rate,
                msg=f"branch_raw_rate mismatch on node '{nd.label}'"
            )
            self.assertEqual(
                nd.spike_effect,
                expected_spike_effect,
                msg=f"spike_effect mismatch on node '{nd.label}' "
                    f"(n_stub={nd.n_stub}, is_sa_lineage={nd.is_sa_lineage})"
            )

            self.assertGreater(nd.branch_raw_rate, 0.0)
            self.assertGreaterEqual(nd.spike_effect, 0.0)
            self.assertGreater(nd.mol_rate, 0.0)

            expected_mol_rate = (
                    clock_rate * nd.branch_raw_rate
                    + nd.spike_effect / nd.edge_length
            )

            self.assertAlmostEqual(
                nd.mol_rate,
                expected_mol_rate,
                places=12
            )

            n_checked_branches += 1

        # make sure the loop above did not silently skip everything
        # (e.g. because the reconstructed tree came back empty)
        self.assertEqual(n_checked_branches, 5)

    def test_gamma_spike_zero_spikes_on_sa_lineage_branch(self):
        """The sibling of an SA with zero stubs gets a fixed spike of 0.0.

        When a branch is the surviving lineage sibling of a sampled
        ancestor (is_sa_lineage == True) and has no stubs of its own,
        the number of gamma-distributed spikes contributing to it is
        zero (see test_gamma_spike for the full rule). A Gamma
        distribution with shape 0.0 is a degenerate point mass at 0.0,
        so this case must be special-cased rather than drawn from
        rng.gamma(), which is what this test checks.
        """

        origin_node = Node(
            taxon=Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node.state = 0
        origin_node.alive = False
        origin_node.sampled = False
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        root_node = Node(
            taxon=Taxon(label="root"), label="root", edge_length=1.0)
        root_node.state = 0
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        origin_node.add_child(root_node)

        # other (unrelated) sampled lineage, so root is a true bifurcation
        other_leaf = Node(
            taxon=Taxon(label="other_leaf"),
            label="other_leaf",
            edge_length=2.0
        )
        other_leaf.state = 0
        other_leaf.alive = True
        other_leaf.sampled = True
        other_leaf.is_sa = False
        other_leaf.is_sa_dummy_parent = False
        other_leaf.is_sa_lineage = False

        root_node.add_child(other_leaf)

        # dummy node marking where the sampled ancestor was collected;
        # directly attached to the root, so it carries zero stubs of its
        # own (no unsampled lineage gets pruned/suppressed above it)
        dummy_node = Node(
            taxon=Taxon(label="dummy"), label="dummy", edge_length=1.0)
        dummy_node.state = 0
        dummy_node.alive = False
        dummy_node.sampled = False
        dummy_node.is_sa = False
        dummy_node.is_sa_dummy_parent = True
        dummy_node.is_sa_lineage = False

        root_node.add_child(dummy_node)

        sa_leaf = Node(
            taxon=Taxon(label="sa_leaf"), label="sa_leaf", edge_length=0.0)
        sa_leaf.state = 0
        sa_leaf.alive = False
        sa_leaf.sampled = True
        sa_leaf.is_sa = True
        sa_leaf.is_sa_dummy_parent = False
        sa_leaf.is_sa_lineage = False

        # terminal branch that has a direct ancestor halfway
        lineage_leaf = Node(
            taxon=Taxon(label="lineage_leaf"),
            label="lineage_leaf",
            edge_length=1.0
        )
        lineage_leaf.state = 0
        lineage_leaf.alive = True
        lineage_leaf.sampled = True
        lineage_leaf.is_sa = False
        lineage_leaf.is_sa_dummy_parent = False
        lineage_leaf.is_sa_lineage = True

        dummy_node.add_child(sa_leaf)
        dummy_node.add_child(lineage_leaf)

        tr_complete = Tree(seed_node=origin_node)
        for nd in (origin_node, root_node, other_leaf, dummy_node,
                   sa_leaf, lineage_leaf):
            tr_complete.taxon_namespace.add_taxon(nd.taxon)

        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            trait=DiscreteTrait(1),
            start_at_origin=True,
            max_age=3.0,
            epsilon=1e-12
        )

        # debugging
        # print(ann_tr.tree.as_string(schema="newick"))

        tr_rec = ann_tr.extract_reconstructed_tree(
            plotting_overhead=False,
            require_obs_both_sides=False
        )

        ann_tr.populate_n_stubs_on_branches_of_rec_tr()

        rec_lineage_leaf = tr_rec.find_node_with_label("lineage_leaf")
        self.assertTrue(rec_lineage_leaf.is_sa_lineage)
        self.assertEqual(rec_lineage_leaf.n_stub, 0)

        # no 'seed' argument -- see 'test_gamma_spike' for why; seed
        # the global generator ourselves instead
        np.random.seed(123)

        # the gradual component is always supplied externally now --
        # see 'test_gamma_spike'; the exact values are irrelevant to
        # this test's own assertions (below), only the count matters
        br_specific_rel_rates = [1.0] * ann_tr.n_rec_branches()

        pjgsc.populate_gamma_spike_mol_rates_on_rec_tr_branches(
            ann_tr,
            clock_rate=0.01,
            spike_mean=0.01,
            spike_shape=2.0,
            br_specific_rel_rates=br_specific_rel_rates
        )

        self.assertEqual(rec_lineage_leaf.spike_effect, 0.0)
        self.assertEqual(
            rec_lineage_leaf.mol_rate,
            0.01 * rec_lineage_leaf.branch_raw_rate
        )

    def test_strict_clock_when_no_stubs(self):
        """'strict_clock_when_no_stubs=True' forces a fixed rate when
        a branch has no stubs, and leaves branches that do have a
        stub with their regular relaxed-clock-plus-spike draw.

        This is a deliberate departure from Douglas et al. (2025),
        used only for illustration (see 'examples/
        fbd_read_gamma_spike.pj'), so it is opt-in and
        defaults to 'False' (see 'test_gamma_spike' and
        'test_gamma_spike_zero_spikes_on_sa_lineage_branch', neither
        of which pass this argument).
        """

        origin_node = Node(
            taxon=Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node.state = 0
        origin_node.alive = False
        origin_node.sampled = False
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        root_node = Node(
            taxon=Taxon(label="root"), label="root", edge_length=1.0)
        root_node.state = 0
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        origin_node.add_child(root_node)

        # sampled outgroup leaf, stub-free: keeps 'root' a true
        # bifurcation (otherwise it, too, would be unifurcation-
        # suppressed away during reconstruction)
        leaf_b_out = Node(
            taxon=Taxon(label="leaf_b_out"),
            label="leaf_b_out",
            edge_length=3.0
        )
        leaf_b_out.state = 0
        leaf_b_out.alive = True
        leaf_b_out.sampled = True
        leaf_b_out.is_sa = False
        leaf_b_out.is_sa_dummy_parent = False
        leaf_b_out.is_sa_lineage = False

        # unsampled internal node: its one sampled child inherits a
        # stub once 'internal' is unifurcation-suppressed away during
        # reconstruction
        internal = Node(
            taxon=Taxon(label="internal"), label="internal", edge_length=1.0)
        internal.state = 0
        internal.alive = False
        internal.sampled = False
        internal.is_sa = False
        internal.is_sa_dummy_parent = False
        internal.is_sa_lineage = False

        root_node.add_child(leaf_b_out)
        root_node.add_child(internal)

        leaf_a_sampled = Node(
            taxon=Taxon(label="leaf_a_sampled"),
            label="leaf_a_sampled",
            edge_length=2.0
        )
        leaf_a_sampled.state = 0
        leaf_a_sampled.alive = True
        leaf_a_sampled.sampled = True
        leaf_a_sampled.is_sa = False
        leaf_a_sampled.is_sa_dummy_parent = False
        leaf_a_sampled.is_sa_lineage = False

        leaf_unsampled = Node(
            taxon=Taxon(label="leaf_unsampled"),
            label="leaf_unsampled",
            edge_length=2.0
        )
        leaf_unsampled.state = 0
        leaf_unsampled.alive = True
        leaf_unsampled.sampled = False
        leaf_unsampled.is_sa = False
        leaf_unsampled.is_sa_dummy_parent = False
        leaf_unsampled.is_sa_lineage = False

        internal.add_child(leaf_a_sampled)
        internal.add_child(leaf_unsampled)

        tr_complete = Tree(seed_node=origin_node)
        for nd in (origin_node, root_node, leaf_b_out, internal,
                   leaf_a_sampled, leaf_unsampled):
            tr_complete.taxon_namespace.add_taxon(nd.taxon)

        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            trait=DiscreteTrait(1),
            start_at_origin=True,
            max_age=4.0,
            epsilon=1e-12
        )

        tr_rec = ann_tr.extract_reconstructed_tree(
            plotting_overhead=False,
            require_obs_both_sides=False
        )

        ann_tr.populate_n_stubs_on_branches_of_rec_tr()

        rec_leaf_b_out = tr_rec.find_node_with_label("leaf_b_out")
        rec_leaf_a_sampled = tr_rec.find_node_with_label("leaf_a_sampled")

        # sanity check on the fixture itself: if 'leaf_a_sampled' ever
        # stops inheriting a stub (e.g. the suppression logic
        # changes), this test silently stops covering the
        # has-a-stub branch below
        self.assertEqual(rec_leaf_b_out.n_stub, 0)
        self.assertEqual(rec_leaf_a_sampled.n_stub, 1)

        clock_rate = 0.01

        np.random.seed(123)

        br_specific_rel_rates = [1.0] * ann_tr.n_rec_branches()

        pjgsc.populate_gamma_spike_mol_rates_on_rec_tr_branches(
            ann_tr,
            clock_rate=clock_rate,
            spike_mean=0.3,
            spike_shape=2.0,
            br_specific_rel_rates=br_specific_rel_rates,
            strict_clock_when_no_stubs=True
        )

        # stub-free branch: forced to the exact strict-clock rate, no
        # relaxed-clock or spike randomness at all
        self.assertEqual(rec_leaf_b_out.branch_raw_rate, 1.0)
        self.assertEqual(rec_leaf_b_out.spike_effect, 0.0)
        self.assertEqual(rec_leaf_b_out.mol_rate, clock_rate)

        # branch with a stub: still gets its own independent
        # relaxed-clock-plus-spike draw, so it is not pinned to
        # 'clock_rate' like the branch above
        self.assertGreater(rec_leaf_a_sampled.branch_raw_rate, 0.0)
        self.assertGreater(rec_leaf_a_sampled.spike_effect, 0.0)
        self.assertNotEqual(rec_leaf_a_sampled.mol_rate, clock_rate)

    def test_origin_with_chained_direct_ancestors_before_root(self):
        """A maximally gnarly tree: origin, chained DAs, stubs, all of it.

        PhyloJunction's reconstructed trees can have a chain of direct
        ancestors between the true origin and the first real
        bifurcation ("root") -- e.g., if the origin itself (or a node
        shortly after it) is chosen for ancestor sampling, possibly
        more than once in a row before anything actually speciates.
        This was never exercised by any existing fixture (every prior
        one has the origin leading directly to a 2-child root, with no
        sampled ancestors above it), so this builds one, combined with
        every other source of complexity a reconstructed tree can
        have, to pin the whole machinery down at once.

        Note that, apart from 'root' itself (a label 'AnnotatedTree'
        requires -- see below), none of the hand-picked node labels
        describe position in a drawing (e.g. "below", "above",
        "upper", "lower") -- a tree can be drawn with the root at the
        top, the bottom, or rotated, so that kind of name is ambiguous
        about what it actually means. Labels instead describe a
        node's *role* (its relationship to the sampled-ancestor/stub
        machinery), which is unambiguous regardless of how the tree is
        drawn.

        - origin -> dummy_1 (is_sa_dummy_parent) -> {DA_1, stub_branch}
        - stub_branch is an ORDINARY, unflagged node -- its own
          "unsampled_side" child goes extinct/unsampled, leaving a
          genuine stub that must survive being unifurcation-suppressed
          away. Its surviving child is 'node_X', which is BOTH
          'is_sa_lineage' (continuing DA_1's lineage) AND
          'is_sa_dummy_parent' (marking DA_2) -- a single node playing
          two roles, matching exactly how PhyloJunction's own
          simulator represents this (a node's edge is extended in
          place until the next real event, rather than a new,
          separate node ever being created for "nothing happens
          here" -- see 'distribution.dn_discrete_sse.py's ancestor-
          sampling code).
        - node_X -> {DA_2, root}. The node labeled 'root' here is
          'AnnotatedTree's own required marker for "the first
          bifurcation after the origin, in the complete tree" --
          'AnnotatedTree.__init__' literally scans for a node labeled
          "root" (or 'alternative_root_label') to compute the origin
          edge length, so this label can't be renamed away; it is
          required terminology, not just a description picked for this
          test. Crucially, THIS IS A DIFFERENT THING from "the
          reconstructed tree's own structural root/seed node"
          (tr_rec.seed_node, the node with no parent) -- here, that
          seed node is 'dummy_1' (see the seed-node assertion below),
          several branches above 'root'. The two notions of "root"
          coincide in every simpler fixture (hence earlier blanket
          statements like "the root always has 0 stubs"), but they are
          not the same thing in general: whenever direct ancestors
          chain above the first bifurcation, like here, 'AnnotatedTree'
          's "root" becomes just an ordinary interior node of the
          reconstructed tree, with a real incoming branch -- so it CAN
          could carry a stub on its own path (see
          'expected_n_stub_by_label' below), the same as any other
          real branch. It is only the reconstructed tree's actual seed
          node that is guaranteed to always have zero stubs and no
          rate (see 'tree.is_real_branch'), not every node named
          "root". 'root' is realistically ALSO 'is_sa_lineage' here
          (it IS, structurally, DA_2's surviving lineage; nothing else
          happens in between), so it must correctly get no baseline
          "+1" spike on its own incoming branch, the same as node_X,
          despite 'root' itself genuinely bifurcating (see
          'calculation.gamma_spike_clock_model's module docstring/
          comments for why this is correct, not a bug).
        - 'stub_branch_2': a SECOND, independent stub (via another
          ordinary, unflagged node playing the same role as
          'stub_branch', whose own child goes extinct), feeding into a
          SECOND, ordinary (non-chained, not touching the origin)
          sampled ancestor ('SA_3'/'dummy_2'/'lineage_3') -- to
          confirm the origin-region machinery and the ordinary,
          everywhere-else machinery agree with each other and don't
          interfere.
        - 3 asymmetric extant tips ('tip_a', 'tip_c', 'tip_d').

        Per the project's own convention for this kind of fixture, no
        randomness is involved in what's checked here: 'n_stub' is
        pure topology bookkeeping, and the "is there a baseline spike
        or not" question only needs checking for being exactly 0.0
        (when there is provably no randomness at all, i.e.
        n_spikes == 0) versus greater than 0.0 (when a random spike
        was genuinely drawn) -- never the exact drawn value itself.
        """

        def new_node(label, edge_length, alive, sampled, is_sa=False,
                     is_sa_dummy_parent=False, is_sa_lineage=False):
            nd = Node(
                taxon=Taxon(label=label), label=label,
                edge_length=edge_length)
            nd.state = 0
            nd.alive = alive
            nd.sampled = sampled
            nd.is_sa = is_sa
            nd.is_sa_dummy_parent = is_sa_dummy_parent
            nd.is_sa_lineage = is_sa_lineage
            return nd

        origin_node = new_node("origin", 0.0, False, False)

        dummy_1 = new_node(
            "dummy_1", 1.0, False, False, is_sa_dummy_parent=True)
        origin_node.add_child(dummy_1)

        DA_1 = new_node("DA_1", 0.0, False, True, is_sa=True)
        # ordinary, unflagged node -- the genuine stub lives here
        stub_branch = new_node("stub_branch", 1.0, False, False)
        dummy_1.add_child(DA_1)
        dummy_1.add_child(stub_branch)

        # extinct/unsampled: 'stub_branch' loses this child during
        # reconstruction, leaving behind a stub on whatever survives
        unsampled_side = new_node("unsampled_side", 0.5, False, False)
        # both is_sa_lineage (continuing DA_1) AND is_sa_dummy_parent
        # (marking DA_2) -- the same node plays both roles, as a real
        # simulated tree would represent it
        node_X = new_node(
            "node_X", 1.0, False, False,
            is_sa_lineage=True, is_sa_dummy_parent=True)
        stub_branch.add_child(unsampled_side)
        stub_branch.add_child(node_X)

        DA_2 = new_node("DA_2", 0.0, False, True, is_sa=True)
        # 'root': required label (see docstring above) for "the first
        # true bifurcation" -- realistically ALSO is_sa_lineage, since
        # nothing else happens between DA_2 and this node's own
        # speciation event
        root_node = new_node(
            "root", 2.0, False, False, is_sa_lineage=True)
        node_X.add_child(DA_2)
        node_X.add_child(root_node)

        # tip_a's edge must bring it to the same final depth as
        # tip_c/tip_d below (root's own depth + 3.7 = 8.7), since
        # every extant tip must share the same (maximal) age
        tip_a = new_node("tip_a", 3.7, True, True)
        root_node.add_child(tip_a)

        # second, independent stub, unrelated to the origin region --
        # plays the same role as 'stub_branch' above
        stub_branch_2 = new_node("stub_branch_2", 1.0, False, False)
        root_node.add_child(stub_branch_2)

        extinct_side = new_node("extinct_side", 0.3, False, False)
        dummy_2 = new_node(
            "dummy_2", 0.7, False, False, is_sa_dummy_parent=True)
        stub_branch_2.add_child(extinct_side)
        stub_branch_2.add_child(dummy_2)

        SA_3 = new_node("SA_3", 0.0, False, True, is_sa=True)
        lineage_3 = new_node(
            "lineage_3", 1.0, False, False, is_sa_lineage=True)
        dummy_2.add_child(SA_3)
        dummy_2.add_child(lineage_3)

        tip_c = new_node("tip_c", 1.0, True, True)
        tip_d = new_node("tip_d", 1.0, True, True)
        lineage_3.add_child(tip_c)
        lineage_3.add_child(tip_d)

        all_nodes = (
            origin_node, dummy_1, DA_1, stub_branch, unsampled_side,
            node_X, DA_2, root_node, tip_a, stub_branch_2,
            extinct_side, dummy_2, SA_3, lineage_3, tip_c, tip_d)

        tr_complete = Tree(seed_node=origin_node)
        for nd in all_nodes:
            tr_complete.taxon_namespace.add_taxon(nd.taxon)

        ann_tr = pjtr.AnnotatedTree(
            tr_complete, trait=DiscreteTrait(1), start_at_origin=True,
            max_age=8.7, epsilon=1e-12)

        tr_rec = ann_tr.extract_reconstructed_tree(
            plotting_overhead=False, require_obs_both_sides=False)

        # sanity check on the fixture itself: 'origin' is
        # unifurcation-suppressed away regardless (as always), and
        # both 'stub_branch' and 'stub_branch_2' -- the two ordinary,
        # unflagged nodes carrying the genuine stubs -- must also be
        # suppressed away (one child each was pruned); if this ever
        # stops being true, the rest of this test is not actually
        # exercising suppression-survival at all
        rec_labels = {
            pjtr.get_node_name(nd) for nd in tr_rec.preorder_node_iter()}
        self.assertNotIn("origin", rec_labels)
        self.assertNotIn("stub_branch", rec_labels)
        self.assertNotIn("stub_branch_2", rec_labels)
        self.assertEqual(tr_rec.seed_node.label, "dummy_1")

        # node_X and lineage_3 must have survived reconstruction with
        # their 'is_sa_lineage' flag intact, despite each being the
        # surviving child of a now-suppressed node -- this is the
        # specific thing that would silently break if a stub were
        # ever (incorrectly) placed directly on an 'is_sa_lineage'
        # node instead of on an ordinary one
        rec_node_X = tr_rec.find_node_with_label("node_X")
        rec_root = tr_rec.find_node_with_label("root")
        rec_lineage_3 = tr_rec.find_node_with_label("lineage_3")
        self.assertTrue(rec_node_X.is_sa_lineage)
        self.assertTrue(rec_root.is_sa_lineage)
        self.assertTrue(rec_lineage_3.is_sa_lineage)

        ann_tr.populate_n_stubs_on_branches_of_rec_tr()

        expected_n_stub_by_label = {
            "dummy_1": 0,
            "DA_1": 0,
            "node_X": 1,
            "DA_2": 0,
            "root": 0,
            "tip_a": 0,
            "dummy_2": 1,
            "SA_3": 0,
            "lineage_3": 0,
            "tip_c": 0,
            "tip_d": 0,
        }

        for nd in tr_rec.preorder_node_iter():
            self.assertEqual(
                nd.n_stub,
                expected_n_stub_by_label[pjtr.get_node_name(nd)],
                msg=f"n_stub mismatch on node '{nd.label}'")

        # 'is_real_branch'/'n_rec_branches' agree on exactly the 7
        # branches that have a real parent and non-zero length
        self.assertEqual(ann_tr.n_rec_branches(), 7)
        real_branch_labels = {
            "node_X", "root", "tip_a", "dummy_2", "lineage_3",
            "tip_c", "tip_d"}
        for nd in tr_rec.preorder_node_iter():
            self.assertEqual(
                pjtr.is_real_branch(nd, ann_tr.epsilon),
                pjtr.get_node_name(nd) in real_branch_labels,
                msg=f"is_real_branch mismatch on node '{nd.label}'")

        np.random.seed(123)

        br_specific_rel_rates = [1.0] * ann_tr.n_rec_branches()

        pjgsc.populate_gamma_spike_mol_rates_on_rec_tr_branches(
            ann_tr,
            clock_rate=0.01,
            spike_mean=0.1,
            spike_shape=2.0,
            br_specific_rel_rates=br_specific_rel_rates
        )

        # not real branches: no rate of their own, full stop
        for label in ("dummy_1", "DA_1", "DA_2", "SA_3"):
            nd = tr_rec.find_node_with_label(label)
            self.assertIsNone(nd.branch_raw_rate)
            self.assertIsNone(nd.spike_effect)
            self.assertIsNone(nd.mol_rate)

        # These continuations have neither stubs on their own segments nor a
        # speciation at their parent. Earlier stubs remain on node_X and dummy_2.
        for node in (rec_root, rec_lineage_3):
            self.assertEqual(node.spike_effect, 0.0)
            self.assertGreater(node.branch_raw_rate, 0.0)

        for label in ("node_X", "tip_a", "dummy_2", "tip_c", "tip_d"):
            nd = tr_rec.find_node_with_label(label)
            self.assertGreater(
                nd.branch_raw_rate, 0.0,
                msg=f"branch_raw_rate should be positive on '{label}'")
            self.assertGreater(
                nd.spike_effect, 0.0,
                msg=f"spike_effect should be strictly positive on "
                    f"'{label}'")


if __name__ == "__main__":
    unittest.main()
