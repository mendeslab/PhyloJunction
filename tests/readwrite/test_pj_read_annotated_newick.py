import phylojunction.calculation.gamma_spike_clock_model as pjgsc
import unittest
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# pj imports
import phylojunction.readwrite.pj_read as pjr
import phylojunction.data.tree as pjtr
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestReadAnnotatedNewick(unittest.TestCase):
    """Test 'read_nwk_tree_str's 'is_complete'/'keep_unsampled_alive'.

    These let a Newick string be read either as an already-
    reconstructed tree (the default, 'is_complete=False' -- every
    leaf is sampled by construction, no annotations needed), or as a
    *complete* tree (e.g. exported from a simulation,
    'is_complete=True') whose leaves must be explicitly annotated
    with '&alive='/'&sampled=' NHX comments, since branch lengths
    alone cannot be trusted to tell a fossil tip (dead, but sampled)
    from an unsampled extinction event for arbitrary external data.
    """

    # Child order must not change which edge continues after an observation; existing
    # examples use only one order. Retain while sampled ancestors use dummy parents.
    def test_continuing_lineage_flags_ignore_child_order(self):
        for children in ("S:0,A:2", "A:2,S:0", "S:0,(T:0,A:1)D:1"):
            with self.subTest(children=children):
                tree = pjr.read_nwk_tree_str("(" + children + ");", in_file=False).tree
                for node in tree.preorder_node_iter():
                    if node.parent_node is not None:
                        self.assertEqual(node.is_sa_lineage,
                                         node.parent_node.is_sa_dummy_parent and not node.is_sa)

    def test_complete_tree_with_fossil_tip_is_reconstructed_correctly(self):
        """A fossil (dead, sampled) tip survives reconstruction.

        'tipa' is an extant sample, 'sa1'/'lineage' is a direct
        ancestor and its surviving lineage, and 'fossiltip' is a
        genuine fossil: dead (not at maximal depth), but explicitly
        annotated 'sampled=True'. All three should survive
        'extract_reconstructed_tree' -- 'fossiltip' is exactly the
        case this feature exists for (see 'tree.AnnotatedTree.
        extract_reconstructed_tree's relaxed 'filter_fn', which no
        longer requires 'alive' for a sampled leaf to be retained).
        """

        nwk = ("(tipa[&alive=True]:3.0,"
               "(sa1:0.0,lineage[&alive=True]:2.0):1.0,"
               "fossiltip[&alive=False,sampled=True]:1.5):0.0;")

        ann_tr = pjr.read_nwk_tree_str(nwk, in_file=False, is_complete=True)

        tipa = ann_tr.tree.find_node_with_taxon_label("tipa")
        sa1 = ann_tr.tree.find_node_with_taxon_label("sa1")
        lineage = ann_tr.tree.find_node_with_taxon_label("lineage")
        fossiltip = ann_tr.tree.find_node_with_taxon_label("fossiltip")

        self.assertTrue(tipa.alive)
        self.assertTrue(tipa.sampled)
        self.assertFalse(tipa.is_sa)

        self.assertTrue(sa1.is_sa)
        self.assertFalse(sa1.alive)
        self.assertTrue(sa1.sampled)

        self.assertTrue(lineage.alive)
        self.assertTrue(lineage.sampled)

        self.assertFalse(fossiltip.alive)
        self.assertTrue(fossiltip.sampled)
        self.assertFalse(fossiltip.is_sa)

        tr_rec = ann_tr.extract_reconstructed_tree(
            plotting_overhead=False, require_obs_both_sides=False)
        rec_leaf_labels = {
            nd.taxon.label for nd in tr_rec.leaf_node_iter()}

        self.assertEqual(
            rec_leaf_labels, {"tipa", "sa1", "lineage", "fossiltip"})

    def test_complete_tree_fossil_only_is_not_discarded(self):
        """A tree with only fossil (no extant, no SA) samples survives.

        Exercises 'extract_reconstructed_tree's early-return guard,
        which must count 'n_extinct_sampled_terminal_nodes' (newly
        populated by '_count_terminal_nodes') -- not just extant
        samples and direct ancestors -- or a fossil-only tree would
        be incorrectly treated as having nothing to reconstruct.
        """

        nwk = ("(fossila[&alive=False,sampled=True]:1.0,"
               "fossilb[&alive=False,sampled=True]:2.0):0.0;")

        ann_tr = pjr.read_nwk_tree_str(nwk, in_file=False, is_complete=True)

        tr_rec = ann_tr.extract_reconstructed_tree(
            plotting_overhead=False, require_obs_both_sides=False)

        self.assertTrue(ann_tr.tree_died)
        self.assertEqual(ann_tr.n_extant_terminal_nodes, 0)
        self.assertEqual(ann_tr.n_extinct_terminal_nodes, 2)
        self.assertEqual(ann_tr.n_extinct_sampled_terminal_nodes, 2)
        self.assertGreater(len(tr_rec.leaf_nodes()), 0)
        rec_leaf_labels = {
            nd.taxon.label for nd in tr_rec.leaf_node_iter()}
        self.assertEqual(rec_leaf_labels, {"fossila", "fossilb"})

    # Depth-based tests cannot protect status when ages are rounded. Check explicit flags
    # and their serialization together; relevant while complete Newick uses these annotations.
    def test_near_present_annotations_survive_roundtrip(self):
        nwk = ("(A[&alive=true,sampled=true]:1,"
               "B[&alive=true,sampled=false]:0.9999,"
               "F[&alive=false,sampled=true]:1,"
               "U[&alive=false,sampled=false]:0.9999);")
        for _ in range(2):
            tree = pjr.read_nwk_tree_str(nwk, in_file=False, is_complete=True)
            self.assertFalse(tree.tree_died)
            self.assertEqual(set(tree.extant_terminal_nodes_labels), {"A", "B"})
            self.assertEqual(set(tree.extant_sampled_terminal_nodes_labels), {"A"})
            self.assertEqual(set(tree.extinct_terminal_nodes_labels), {"F", "U"})
            self.assertEqual(set(tree.extinct_sampled_terminal_nodes_labels), {"F"})
            rec = tree.extract_reconstructed_tree(plotting_overhead=False)
            self.assertEqual({nd.taxon.label for nd in rec.leaf_node_iter()}, {"A", "F"})
            nwk = tree.tree.as_string(schema="newick", suppress_annotations=False)

    def test_complete_tree_gamma_spike_end_to_end(self):
        """Gamma-spike runs cleanly on a complete tree with a fossil tip.

        No specific rate values are checked (those are covered
        elsewhere, e.g. 'test_gamma_spike_clock_model.py') -- this
        just confirms the whole pipeline (stub-counting, then rate
        population) does not choke on a fossil-tip-bearing tree read
        from an annotated Newick string, and that the fossil tip gets
        a real rate of its own, same as any other real branch.
        """

        nwk = ("(tipa[&alive=True]:3.0,"
               "(sa1:0.0,lineage[&alive=True]:2.0):1.0,"
               "fossiltip[&alive=False,sampled=True]:1.5):0.0;")

        ann_tr = pjr.read_nwk_tree_str(nwk, in_file=False, is_complete=True)

        np.random.seed(42)
        br_specific_rel_rates = [1.0] * ann_tr.n_rec_branches()
        pjgsc.populate_gamma_spike_mol_rates_on_rec_tr_branches(
            ann_tr, clock_rate=0.01, spike_mean=0.1,
            spike_shape=2.0,
            br_specific_rel_rates=br_specific_rel_rates)

        fossiltip = ann_tr.tree_reconstructed.find_node_with_taxon_label(
            "fossiltip")
        self.assertIsNotNone(fossiltip.mol_rate)
        self.assertGreater(fossiltip.mol_rate, 0.0)

    def test_node_names_attribute_root_label_with_no_origin(self):
        """'alternative_root_label' is the true root when there is no origin.

        'node_names_attribute' (used to cross-reference externally-
        computed stochastic maps by node index) relabels every node,
        including the root, by its index-derived name (e.g. "nd0")
        instead of the literal string "root". When the tree has no
        separate origin node (the seed node IS the root, i.e. it has
        two children), the seed node is never compared against the
        literal "root" string, so 'seen_root' used to stay False for
        the rest of parsing -- the *next* non-dummy-parent internal
        node visited (here, the root's own left child) would then be
        wrongly adopted as "the root" instead, and every later such
        node would keep overwriting 'alternative_root_label' again,
        leaving it pointing at whichever one was visited last ("nd3"
        here) rather than at the true root ("nd0"). This, in turn,
        broke 'extract_reconstructed_tree's fallback search for the
        root when 'map_attr(fun="smap", ...)' is combined with
        'gamma_spike_clock' downstream.
        """

        nwk = ("((tipa:1.0[&index=1],tipb:1.0[&index=2]):0.5[&index=3],"
               "tipc:1.5[&index=4]):0.0[&index=0];")

        ann_tr = pjr.read_nwk_tree_str(
            nwk, in_file=False, node_names_attribute="index")

        self.assertEqual(ann_tr.tree.seed_node.label, "nd0")
        self.assertEqual(ann_tr.alternative_root_label, "nd0")

        tr_rec = ann_tr.extract_reconstructed_tree(
            plotting_overhead=False, require_obs_both_sides=False)

        self.assertEqual(tr_rec.seed_node.label, "nd0")
        rec_leaf_labels = {nd.label for nd in tr_rec.leaf_node_iter()}
        self.assertEqual(rec_leaf_labels, {"nd1", "nd2", "nd4"})

    def test_missing_alive_annotation_raises(self):
        """A dead leaf with no '&alive=' annotation is rejected.

        'is_complete=True' means branch lengths/depth are not trusted
        at all -- there is no silent fallback here (unlike the
        'is_complete=False' default).
        """

        nwk = "(tipa[&alive=True]:2.0,fossiltip[&sampled=True]:1.0):0.0;"

        with self.assertRaises(
                ec.AnnotatedTreeMissingNodeAnnotationError) as exc:
            pjr.read_nwk_tree_str(nwk, in_file=False, is_complete=True)

        self.assertIn("fossiltip", str(exc.exception))
        self.assertIn("alive", str(exc.exception))

    def test_missing_sampled_annotation_on_dead_leaf_raises(self):
        """A leaf explicitly marked dead, but with no '&sampled=', errors.

        'alive=False' alone is not enough -- PhyloJunction cannot
        tell a fossil tip from an unsampled extinction event without
        an explicit '&sampled=' annotation too.
        """

        nwk = "(tipa[&alive=True]:2.0,deadtip[&alive=False]:1.0):0.0;"

        with self.assertRaises(
                ec.AnnotatedTreeMissingNodeAnnotationError) as exc:
            pjr.read_nwk_tree_str(nwk, in_file=False, is_complete=True)

        self.assertIn("deadtip", str(exc.exception))
        self.assertIn("sampled", str(exc.exception))

    def test_internal_node_annotated_alive_raises(self):
        """An internal node explicitly annotated '&alive=True' errors.

        Internal nodes are never alive nor sampled -- this is a
        structural invariant, not something inferred from the
        annotation, so an explicit, contradictory annotation on one
        must be rejected rather than silently ignored.
        """

        nwk = ("(tipa[&alive=True]:2.0,"
               "(tipb[&alive=True]:1.0,tipc[&alive=True]:1.0)"
               "[&alive=True]:1.0):0.0;")

        with self.assertRaises(
                ec.AnnotatedTreeIncorrectAnnotationError) as exc:
            pjr.read_nwk_tree_str(nwk, in_file=False, is_complete=True)

        self.assertIn("internal nodes are never alive", str(exc.exception))

    def test_internal_node_annotated_sampled_raises(self):
        """Symmetric check: an internal node '&sampled=True' errors too."""

        nwk = ("(tipa[&alive=True]:2.0,"
               "(tipb[&alive=True]:1.0,tipc[&alive=True]:1.0)"
               "[&sampled=True]:1.0):0.0;")

        with self.assertRaises(
                ec.AnnotatedTreeIncorrectAnnotationError) as exc:
            pjr.read_nwk_tree_str(nwk, in_file=False, is_complete=True)

        self.assertIn(
            "internal nodes are never sampled", str(exc.exception))

    def test_direct_ancestor_annotated_alive_raises(self):
        """A direct ancestor (zero-length leaf) '&alive=True' errors.

        Direct (sampled) ancestors are mid-lineage sampling events,
        never present-day tips -- this holds structurally, regardless
        of any annotation saying otherwise.
        """

        nwk = "(tipa[&alive=True]:2.0,sa1[&alive=True]:0.0):0.0;"

        with self.assertRaises(
                ec.AnnotatedTreeIncorrectAnnotationError) as exc:
            pjr.read_nwk_tree_str(nwk, in_file=False, is_complete=True)

        self.assertIn(
            "direct (sampled) ancestors are never alive",
            str(exc.exception))

    def test_direct_ancestor_annotated_unsampled_raises(self):
        """A direct ancestor '&sampled=False' errors too.

        Direct (sampled) ancestors are always sampled, by definition
        -- an explicit contradiction is rejected the same way as for
        'alive' above.
        """

        nwk = "(tipa[&alive=True]:2.0,sa1[&sampled=False]:0.0):0.0;"

        with self.assertRaises(
                ec.AnnotatedTreeIncorrectAnnotationError) as exc:
            pjr.read_nwk_tree_str(nwk, in_file=False, is_complete=True)

        self.assertIn(
            "direct (sampled) ancestors are always sampled",
            str(exc.exception))

    def test_keep_unsampled_alive_false_honors_annotation(self):
        """Default: an alive leaf's explicit '&sampled=False' is honored.

        'keep_unsampled_alive=False' (the default) treats such a leaf
        the same way an unsampled extinction would be -- it is pruned
        away when the reconstructed tree is extracted.
        """

        nwk = ("(tipa[&alive=True,sampled=False]:2.0,"
               "tipb[&alive=True]:2.0):0.0;")

        ann_tr = pjr.read_nwk_tree_str(
            nwk, in_file=False, is_complete=True,
            keep_unsampled_alive=False)

        tipa = ann_tr.tree.find_node_with_taxon_label("tipa")
        self.assertFalse(tipa.sampled)

        tr_rec = ann_tr.extract_reconstructed_tree(
            plotting_overhead=False, require_obs_both_sides=False)
        rec_leaf_labels = {
            nd.taxon.label for nd in tr_rec.leaf_node_iter()}
        self.assertNotIn("tipa", rec_leaf_labels)
        self.assertIn("tipb", rec_leaf_labels)

    def test_keep_unsampled_alive_true_overrides_annotation(self):
        """'keep_unsampled_alive=True' keeps an unsampled living leaf.

        Represents a known-but-unsequenced living taxon: despite its
        explicit '&sampled=False' annotation, it is forced back to
        'sampled=True' and so survives reconstruction, rather than
        being pruned away like an unsampled extinction would be.
        """

        nwk = ("(tipa[&alive=True,sampled=False]:2.0,"
               "tipb[&alive=True]:2.0):0.0;")

        ann_tr = pjr.read_nwk_tree_str(
            nwk, in_file=False, is_complete=True,
            keep_unsampled_alive=True)

        tipa = ann_tr.tree.find_node_with_taxon_label("tipa")
        self.assertTrue(tipa.sampled)

        tr_rec = ann_tr.extract_reconstructed_tree(
            plotting_overhead=False, require_obs_both_sides=False)
        rec_leaf_labels = {
            nd.taxon.label for nd in tr_rec.leaf_node_iter()}
        self.assertIn("tipa", rec_leaf_labels)

    def test_keep_unsampled_alive_true_no_effect_when_absent(self):
        """'keep_unsampled_alive=True' with no annotation still samples.

        The flag only changes behavior when a contrary ('&sampled=
        False') annotation is actually given -- an alive leaf with no
        '&sampled=' at all defaults to sampled regardless of the flag.
        """

        nwk = "(tipa[&alive=True]:2.0,tipb[&alive=True]:2.0):0.0;"

        ann_tr = pjr.read_nwk_tree_str(
            nwk, in_file=False, is_complete=True,
            keep_unsampled_alive=True)

        tipa = ann_tr.tree.find_node_with_taxon_label("tipa")
        self.assertTrue(tipa.sampled)

    def test_is_complete_default_treats_dead_leaf_as_sampled(self):
        """Default ('is_complete=False'): a dead leaf is still sampled.

        The default means "this Newick string is already the
        reconstructed tree" -- every leaf is sampled by construction,
        whether it is alive or a fossil, unlike the old (pre-
        'is_complete') heuristic, which only ever forced 'sampled=True'
        for leaves reaching the tree's maximum depth. No annotations
        are needed or consulted in this mode.
        """

        nwk = "(tipa:2.0,deadtip:1.0):0.0;"

        ann_tr = pjr.read_nwk_tree_str(nwk, in_file=False)

        tipa = ann_tr.tree.find_node_with_taxon_label("tipa")
        deadtip = ann_tr.tree.find_node_with_taxon_label("deadtip")

        self.assertTrue(tipa.alive)
        self.assertTrue(tipa.sampled)

        self.assertFalse(deadtip.alive)
        self.assertTrue(deadtip.sampled)

        tr_rec = ann_tr.extract_reconstructed_tree(
            plotting_overhead=False, require_obs_both_sides=False)
        rec_leaf_labels = {
            nd.taxon.label for nd in tr_rec.leaf_node_iter()}
        self.assertEqual(rec_leaf_labels, {"tipa", "deadtip"})

    def test_is_complete_default_ignores_internal_node_alive_annotation(self):
        """Default mode still never checks an *internal* node's '&alive='.

        Unlike an ordinary leaf's '&alive=' (now rejected outright --
        see 'test_is_complete_default_rejects_alive_annotation'
        below), an internal node's own annotations are not consulted
        at all in this mode, so one that would be invalid under
        'is_complete=True' (e.g. '&alive=True' on an internal node)
        does not raise here.
        """

        nwk = "(tipa:2.0,(tipb:1.0,tipc:1.0)[&alive=True]:1.0):0.0;"

        # must not raise
        ann_tr = pjr.read_nwk_tree_str(nwk, in_file=False)

        tipb = ann_tr.tree.find_node_with_taxon_label("tipb")
        self.assertTrue(tipb.sampled)

    def test_is_complete_default_rejects_sampled_false_on_leaf(self):
        """An ordinary leaf's explicit '&sampled=False' is rejected.

        'is_complete=False' (the default) means the Newick string is
        already the *reconstructed* tree, so every leaf is sampled by
        construction -- a leaf explicitly annotated '&sampled=False'
        contradicts that, and is rejected rather than silently
        ignored. A fossil tip needs no annotation at all here -- it
        is detected purely from depth, same as any other leaf; see
        'test_is_complete_default_treats_dead_leaf_as_sampled' above,
        and 'test_is_complete_default_rejects_alive_annotation' below
        for the symmetric '&alive=' check.
        """

        nwk = "(tipa:2.0,deadtip[&sampled=False]:1.0):0.0;"

        with self.assertRaises(
                ec.AnnotatedTreeIncorrectAnnotationError) as exc:
            pjr.read_nwk_tree_str(nwk, in_file=False)

        self.assertIn(
            "every leaf must be sampled, by definition",
            str(exc.exception))

    def test_is_complete_default_rejects_alive_annotation(self):
        """An ordinary leaf's explicit '&alive=' is rejected outright.

        Under 'is_complete=False' (the default), 'alive' is always
        deduced from depth -- an '&alive=' annotation is never
        consulted at all, so one being present is rejected outright,
        regardless of whether its value would have agreed with the
        depth-based inference or not (unlike '&sampled=', where only
        the one value that can never be reconciled -- 'False' -- is
        rejected; see
        'test_is_complete_default_rejects_sampled_false_on_leaf'
        above). Both exist to catch a script that most likely meant
        to pass 'is_complete="true"' instead.
        """

        nwk = "(tipa:2.0,deadtip[&alive=False]:1.0):0.0;"

        with self.assertRaises(
                ec.AnnotatedTreeIncorrectAnnotationError) as exc:
            pjr.read_nwk_tree_str(nwk, in_file=False)

        self.assertIn(
            "never consults it -- 'alive' is always deduced from "
            "depth instead", str(exc.exception))

        # even a value that *would* have agreed with depth-based
        # inference (here, 'tipa' genuinely is at maximal depth) is
        # still rejected outright
        nwk_matching = "(tipa[&alive=True]:2.0,tipb:2.0):0.0;"

        with self.assertRaises(
                ec.AnnotatedTreeIncorrectAnnotationError):
            pjr.read_nwk_tree_str(nwk_matching, in_file=False)

    def test_is_complete_default_rejects_sampled_false_on_direct_ancestor(self):
        """A direct ancestor's explicit '&sampled=False' is rejected too.

        Same reasoning as the ordinary-leaf case above, just for a
        zero-length (direct ancestor) leaf instead -- it is always
        sampled, by definition, under 'is_complete=False' as well.
        """

        nwk = "(tipa:2.0,sa1[&sampled=False]:0.0):0.0;"

        with self.assertRaises(
                ec.AnnotatedTreeIncorrectAnnotationError) as exc:
            pjr.read_nwk_tree_str(nwk, in_file=False)

        self.assertIn(
            "every leaf -- including a direct (sampled) ancestor -- "
            "must be sampled, by definition",
            str(exc.exception))

    def test_seed_node_with_direct_ancestor_child_is_dummy_parent(self):
        """A direct ancestor branching off right at the seed is detected.

        The seed node (origin or root) used to always get
        'is_sa_dummy_parent=False' hardcoded, unlike every other
        internal node (which scans its own children for a zero-length
        one). A direct ancestor can perfectly well branch off right
        at the seed itself (e.g. a fossil sampled very early on, with
        the rest of the tree as its only sibling) -- when that
        happens, the seed ends up with exactly one non-SA child, which
        'tree.get_y_coord_from_n_obs_nodes' has no branch for unless
        the seed is correctly recognized as a dummy parent (it is
        neither the origin-average case, nor a two-real-children
        bifurcating node) -- silently leaving it out of the returned
        y-coordinates dict entirely, surfacing downstream as a
        'KeyError' in 'tree._draw_clade' when plotting the complete
        tree. This checks the fix directly (the attribute itself),
        and then that plotting no longer raises.
        """

        nwk = "(sa1:0.0,(tipa:1.0,tipb:1.0):1.0):0.0;"

        ann_tr = pjr.read_nwk_tree_str(nwk, in_file=False)

        self.assertTrue(ann_tr.tree.seed_node.is_sa_dummy_parent)

        fig = plt.figure()
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])

        # must not raise
        pjtr.plot_ann_tree(ann_tr, axes, draw_reconstructed=False)


if __name__ == "__main__":
    unittest.main()
