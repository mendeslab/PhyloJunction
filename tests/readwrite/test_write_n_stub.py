import unittest
import pandas as pd
from dendropy import Node, Taxon, Tree

# pj imports
import phylojunction.data.tree as pjtr
import phylojunction.readwrite.pj_write as pjw

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestWriteNStub(unittest.TestCase):

    @staticmethod
    def _new_node(label: str, edge_length: float, alive: bool,
                  sampled: bool) -> Node:
        nd = Node(taxon=Taxon(label=label), label=label,
                  edge_length=edge_length)
        nd.state = 0
        nd.alive = alive
        nd.sampled = sampled
        nd.is_sa = False
        nd.is_sa_dummy_parent = False
        nd.is_sa_lineage = False
        return nd

    def _tree_with_one_stub(self) -> pjtr.AnnotatedTree:
        """Two-tip reconstructed tree where one tip inherits a stub.

        Same topology used in 'test_gamma_spike_clock_model.
        test_strict_clock_when_no_stubs': an unsampled sibling lineage
        is unifurcation-suppressed away during reconstruction, and the
        stub left behind is transferred onto 'leaf_a_sampled'.
        """

        origin_node = self._new_node("origin", 0.0, False, False)
        root_node = self._new_node("root", 1.0, False, False)
        origin_node.add_child(root_node)

        leaf_b_out = self._new_node("leaf_b_out", 3.0, True, True)
        internal = self._new_node("internal", 1.0, False, False)
        root_node.add_child(leaf_b_out)
        root_node.add_child(internal)

        leaf_a_sampled = self._new_node("leaf_a_sampled", 2.0, True, True)
        leaf_unsampled = self._new_node("leaf_unsampled", 2.0, True, False)
        internal.add_child(leaf_a_sampled)
        internal.add_child(leaf_unsampled)

        tr_complete = Tree(seed_node=origin_node)
        for nd in (origin_node, root_node, leaf_b_out, internal,
                   leaf_a_sampled, leaf_unsampled):
            tr_complete.taxon_namespace.add_taxon(nd.taxon)

        return pjtr.AnnotatedTree(
            tr_complete, total_state_count=1, start_at_origin=True,
            max_age=4.0, epsilon=1e-12)

    def _tree_with_no_stubs(self) -> pjtr.AnnotatedTree:
        """A plain, fully-sampled two-tip cherry: no stubs anywhere."""

        origin_node = self._new_node("origin", 0.0, False, False)
        root_node = self._new_node("root", 1.0, False, False)
        origin_node.add_child(root_node)

        leaf_x = self._new_node("leaf_x", 1.0, True, True)
        leaf_y = self._new_node("leaf_y", 1.0, True, True)
        root_node.add_child(leaf_x)
        root_node.add_child(leaf_y)

        tr_complete = Tree(seed_node=origin_node)
        for nd in (origin_node, root_node, leaf_x, leaf_y):
            tr_complete.taxon_namespace.add_taxon(nd.taxon)

        return pjtr.AnnotatedTree(
            tr_complete, total_state_count=1, start_at_origin=True,
            max_age=2.0, epsilon=1e-12)

    def test_prep_n_stub_df(self) -> None:
        """One row per tree, one column per node label seen in ANY tree."""

        tree1 = self._tree_with_one_stub()
        tree2 = self._tree_with_no_stubs()

        n_stub_df = pjw.prep_n_stub_df(
            [tree1, tree2], sample_size=2, n_repl=1)

        self.assertEqual(list(n_stub_df["sample"]), [1, 2])
        self.assertEqual(list(n_stub_df["replicate"]), [1, 1])

        # column set is the UNION of both trees' node labels, even
        # though they do not overlap at all (different topologies)
        for label in ("leaf_b_out", "leaf_a_sampled", "leaf_x", "leaf_y"):
            self.assertIn(label, n_stub_df.columns)

        # tree 1 (row 0): 'leaf_a_sampled' has the one real stub,
        # 'leaf_b_out' has none
        self.assertEqual(n_stub_df.at[0, "leaf_a_sampled"], 1)
        self.assertEqual(n_stub_df.at[0, "leaf_b_out"], 0)

        # tree 1 has no node labeled 'leaf_x' -- that cell is empty,
        # not (incorrectly) 0
        self.assertTrue(pd.isna(n_stub_df.at[0, "leaf_x"]))

        # tree 2 (row 1): no stubs anywhere
        self.assertEqual(n_stub_df.at[1, "leaf_x"], 0)
        self.assertEqual(n_stub_df.at[1, "leaf_y"], 0)
        self.assertTrue(pd.isna(n_stub_df.at[1, "leaf_a_sampled"]))

    def test_prep_n_stub_df_n_repl_greater_than_one(self) -> None:
        """2 samples x 2 replicates: rows/'sample'/'replicate' line up.

        'node_val' is expected in (sample 1 repl 1, sample 1 repl 2,
        sample 2 repl 1, sample 2 repl 2, ...) order -- the same
        convention every other per-(sample, replicate) table in
        PhyloJunction uses (e.g. 'initialize_tree_dataframe', which
        'prep_n_stub_df' itself builds its "sample"/"replicate"
        columns from). This pins that down explicitly for
        'prep_n_stub_df', alternating between the two fixture trees
        so each row's stub counts are also independently checkable.
        """

        trees = [
            self._tree_with_one_stub(),   # sample 1, repl 1
            self._tree_with_no_stubs(),   # sample 1, repl 2
            self._tree_with_no_stubs(),   # sample 2, repl 1
            self._tree_with_one_stub(),   # sample 2, repl 2
        ]

        n_stub_df = pjw.prep_n_stub_df(trees, sample_size=2, n_repl=2)

        self.assertEqual(list(n_stub_df["sample"]), [1, 1, 2, 2])
        self.assertEqual(list(n_stub_df["replicate"]), [1, 2, 1, 2])

        # row 0 (sample 1, repl 1): the one-stub tree
        self.assertEqual(n_stub_df.at[0, "leaf_a_sampled"], 1)
        self.assertEqual(n_stub_df.at[0, "leaf_b_out"], 0)
        self.assertTrue(pd.isna(n_stub_df.at[0, "leaf_x"]))

        # row 1 (sample 1, repl 2): the no-stubs tree
        self.assertEqual(n_stub_df.at[1, "leaf_x"], 0)
        self.assertEqual(n_stub_df.at[1, "leaf_y"], 0)
        self.assertTrue(pd.isna(n_stub_df.at[1, "leaf_a_sampled"]))

        # row 2 (sample 2, repl 1): the no-stubs tree again
        self.assertEqual(n_stub_df.at[2, "leaf_x"], 0)
        self.assertEqual(n_stub_df.at[2, "leaf_y"], 0)
        self.assertTrue(pd.isna(n_stub_df.at[2, "leaf_a_sampled"]))

        # row 3 (sample 2, repl 2): the one-stub tree again
        self.assertEqual(n_stub_df.at[3, "leaf_a_sampled"], 1)
        self.assertEqual(n_stub_df.at[3, "leaf_b_out"], 0)
        self.assertTrue(pd.isna(n_stub_df.at[3, "leaf_x"]))

    def test_prep_n_stub_df_all_zero_stubs_has_no_na(self) -> None:
        """When every tree shares the same (stub-free) topology, no NA.

        Complements 'test_prep_n_stub_df', which only checks that a
        column gets NA where a node genuinely does not exist in a
        given tree. This checks the opposite, equally important case:
        when every tree in the batch shares the exact same topology
        (same node labels throughout), nothing should ever be NA --
        every cell should be a real, explicit stub count (0, here).
        """

        trees = [self._tree_with_no_stubs() for _ in range(3)]

        n_stub_df = pjw.prep_n_stub_df(trees, sample_size=3, n_repl=1)

        value_columns = [
            col for col in n_stub_df.columns
            if col not in ("sample", "replicate")]

        self.assertEqual(set(value_columns), {"root", "leaf_x", "leaf_y"})
        self.assertFalse(n_stub_df[value_columns].isna().any().any())
        self.assertTrue((n_stub_df[value_columns] == 0).all().all())


if __name__ == "__main__":
    unittest.main()
