from phylojunction.data.trait import DiscreteTrait
import unittest
import copy
import matplotlib
from dendropy import Tree, Node, Taxon

# pj imports
import phylojunction.data.tree as pjtr
import phylojunction.data.sampled_ancestor as pjsa
import phylojunction.data.attribute_transition as pjat

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestExtractReconstructedTree(unittest.TestCase):

    def test_obs_both_sides_root_with_origin_true(self):
        """
        Test method .is_extant_or_sa_on_both_sides_complete_tr_root().
        
        Tests if AnnotatedTree's method verifies that complete tree
        has one observed taxon on both sides of its root.

        To see it on icytree:
        ((sp2:1.0,(sp3:0.5,sp4:0.5)nd1:0.5)root:1.0)origin:0.0;
        """

        origin_node = Node(taxon=Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node.state = 0
        origin_node.alive = False
        origin_node.sampled = False
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        root_node = Node(taxon=Taxon(label="root"), label="root", edge_length=1.0)
        root_node.state = 0
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        origin_node.add_child(root_node)

        # left child of root node
        extant_sp1 = Node(taxon=Taxon(label="sp2"), label="sp2", edge_length=1.0)
        extant_sp1.state = 0
        extant_sp1.alive = True
        extant_sp1.sampled = True
        extant_sp1.is_sa = False
        extant_sp1.is_sa_dummy_parent = False
        extant_sp1.is_sa_lineage = False

        # right child of root_node
        internal_node1 = Node(taxon=Taxon(label="nd1"), label="nd1", edge_length=0.5)
        internal_node1.state = 0
        internal_node1.alive = False
        internal_node1.sampled = False
        internal_node1.is_sa = False
        internal_node1.is_sa_dummy_parent = False
        internal_node1.is_sa_lineage = False

        root_node.add_child(extant_sp1)
        root_node.add_child(internal_node1)

        # left child of nd2
        extant_sp2 = Node(taxon=Taxon(label="sp3"), label="sp3", edge_length=0.5)
        extant_sp2.state = 0
        extant_sp2.alive = True
        extant_sp2.sampled = True
        extant_sp2.is_sa = False
        extant_sp2.is_sa_dummy_parent = False
        extant_sp2.is_sa_lineage = False

        # right child of nd2
        extant_sp3 = Node(taxon=Taxon(label="sp4"), label="sp4", edge_length=0.5)
        extant_sp3.state = 0
        extant_sp3.alive = True
        extant_sp3.sampled = True
        extant_sp3.is_sa = False
        extant_sp3.is_sa_dummy_parent = False
        extant_sp3.is_sa_lineage = False

        internal_node1.add_child(extant_sp2)
        internal_node1.add_child(extant_sp3)

        tr_complete = Tree(seed_node=origin_node)
        tr_complete.taxon_namespace.add_taxon(origin_node.taxon)
        tr_complete.taxon_namespace.add_taxon(root_node.taxon)
        tr_complete.taxon_namespace.add_taxon(internal_node1)
        tr_complete.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp2.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp3.taxon)

        total_state_count = 1
        
        max_age = 2.0

        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            max_age=max_age,
            epsilon=1e-12)

        # print(ann_tr.tree.as_string(
        #     schema="newick",
        #     suppress_annotations=True,
        #     suppress_internal_taxon_labels=True))

        self.assertTrue(ann_tr.is_extant_or_sa_on_both_sides_complete_tr_root(root_node))

    def test_obs_both_sides_root_with_origin_true_larger(self):
        """
        Test method .is_extant_or_sa_on_both_sides_complete_tr_root().
        
        Tests if AnnotatedTree's method verifies that complete tree
        has one observed taxon on both sides of its root.

        To see it on icytree:
        (((sp1:0.5,sp2:0.5)nd1:1.0,(sp3:0.4,(sp4:0.25,sp5:0.25)nd3:0.25)nd2:1.0)root:1.0)origin:0.0;
        """

        origin_node = Node(taxon=Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node.state = 0
        origin_node.alive = False
        origin_node.sampled = False
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        root_node = Node(taxon=Taxon(label="root"), label="root", edge_length=1.0)
        root_node.state = 0
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        origin_node.add_child(root_node)

        # left child of root
        internal_node1 = Node(taxon=Taxon(label="nd1"), label="nd1", edge_length=1.0)
        internal_node1.state = 0
        internal_node1.alive = False
        internal_node1.sampled = False
        internal_node1.is_sa = False
        internal_node1.is_sa_dummy_parent = False
        internal_node1.is_sa_lineage = False

        # left child of internal_node1
        extant_sp1 = Node(taxon=Taxon(label="sp1"), label="sp1", edge_length=0.5)
        extant_sp1.state = 0
        extant_sp1.alive = True
        extant_sp1.sampled = True
        extant_sp1.is_sa = False
        extant_sp1.is_sa_dummy_parent = False
        extant_sp1.is_sa_lineage = False

        # right child of internal_node1
        extant_sp2 = Node(taxon=Taxon(label="sp2"), label="sp2", edge_length=0.5)
        extant_sp2.state = 0
        extant_sp2.alive = False
        extant_sp2.sampled = False
        extant_sp2.is_sa = False
        extant_sp2.is_sa_dummy_parent = False
        extant_sp2.is_sa_lineage = False

        internal_node1.add_child(extant_sp1)
        internal_node1.add_child(extant_sp2)

        # right child of root
        internal_node2 = Node(taxon=Taxon(label="nd2"), label="nd2", edge_length=1.0)
        internal_node2.state = 0
        internal_node2.alive = False
        internal_node2.sampled = False
        internal_node2.is_sa = False
        internal_node2.is_sa_dummy_parent = False
        internal_node2.is_sa_lineage = False

        # left child of internal_node2
        extant_sp3 = Node(taxon=Taxon(label="sp3"), label="sp3", edge_length=0.4)
        extant_sp3.state = 0
        extant_sp3.alive = False
        extant_sp3.sampled = False
        extant_sp3.is_sa = False
        extant_sp3.is_sa_dummy_parent = False
        extant_sp3.is_sa_lineage = False

        # right child of internal_node2
        internal_node3 = Node(taxon=Taxon(label="nd3"), label="nd3", edge_length=0.25)
        internal_node3.state = 0
        internal_node3.alive = False
        internal_node3.sampled = False
        internal_node3.is_sa = False
        internal_node3.is_sa_dummy_parent = False
        internal_node3.is_sa_lineage = False

        # left child of internal_node3
        extant_sp4 = Node(taxon=Taxon(label="sp4"), label="sp4", edge_length=0.25)
        extant_sp4.state = 0
        extant_sp4.alive = True
        extant_sp4.sampled = True
        extant_sp4.is_sa = False
        extant_sp4.is_sa_dummy_parent = False
        extant_sp4.is_sa_lineage = False

        # right child of internal_node3
        extant_sp5 = Node(taxon=Taxon(label="sp5"), label="sp5", edge_length=0.25)
        extant_sp5.state = 0
        extant_sp5.alive = True
        extant_sp5.sampled = True
        extant_sp5.is_sa = False
        extant_sp5.is_sa_dummy_parent = False
        extant_sp5.is_sa_lineage = False

        internal_node3.add_child(extant_sp4)
        internal_node3.add_child(extant_sp5)

        internal_node2.add_child(extant_sp3)
        internal_node2.add_child(internal_node3)

        root_node.add_child(internal_node1)
        root_node.add_child(internal_node2)

        tr_complete = Tree(seed_node=origin_node)
        tr_complete.taxon_namespace.add_taxon(origin_node.taxon)
        tr_complete.taxon_namespace.add_taxon(root_node.taxon)
        tr_complete.taxon_namespace.add_taxon(internal_node1)
        tr_complete.taxon_namespace.add_taxon(internal_node2)
        tr_complete.taxon_namespace.add_taxon(internal_node3)
        tr_complete.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp2.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp3.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp4.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp5.taxon)

        total_state_count = 1
        
        max_age = 2.0

        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            max_age=max_age,
            epsilon=1e-12)

        # print(
        #     ann_tr.tree.as_string(schema="newick",
        #                           suppress_annotations=True,
        #                           suppress_internal_taxon_labels=True))

        self.assertTrue(ann_tr.is_extant_or_sa_on_both_sides_complete_tr_root(root_node))

    def test_obs_both_sides_root_with_origin_false(self):
        """
        Test method .is_extant_or_sa_on_both_sides_complete_tr_root().

        Tests if AnnotatedTree's method verifies that complete tree
        does not have one observed taxon on both sides of its root.

        To see it on icytree:
        ((sp2:1.0,(sp3:0.5,sp4:0.5)nd1:0.5)root:1.0)origin:0.0;
        """

        origin_node = Node(taxon=Taxon(label="origin"),
                           label="origin",
                           edge_length=0.0)
        origin_node.state = 0
        origin_node.alive = False
        origin_node.sampled = False
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        root_node = Node(taxon=Taxon(label="root"),
                         label="root",
                         edge_length=1.0)
        root_node.state = 0
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        origin_node.add_child(root_node)

        # left child of root node
        extant_sp1 = Node(taxon=Taxon(label="sp2"),
                          label="sp2",
                          edge_length=1.0)
        extant_sp1.state = 0
        extant_sp1.alive = False
        extant_sp1.sampled = False
        extant_sp1.is_sa = False
        extant_sp1.is_sa_dummy_parent = False
        extant_sp1.is_sa_lineage = False

        # right child of root_node
        internal_node1 = Node(taxon=Taxon(label="nd1"),
                              label="nd1",
                              edge_length=0.5)
        internal_node1.state = 0
        internal_node1.alive = False
        internal_node1.sampled = False
        internal_node1.is_sa = False
        internal_node1.is_sa_dummy_parent = False
        internal_node1.is_sa_lineage = False

        root_node.add_child(extant_sp1)
        root_node.add_child(internal_node1)

        # left child of nd2
        extant_sp2 = Node(taxon=Taxon(label="sp3"),
                          label="sp3",
                          edge_length=0.5)
        extant_sp2.state = 0
        extant_sp2.alive = True
        extant_sp2.sampled = True
        extant_sp2.is_sa = False
        extant_sp2.is_sa_dummy_parent = False
        extant_sp2.is_sa_lineage = False

        # right child of nd2
        extant_sp3 = Node(taxon=Taxon(label="sp4"),
                          label="sp4",
                          edge_length=0.5)
        extant_sp3.state = 0
        extant_sp3.alive = True
        extant_sp3.sampled = True
        extant_sp3.is_sa = False
        extant_sp3.is_sa_dummy_parent = False
        extant_sp3.is_sa_lineage = False

        internal_node1.add_child(extant_sp2)
        internal_node1.add_child(extant_sp3)

        tr_complete = Tree(seed_node=origin_node)
        tr_complete.taxon_namespace.add_taxon(origin_node.taxon)
        tr_complete.taxon_namespace.add_taxon(root_node.taxon)
        tr_complete.taxon_namespace.add_taxon(internal_node1)
        tr_complete.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp2.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp3.taxon)

        total_state_count = 1
        
        max_age = 2.0

        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            max_age=max_age,
            epsilon=1e-12)

        # print(ann_tr.tree.as_string(
        #     schema="newick",
        #     suppress_annotations=True,
        #     suppress_internal_taxon_labels=True))

        self.assertFalse(ann_tr.\
                         is_extant_or_sa_on_both_sides_complete_tr_root(
                             root_node))

    def test_extract_reconstructed_tree_two_extant_one_extinct_one_sa(self):
        """
        Test extraction of reconstructed tree from complete tree (no origin).

        Tree has two extant taxa (sp1 and sp3), one extinct taxon (sp2), one sampled ancestor.
        There are observable taxa on both sides of the tree, so conditioning on that
        or not is irrelevant.
        The reconstructed-tree string should always have something in it.

        To see it on icytree: ((sa1:0.0,sp1:1.0)dummy1:1.0,(sp2:1.0,sp3:1.25)nd1:0.75)root:0.0;
        """

        root_node = Node(taxon=Taxon(label="root"), label="root", edge_length=0.0)
        root_node.state = 0
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        # left child of root node after ancestor sampling happens on who would have been the left child ("sp1")
        dummy_node = Node(taxon=Taxon(label="dummy1"), label="dummy1", edge_length=1.0)
        dummy_node.state = 0
        dummy_node.alive = False
        dummy_node.sampled = False
        dummy_node.is_sa = False
        dummy_node.is_sa_dummy_parent = True
        dummy_node.is_sa_lineage = False

        # right child of root
        int_node = Node(taxon=Taxon(label="nd1"), label="nd1", edge_length=0.75)
        int_node.state = 0
        int_node.alive = False
        int_node.sampled = False
        int_node.is_sa = False
        int_node.is_sa_dummy_parent = False
        int_node.is_sa_lineage = False

        extinct_sp2 = Node(taxon=Taxon(label="sp2"), label="sp2", edge_length=1.0)
        extinct_sp2.state = 0
        extinct_sp2.alive = False
        extinct_sp2.sampled = False
        extinct_sp2.is_sa = False
        extinct_sp2.is_sa_dummy_parent = False
        extinct_sp2.is_sa_lineage = False

        extant_sp3 = Node(taxon=Taxon(label="sp3"), label="sp3", edge_length=1.25)
        extant_sp3.state = 0
        extant_sp3.alive = True
        extant_sp3.sampled = True
        extant_sp3.is_sa = False
        extant_sp3.is_sa_dummy_parent = False
        extant_sp3.is_sa_lineage = False

        int_node.add_child(extinct_sp2)
        int_node.add_child(extant_sp3)
        
        # right child of dummy_node
        sa_node = Node(taxon=Taxon(label="sa1"), label="sa1", edge_length=0.0)
        sa_node.state = 0
        sa_node.alive = False
        sa_node.sampled = True
        sa_node.is_sa = True
        sa_node.is_sa_dummy_parent = False
        sa_node.is_sa_lineage = False

        # left child of dummy node
        # left child of root node
        extant_sp1 = Node(taxon=Taxon(label="sp1"), label="sp1", edge_length=1.0)
        extant_sp1.state = 0
        extant_sp1.alive = True
        extant_sp1.sampled = True
        extant_sp1.is_sa = False
        extant_sp1.is_sa_dummy_parent = False
        extant_sp1.is_sa_lineage = True
        
        dummy_node.add_child(sa_node)
        dummy_node.add_child(extant_sp1)

        root_node.add_child(dummy_node)
        root_node.add_child(int_node)

        tr_complete = Tree(seed_node=root_node)
        tr_complete.taxon_namespace.add_taxon(root_node.taxon)
        tr_complete.taxon_namespace.add_taxon(dummy_node.taxon)
        tr_complete.taxon_namespace.add_taxon(sa_node.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp2.taxon)

        total_state_count = 1

        sa_global_time = 1.0
        time_to_sa_lineage_node = 1.0
        sa = pjsa.SampledAncestor(
            "sa1",
            "sp1",
            sa_global_time,
            time_to_lineage_node=time_to_sa_lineage_node)
        sa_lineage_dict = { "sp1": [sa] }
        
        max_age = 2.0

        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(total_state_count),
            start_at_origin=False,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        # print(ann_tr.tree.as_string(schema="newick", suppress_annotations=True))

        tr_rec1 = \
            ann_tr.extract_reconstructed_tree(plotting_overhead=False,
                                              require_obs_both_sides=False)
        tr_rec2 = \
            ann_tr.extract_reconstructed_tree(plotting_overhead=False,
                                              require_obs_both_sides=True)

        tr_rec_str1 = \
            tr_rec1.as_string(schema="newick",
                              suppress_internal_taxon_labels=True,
                              suppress_internal_node_labels=True,
                              suppress_rooting=True)
        tr_rec_str2 = \
            tr_rec2.as_string(schema="newick",
                              suppress_internal_taxon_labels=True,
                              suppress_internal_node_labels=True,
                              suppress_rooting=True)

        self.assertEqual(
            ann_tr.tree.as_string(schema="newick", suppress_annotations=True),
            ("((sa1:0.0,sp1:1.0)dummy1_dummy1:1.0,(sp2:1.0,sp3:1.25)nd1_nd1:0.75)"
             "root_root:0.0;\n"))
        self.assertEqual(
            tr_rec_str1,
            "((sa1:0.0,sp1:1.0):1.0,sp3:2.0):0.0;\n")
        self.assertEqual(
            tr_rec_str2,
            "((sa1:0.0,sp1:1.0):1.0,sp3:2.0):0.0;\n")

    def test_extract_reconstructed_tree_one_extant_two_extinct_one_sa(self):
        """
        Test extraction of rec tree from complete tree (no origin).

        Tree has two extant taxa (sp1 and sp3), one extinct taxon (sp2),
        one sampled ancestor.
        
        Without conditioning on having observed taxa on both sides of
        the tree, we get a stich with a bead in the middle.
        
        If conditioning on having observed taxa on both sides of the
        tree, the reconstructed-tree string should be empty.

        To see it on icytree: ((sa1:0.0,sp1:0.75):1.0,(sp2:0.8,sp3:0.8):0.75):0.0;
        """

        root_node = Node(taxon=Taxon(label="root"), label="root", edge_length=0.0)
        root_node.state = 0
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        # left child of root node after ancestor sampling happens on who would have been the left child ("sp1")
        dummy_node = Node(taxon=Taxon(label="dummy1"), label="dummy1", edge_length=1.0)
        dummy_node.state = 0
        dummy_node.alive = False
        dummy_node.sampled = False
        dummy_node.is_sa = False
        dummy_node.is_sa_dummy_parent = True
        dummy_node.is_sa_lineage = False

        # right child of root
        int_node = Node(taxon=Taxon(label="nd1"), label="nd1", edge_length=0.75)
        int_node.state = 0
        int_node.alive = False
        int_node.sampled = False
        int_node.is_sa = False
        int_node.is_sa_dummy_parent = False
        int_node.is_sa_lineage = False

        extinct_sp2 = Node(taxon=Taxon(label="sp2"), label="sp2", edge_length=0.8)
        extinct_sp2.state = 0
        extinct_sp2.alive = False
        extinct_sp2.sampled = False
        extinct_sp2.is_sa = False
        extinct_sp2.is_sa_dummy_parent = False
        extinct_sp2.is_sa_lineage = False

        extant_sp3 = Node(taxon=Taxon(label="sp3"), label="sp3", edge_length=0.8)
        extant_sp3.state = 0
        extant_sp3.alive = False
        extant_sp3.sampled = False
        extant_sp3.is_sa = False
        extant_sp3.is_sa_dummy_parent = False
        extant_sp3.is_sa_lineage = False

        int_node.add_child(extinct_sp2)
        int_node.add_child(extant_sp3)
        
        # right child of dummy_node
        sa_node = Node(taxon=Taxon(label="sa1"), label="sa1", edge_length=0.0)
        sa_node.state = 0
        sa_node.alive = False
        sa_node.sampled = False
        sa_node.is_sa = True
        sa_node.is_sa_dummy_parent = False
        sa_node.is_sa_lineage = False

        # left child of dummy node
        # left child of root node
        extant_sp1 = Node(taxon=Taxon(label="sp1"), label="sp1", edge_length=0.75)
        extant_sp1.state = 0
        extant_sp1.alive = True
        extant_sp1.sampled = True
        extant_sp1.is_sa = False
        extant_sp1.is_sa_dummy_parent = False
        extant_sp1.is_sa_lineage = False
        
        dummy_node.add_child(sa_node)
        dummy_node.add_child(extant_sp1)

        root_node.add_child(dummy_node)
        root_node.add_child(int_node)

        tr_complete = Tree(seed_node=root_node)
        tr_complete.taxon_namespace.add_taxon(root_node.taxon)
        tr_complete.taxon_namespace.add_taxon(dummy_node.taxon)
        tr_complete.taxon_namespace.add_taxon(sa_node.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp2.taxon)

        total_state_count = 1

        sa_global_time = 1.0
        time_to_sa_lineage_node = 1.0
        sa = pjsa.SampledAncestor(
            "sa1",
            "sp1",
            sa_global_time,
            time_to_lineage_node=time_to_sa_lineage_node)
        sa_lineage_dict = { "sp1": [sa] }
        
        max_age = 2.0

        ann_tr1 = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(total_state_count),
            start_at_origin=False,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        ann_tr2 = copy.deepcopy(ann_tr1)

        tr_rec1 = \
            ann_tr1.extract_reconstructed_tree(plotting_overhead=False,
                                               require_obs_both_sides=False)
        tr_rec2 = \
            ann_tr2.extract_reconstructed_tree(plotting_overhead=False,
                                               require_obs_both_sides=True)

        tr_rec_str1 = tr_rec1.as_string(
            schema="newick",
            suppress_internal_taxon_labels=True,
            suppress_internal_node_labels=True,
            suppress_rooting=True)

        tr_rec_str2 = tr_rec2.as_string(
            schema="newick",
            suppress_internal_taxon_labels=True,
            suppress_internal_node_labels=True,
            suppress_rooting=True)

        self.assertEqual(
            ann_tr1.tree.as_string(schema="newick",
                                   suppress_annotations=True),
                                   ("((sa1:0.0,sp1:0.75)dummy1_dummy1:1.0,"
                                    "(sp2:0.8,sp3:0.8)nd1_nd1:0.75)root_root"
                                    ":0.0;\n"))
        self.assertEqual(tr_rec_str1, "(sa1:0.0,sp1:0.75):0.0;\n")
        self.assertEqual(tr_rec_str2, ";\n")

    def test_extract_reconstructed_tree_all_extinct_one_sa(self):
        """
        Test extraction of rec. tree from complete tree (no origin).

        Tree has one observable taxon, on one side of the root, a
        sampled ancestor. Test checks that the reconstructed tree
        is either a single taxon string, or an empty string (when
        taxa are required on both sides of the root).

        To see it on icytree:
        ((sa1:0.0,sp1:0.5)dummy1:1.0,sp2:1.25)root:0.0;
        """

        root_node = Node(taxon=Taxon(label="root"), label="root", edge_length=0.0)
        root_node.state = 0
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        # left child of root node after ancestor sampling happens on who would have been the left child ("sp1")
        dummy_node = Node(taxon=Taxon(label="dummy1"), label="dummy1", edge_length=1.0)
        dummy_node.state = 0
        dummy_node.alive = False
        dummy_node.sampled = False
        dummy_node.is_sa = False
        dummy_node.is_sa_dummy_parent = True
        dummy_node.is_sa_lineage = False

        # right child of root node, goes extinct
        extant_sp2 = Node(taxon=Taxon(label="sp2"), label="sp2", edge_length=1.25)
        extant_sp2.state = 0
        extant_sp2.alive = False
        extant_sp2.sampled = False
        extant_sp2.is_sa = False
        extant_sp2.is_sa_dummy_parent = False
        extant_sp2.is_sa_lineage = False

        root_node.add_child(dummy_node)
        root_node.add_child(extant_sp2)
        
        # right child of dummy_node
        sa_node = Node(taxon=Taxon(label="sa1"), label="sa1", edge_length=0.0)
        sa_node.state = 0
        sa_node.alive = False
        sa_node.sampled = True
        sa_node.is_sa = True
        sa_node.is_sa_dummy_parent = False
        sa_node.is_sa_lineage = False

        # left child of dummy node
        # left child of root node
        extant_sp1 = Node(taxon=Taxon(label="sp1"), label="sp1", edge_length=0.5)
        extant_sp1.state = 0
        extant_sp1.alive = False
        extant_sp1.sampled = False
        extant_sp1.is_sa = False
        extant_sp1.is_sa_dummy_parent = False
        extant_sp1.is_sa_lineage = False
        
        dummy_node.add_child(sa_node)
        dummy_node.add_child(extant_sp1)

        tr_complete = Tree(seed_node=root_node)
        tr_complete.taxon_namespace.add_taxon(root_node.taxon)
        tr_complete.taxon_namespace.add_taxon(dummy_node.taxon)
        tr_complete.taxon_namespace.add_taxon(sa_node.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp2.taxon)

        total_state_count = 1

        sa_global_time = 1.0
        time_to_sa_lineage_node = 1.0
        sa = pjsa.SampledAncestor(
            "sa1",
            "sp1",
            sa_global_time,
            time_to_lineage_node=time_to_sa_lineage_node)
        sa_lineage_dict = { "sp1": [sa] }
        
        max_age = 2.0

        ann_tr1 = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(total_state_count),
            start_at_origin=False,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        ann_tr2 = copy.deepcopy(ann_tr1)

        tr_rec1 = \
            ann_tr1.extract_reconstructed_tree(plotting_overhead=False,
                                               require_obs_both_sides=False)
        tr_rec2 = \
            ann_tr2.extract_reconstructed_tree(plotting_overhead=False,
                                               require_obs_both_sides=True)

        tr_rec_str1 = tr_rec1.as_string(
            schema="newick",
            suppress_internal_taxon_labels=True,
            suppress_internal_node_labels=True)

        # debugging
        # print(tr_rec_str1)

        tr_rec_str2 = tr_rec2.as_string(
            schema="newick",
            suppress_internal_taxon_labels=True,
            suppress_internal_node_labels=True)
        # print(tr_rec_str2)

        self.assertEqual(
            ann_tr1.tree.as_string(
                schema="newick",
                suppress_annotations=True),
                ("((sa1:0.0,sp1:0.5)dummy1_dummy1:1.0,sp2:1.25)"
                 "root_root:0.0;\n"))
        self.assertEqual(tr_rec_str1, "(sa1:1.0):0.0;\n")
        self.assertEqual(tr_rec_str2, ";\n")

    def test_extract_reconstructed_tree_all_extinct_no_sa(self):
        """
        Test extraction of rec. tree from complete tree (no origin).

        Tree has died, string for reconstructed tree should be empty.
        """

        root_node = Node(taxon=Taxon(label="root"),
                         label="root",
                         edge_length=0.0)
        root_node.state = 0
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        # left child of root node
        extant_sp1 = Node(taxon=Taxon(label="sp1"),
                          label="sp1",
                          edge_length=1.0)
        extant_sp1.state = 0
        extant_sp1.alive = False
        extant_sp1.sampled = False
        extant_sp1.is_sa = False
        extant_sp1.is_sa_dummy_parent = False
        extant_sp1.is_sa_lineage = False

        # right child of root node, goes extinct
        extant_sp2 = Node(taxon=Taxon(label="sp2"),
                          label="sp2",
                          edge_length=1.0)
        extant_sp2.state = 0
        extant_sp2.alive = False
        extant_sp2.sampled = False
        extant_sp2.is_sa = False
        extant_sp2.is_sa_dummy_parent = False
        extant_sp2.is_sa_lineage = False

        root_node.add_child(extant_sp1)
        root_node.add_child(extant_sp2)

        tr_complete = Tree(seed_node=root_node)
        tr_complete.taxon_namespace.add_taxon(root_node.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp2.taxon)

        total_state_count = 1

        max_age = 2.0

        ann_tr1 = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(total_state_count),
            start_at_origin=False,
            max_age=max_age,
            epsilon=1e-12)

        ann_tr2 = copy.deepcopy(ann_tr1)

        tr_rec1 = ann_tr1.extract_reconstructed_tree(
            plotting_overhead=False,
            require_obs_both_sides=False)
        tr_rec2 = ann_tr2.extract_reconstructed_tree(
            plotting_overhead=False,
            require_obs_both_sides=True)

        tr_rec_str1 = tr_rec1.as_string(
            schema="newick",
            suppress_internal_taxon_labels=True,
            suppress_internal_node_labels=True)
        tr_rec_str2 = tr_rec2.as_string(
            schema="newick",
            suppress_internal_taxon_labels=True,
            suppress_internal_node_labels=True)
        
        self.assertEqual(
            ann_tr1.tree.as_string(
                schema="newick",
                suppress_annotations=True),
                "(sp1:1.0,sp2:1.0)root_root:0.0;\n")
        self.assertEqual(tr_rec_str1, ";\n")
        self.assertEqual(tr_rec_str2, ";\n")

    def test_extract_reconstructed_tree_one_root_side_dies_three_survive(self):
        """
        Test extraction of rec. tree from complete tree (no origin).

        The complete tree has three surviving taxa on only one side of
        the tree. We re-root the reconstructed tree at the MRCA of
        those three taxa.

        To see it on icytree:
        (sp1:0.8,(sp2:1.0,(sp3:0.5,sp4:0.5)nd2:0.5)nd1:1.0)root:0.0;
        """

        root_node = Node(taxon=Taxon(label="root"),
                         label="root",
                         edge_length=0.0)
        root_node.state = 0
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        # left child of root_node
        extinct_sp1 = Node(taxon=Taxon(label="sp1"),
                           label="sp1",
                           edge_length=0.8)
        extinct_sp1.state = 0
        extinct_sp1.alive = False
        extinct_sp1.sampled = False
        extinct_sp1.is_sa = False
        extinct_sp1.is_sa_dummy_parent = False
        extinct_sp1.is_sa_lineage = False

        # right child of root_node
        internal_node1 = Node(taxon=Taxon(label="nd1"),
                              label="nd1",
                              edge_length=1.0)
        internal_node1.state = 0
        internal_node1.alive = True
        internal_node1.sampled = True
        internal_node1.is_sa = False
        internal_node1.is_sa_dummy_parent = False
        internal_node1.is_sa_lineage = False

        root_node.add_child(extinct_sp1)
        root_node.add_child(internal_node1)

        # left child of nd1
        extant_sp1 = Node(taxon=Taxon(label="sp2"),
                          label="sp2",
                          edge_length=1.0)
        extant_sp1.state = 0
        extant_sp1.alive = True
        extant_sp1.sampled = True
        extant_sp1.is_sa = False
        extant_sp1.is_sa_dummy_parent = False
        extant_sp1.is_sa_lineage = False

        # right child of nd1
        internal_node2 = Node(taxon=Taxon(label="nd2"),
                              label="nd2",
                              edge_length=0.5)
        internal_node2.state = 0
        internal_node2.alive = False
        internal_node2.sampled = False
        internal_node2.is_sa = False
        internal_node2.is_sa_dummy_parent = False
        internal_node2.is_sa_lineage = False

        internal_node1.add_child(extant_sp1)
        internal_node1.add_child(internal_node2)

        # left child of nd2
        extant_sp2 = Node(taxon=Taxon(label="sp3"),
                          label="sp3",
                          edge_length=0.5)
        extant_sp2.state = 0
        extant_sp2.alive = True
        extant_sp2.sampled = True
        extant_sp2.is_sa = False
        extant_sp2.is_sa_dummy_parent = False
        extant_sp2.is_sa_lineage = False

        # right child of nd2
        extant_sp3 = Node(taxon=Taxon(label="sp4"),
                          label="sp4",
                          edge_length=0.5)
        extant_sp3.state = 0
        extant_sp3.alive = True
        extant_sp3.sampled = True
        extant_sp3.is_sa = False
        extant_sp3.is_sa_dummy_parent = False
        extant_sp3.is_sa_lineage = False

        internal_node2.add_child(extant_sp2)
        internal_node2.add_child(extant_sp3)

        tr_complete = Tree(seed_node=root_node)
        tr_complete.taxon_namespace.add_taxon(root_node.taxon)
        tr_complete.taxon_namespace.add_taxon(internal_node1)
        tr_complete.taxon_namespace.add_taxon(internal_node2)
        tr_complete.taxon_namespace.add_taxon(extinct_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp2.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp3.taxon)

        total_state_count = 1
        
        max_age = 2.0

        ann_tr1 = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(total_state_count),
            start_at_origin=False,
            max_age=max_age,
            epsilon=1e-12)

        ann_tr2 = copy.deepcopy(ann_tr1)

        tr_rec1 = ann_tr1.extract_reconstructed_tree(
            plotting_overhead=False,
            require_obs_both_sides=False)
        tr_rec2 = ann_tr2.extract_reconstructed_tree(
            plotting_overhead=False,
            require_obs_both_sides=True)

        tr_rec_str1 = tr_rec1.as_string(
            schema="newick",
            suppress_internal_taxon_labels=True,
            suppress_internal_node_labels=False,
            suppress_rooting=True)
        
        tr_rec_str2 = tr_rec2.as_string(
            schema="newick",
            suppress_internal_taxon_labels=True,
            suppress_internal_node_labels=False)

        self.assertEqual(
            ann_tr1.tree.as_string(
                schema="newick",
                suppress_annotations=True,
                suppress_internal_taxon_labels=True),
                ("(sp1:0.8,(sp2:1.0,(sp3:0.5,sp4:0.5)nd2:0.5)nd1:1.0)"
                 "root:0.0;\n"))
        self.assertEqual(tr_rec_str1,
                         "(sp2:1.0,(sp3:0.5,sp4:0.5)nd2:0.5)nd1:0.0;\n")
        self.assertEqual(tr_rec_str2, ";\n")

    def test_extract_reconstructed_tree_origin_sa_before_root_tree_dies(self):
        """
        Test extraction of rec. tree from complete tree (has origin).

        The complete tree has only one observable taxon (an SA) and then dies.
        Tree has one observable taxon (an SA) before the root, and then
        dies. Test checks that the reconstructed tree is either a single
        taxon string, or an empty string (when taxa are required on both
        sides of the complete tree's root).
        """

        origin_node = Node(taxon=Taxon(label="origin"),
                           label="origin",
                           edge_length=0.0)
        origin_node.state = 0
        origin_node.alive = False
        origin_node.sampled = False
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        # only child of origin node
        dummy_node = Node(taxon=Taxon(label="dummy1"),
                          label="dummy1",
                          edge_length=1.0)
        dummy_node.state = 0
        dummy_node.alive = False
        dummy_node.sampled = False
        dummy_node.is_sa = False
        dummy_node.is_sa_dummy_parent = True
        dummy_node.is_sa_lineage = False

        origin_node.add_child(dummy_node)
        
        # left child of dummy_node
        sa_node = Node(taxon=Taxon(label="sa1"),
                       label="sa1",
                       edge_length=0.0)
        sa_node.state = 0
        sa_node.alive = False
        sa_node.sampled = False
        sa_node.is_sa = True
        sa_node.is_sa_dummy_parent = False
        sa_node.is_sa_lineage = False

        # right child of dummy node, goes extinct
        brosc_node = Node(taxon=Taxon(label="brosc"),
                          label="brosc",
                          edge_length=0.75)
        brosc_node.state = 0
        brosc_node.alive = False
        brosc_node.sampled = False
        brosc_node.is_sa = False
        brosc_node.is_sa_dummy_parent = False
        brosc_node.is_sa_lineage = True

        dummy_node.add_child(sa_node)
        dummy_node.add_child(brosc_node)

        tr_complete = Tree(seed_node=origin_node)
        tr_complete.taxon_namespace.add_taxon(origin_node.taxon)
        tr_complete.taxon_namespace.add_taxon(dummy_node.taxon)
        tr_complete.taxon_namespace.add_taxon(sa_node.taxon)
        tr_complete.taxon_namespace.add_taxon(brosc_node.taxon)

        total_state_count = 1

        sa_global_time = 1.0
        time_to_sa_lineage_node = 1.0
        sa = pjsa.SampledAncestor(
            "sa1",
            "brosc",
            sa_global_time,
            time_to_lineage_node=time_to_sa_lineage_node)
        sa_lineage_dict = { "brosc": [sa] }
        
        max_age = 2.0

        ann_tr1 = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        ann_tr2 = copy.deepcopy(ann_tr1)

        tr_rec1 = ann_tr1.extract_reconstructed_tree(
            plotting_overhead=False,
            require_obs_both_sides=False)
        tr_rec2 = ann_tr2.extract_reconstructed_tree(
            plotting_overhead=False,
            require_obs_both_sides=True)

        tr_rec_str1 = tr_rec1.as_string(
            schema="newick",
            suppress_internal_taxon_labels=True,
            suppress_internal_node_labels=True)
        tr_rec_str2 = tr_rec2.as_string(
            schema="newick",
            suppress_internal_taxon_labels=True,
            suppress_internal_node_labels=True)

        self.assertEqual(
            ann_tr1.tree.as_string(
                schema="newick",
                suppress_annotations=True,
                suppress_internal_taxon_labels=True),
                "((sa1:0.0,brosc:0.75)dummy1:1.0)origin:0.0;\n")
        self.assertEqual(tr_rec_str1, "(sa1:1.0):0.0;\n")
        self.assertEqual(tr_rec_str2, ";\n")

    def test_extract_reconstructed_tree_origin_sa_before_root_tree_one_survives(self):
        """
        Test extraction of rec. tree from complete tree (has origin).

        The complete tree has two observable taxa, an SA, and then ONE
        surviving extant taxon sp1.

        To see it on icytree:
        ((sa1:0.0,brosc:1.0)dummy1_dummy1:1.0)origin_origin:0.0;
        """
        
        origin_node = Node(taxon=Taxon(label="origin"),
                           label="origin",
                           edge_length=0.0)
        origin_node.state = 0
        origin_node.alive = False
        origin_node.sampled = False
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        # only child of origin node
        dummy_node = Node(taxon=Taxon(label="dummy1"),
                          label="dummy1",
                          edge_length=1.0)
        dummy_node.state = 0
        dummy_node.alive = False
        dummy_node.sampled = False
        dummy_node.is_sa = False
        dummy_node.is_sa_dummy_parent = True
        dummy_node.is_sa_lineage = False

        origin_node.add_child(dummy_node)
        
        # left child of dummy_node
        sa_node = Node(taxon=Taxon(label="sa1"),
                       label="sa1",
                       edge_length=0.0)
        sa_node.state = 0
        sa_node.alive = False
        sa_node.sampled = True
        sa_node.is_sa = True
        sa_node.is_sa_dummy_parent = False
        sa_node.is_sa_lineage = False

        # right child of dummy node, goes extinct
        brosc_node = Node(taxon=Taxon(label="brosc"),
                          label="brosc",
                          edge_length=1.0)
        brosc_node.state = 0
        brosc_node.alive = True
        brosc_node.sampled = True
        brosc_node.is_sa = False
        brosc_node.is_sa_dummy_parent = False
        brosc_node.is_sa_lineage = True

        dummy_node.add_child(sa_node)
        dummy_node.add_child(brosc_node)

        tr_complete = Tree(seed_node=origin_node)
        tr_complete.taxon_namespace.add_taxon(origin_node.taxon)
        tr_complete.taxon_namespace.add_taxon(dummy_node.taxon)
        tr_complete.taxon_namespace.add_taxon(sa_node.taxon)
        tr_complete.taxon_namespace.add_taxon(brosc_node.taxon)

        total_state_count = 1

        sa_global_time = 1.0
        time_to_sa_lineage_node = 1.0
        sa = pjsa.SampledAncestor(
            "sa1",
            "brosc",
            sa_global_time,
            time_to_lineage_node=time_to_sa_lineage_node)
        sa_lineage_dict = { "brosc": [sa] }
        
        max_age = 2.0

        ann_tr1 = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        ann_tr2 = copy.deepcopy(ann_tr1)

        # print(ann_tr1.tree.as_string(schema="newick",
        #                              suppress_annotations=True,
        #                              suppress_internal_node_labels=False,
        #                              suppress_internal_taxon_labels=False))

        tr_rec1 = ann_tr1.extract_reconstructed_tree(
            plotting_overhead=False,
            require_obs_both_sides=False)
        tr_rec2 = ann_tr2.extract_reconstructed_tree(
            plotting_overhead=False,
            require_obs_both_sides=True)

        tr_rec_str1 = tr_rec1.as_string(
            schema="newick",
            suppress_internal_taxon_labels=True,
            suppress_internal_node_labels=False,
            suppress_rooting=True)
        tr_rec_str2 = tr_rec2.as_string(
            schema="newick",
            suppress_internal_taxon_labels=True,
            suppress_internal_node_labels=False,
            suppress_rooting=True)

        # print(tr_rec_str1)
        # print(tr_rec_str2)

        self.assertEqual(
            ann_tr1.tree.as_string(
                schema="newick",
                suppress_annotations=True,
                suppress_internal_taxon_labels=True),
                "((sa1:0.0,brosc:1.0)dummy1:1.0)origin:0.0;\n")
        self.assertEqual(tr_rec_str1, "(sa1:0.0,brosc:1.0)dummy1:0.0;\n")
        self.assertEqual(tr_rec_str2, ";\n")

    def test_extract_reconstructed_tree_origin_sa_before_root_tree_two_survive(self):
        """
        Test extraction of rec. tree from complete tree (has origin).

        The complete tree has two observable taxa, an SA, and then TWO
        surviving extant taxa. We re-root the reconstructed tree at
        the dummy node above the root.

        To see it on icytree:
        ((sa1:0.0,(sp1:1.0,sp2:1.0)root:1.0)dummy1:1.0)origin:0.0;
        """
        
        origin_node = Node(taxon=Taxon(label="origin"),
                           label="origin",
                           edge_length=0.0)
        origin_node.state = 0
        origin_node.alive = False
        origin_node.sampled = False
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        # only child of origin node
        dummy_node = Node(taxon=Taxon(label="dummy1"),
                          label="dummy1",
                          edge_length=1.0)
        dummy_node.state = 0
        dummy_node.alive = False
        dummy_node.sampled = False
        dummy_node.is_sa = False
        dummy_node.is_sa_dummy_parent = True
        dummy_node.is_sa_lineage = False
        
        # left child of dummy_node
        sa_node = Node(taxon=Taxon(label="sa1"),
                       label="sa1",
                       edge_length=0.0)
        sa_node.state = 0
        sa_node.alive = False
        sa_node.sampled = True
        sa_node.is_sa = True
        sa_node.is_sa_dummy_parent = False
        sa_node.is_sa_lineage = False

        # right child of dummy node
        root_node = Node(taxon=Taxon(label="root"),
                         label="root",
                         edge_length=1.0)
        root_node.state = 0
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        # left child of nd1
        extant_sp1 = Node(taxon=Taxon(label="sp1"),
                          label="sp1",
                          edge_length=1.0)
        extant_sp1.state = 0
        extant_sp1.alive = True
        extant_sp1.sampled = True
        extant_sp1.is_sa = False
        extant_sp1.is_sa_dummy_parent = False
        extant_sp1.is_sa_lineage = False

        # right child of nd1
        extant_sp2 = Node(taxon=Taxon(label="sp2"),
                          label="sp2",
                          edge_length=1.0)
        extant_sp2.state = 0
        extant_sp2.alive = True
        extant_sp2.sampled = True
        extant_sp2.is_sa = False
        extant_sp2.is_sa_dummy_parent = False
        extant_sp2.is_sa_lineage = False

        root_node.add_child(extant_sp1)
        root_node.add_child(extant_sp2)

        dummy_node.add_child(sa_node)
        dummy_node.add_child(root_node)

        origin_node.add_child(dummy_node)

        tr_complete = Tree(seed_node=origin_node)
        tr_complete.taxon_namespace.add_taxon(origin_node.taxon)
        tr_complete.taxon_namespace.add_taxon(dummy_node.taxon)
        tr_complete.taxon_namespace.add_taxon(sa_node.taxon)
        tr_complete.taxon_namespace.add_taxon(root_node.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp2.taxon)

        total_state_count = 1

        sa_global_time = 1.0
        time_to_sa_lineage_node = 1.0
        sa = pjsa.SampledAncestor(
            "sa1",
            "root",
            sa_global_time,
            time_to_lineage_node=time_to_sa_lineage_node)
        sa_lineage_dict = { "root": [sa] }
        
        max_age = 2.0

        ann_tr1 = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        ann_tr2 = copy.deepcopy(ann_tr1)

        # looking at it
        # print(ann_tr1.tree.as_string(
        #     schema="newick",
        #     suppress_annotations=True,
        #     suppress_internal_node_labels=False,
        #     suppress_internal_taxon_labels=True))

        tr_rec1 = ann_tr1.extract_reconstructed_tree(
            plotting_overhead=False,
            require_obs_both_sides=False)
        tr_rec2 = ann_tr2.extract_reconstructed_tree(
            plotting_overhead=False,
            require_obs_both_sides=True)

        tr_rec_str1 = tr_rec1.as_string(
            schema="newick",
            suppress_internal_taxon_labels=True,
            suppress_internal_node_labels=False,
            suppress_rooting=True)
        # print(tr_rec_str1)
        
        tr_rec_str2 = tr_rec2.as_string(
            schema="newick",
            suppress_internal_taxon_labels=True,
            suppress_internal_node_labels=False,
            suppress_rooting=True)
        # print(tr_rec_str2)

        self.assertEqual(
            ann_tr1.tree.as_string(
                schema="newick",
                suppress_annotations=True,
                suppress_internal_taxon_labels=True),
                ("((sa1:0.0,(sp1:1.0,sp2:1.0)root:1.0)dummy1:1.0)"
                 "origin:0.0;\n"))
        self.assertEqual(tr_rec_str1,
                         "(sa1:0.0,(sp1:1.0,sp2:1.0)root:1.0)dummy1:0.0;\n")
        self.assertEqual(tr_rec_str2,
                         "(sa1:0.0,(sp1:1.0,sp2:1.0)root:1.0)dummy1:0.0;\n")

    def test_extract_reconstructed_tree_origin_one_root_side_dies_three_survive(self):
        """
        Test extraction of rec. tree from complete tree (has origin).

        The complete tree has three surviving taxa on only one side of the tree.
        We re-root the reconstructed tree at the MRCA of those three taxa.

        To see it on icytree:
        ((sp1:0.8,(sp2:1.0,(sp3:0.5,sp4:0.5)nd2:0.5)nd1:1.0)root:1.0)origin:0.0;
        """
        
        origin_node = Node(taxon=Taxon(label="origin"),
                           label="origin",
                           edge_length=0.0)
        origin_node.state = 0
        origin_node.alive = False
        origin_node.sampled = False
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        root_node = Node(taxon=Taxon(label="root"),
                         label="root",
                         edge_length=1.0)
        root_node.state = 0
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        origin_node.add_child(root_node)

        # left child of root_node
        extinct_sp1 = Node(taxon=Taxon(label="sp1"),
                           label="sp1",
                           edge_length=0.8)
        extinct_sp1.state = 0
        extinct_sp1.alive = False
        extinct_sp1.sampled = False
        extinct_sp1.is_sa = False
        extinct_sp1.is_sa_dummy_parent = False
        extinct_sp1.is_sa_lineage = False

        # right child of root_node
        internal_node1 = Node(taxon=Taxon(label="nd1"),
                              label="nd1",
                              edge_length=1.0)
        internal_node1.state = 0
        internal_node1.alive = True
        internal_node1.sampled = True
        internal_node1.is_sa = False
        internal_node1.is_sa_dummy_parent = False
        internal_node1.is_sa_lineage = False

        root_node.add_child(extinct_sp1)
        root_node.add_child(internal_node1)

        # left child of nd1
        extant_sp1 = Node(taxon=Taxon(label="sp2"),
                          label="sp2",
                          edge_length=1.0)
        extant_sp1.state = 0
        extant_sp1.alive = True
        extant_sp1.sampled = True
        extant_sp1.is_sa = False
        extant_sp1.is_sa_dummy_parent = False
        extant_sp1.is_sa_lineage = False

        # right child of nd1
        internal_node2 = Node(taxon=Taxon(label="nd2"),
                              label="nd2",
                              edge_length=0.5)
        internal_node2.state = 0
        internal_node2.alive = False
        internal_node2.sampled = False
        internal_node2.is_sa = False
        internal_node2.is_sa_dummy_parent = False
        internal_node2.is_sa_lineage = False

        internal_node1.add_child(extant_sp1)
        internal_node1.add_child(internal_node2)

        # left child of nd2
        extant_sp2 = Node(taxon=Taxon(label="sp3"),
                          label="sp3",
                          edge_length=0.5)
        extant_sp2.state = 0
        extant_sp2.alive = True
        extant_sp2.sampled = True
        extant_sp2.is_sa = False
        extant_sp2.is_sa_dummy_parent = False
        extant_sp2.is_sa_lineage = False

        # right child of nd2
        extant_sp3 = Node(taxon=Taxon(label="sp4"),
                          label="sp4",
                          edge_length=0.5)
        extant_sp3.state = 0
        extant_sp3.alive = True
        extant_sp3.sampled = True
        extant_sp3.is_sa = False
        extant_sp3.is_sa_dummy_parent = False
        extant_sp3.is_sa_lineage = False

        internal_node2.add_child(extant_sp2)
        internal_node2.add_child(extant_sp3)

        tr_complete = Tree(seed_node=origin_node)
        tr_complete.taxon_namespace.add_taxon(origin_node.taxon)
        tr_complete.taxon_namespace.add_taxon(root_node.taxon)
        tr_complete.taxon_namespace.add_taxon(internal_node1)
        tr_complete.taxon_namespace.add_taxon(internal_node2)
        tr_complete.taxon_namespace.add_taxon(extinct_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp2.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp3.taxon)

        total_state_count = 1
        
        max_age = 2.0

        ann_tr1 = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            max_age=max_age,
            epsilon=1e-12)

        ann_tr2 = copy.deepcopy(ann_tr1)

        # print(ann_tr1.tree.as_string(schema="newick",
        #                              suppress_annotations=True,
        #                              suppress_internal_node_labels=False,
        #                              suppress_internal_taxon_labels=True))

        tr_rec1 = ann_tr1.extract_reconstructed_tree(
            plotting_overhead=False,
            require_obs_both_sides=False)
        tr_rec2 = ann_tr2.extract_reconstructed_tree(
            plotting_overhead=False,
            require_obs_both_sides=True)

        tr_rec_str1 = tr_rec1.as_string(
            schema="newick",
            suppress_internal_taxon_labels=True,
            suppress_internal_node_labels=False,
            suppress_rooting=True)
        
        tr_rec_str2 = tr_rec2.as_string(
            schema="newick",
            suppress_internal_taxon_labels=True,
            suppress_internal_node_labels=False)

        self.assertEqual(
            ann_tr1.tree.as_string(
                schema="newick",
                suppress_annotations=True,
                suppress_internal_taxon_labels=True),
                ("((sp1:0.8,(sp2:1.0,(sp3:0.5,sp4:0.5)nd2:0.5)nd1:1.0)"
                 "root:1.0)origin:0.0;\n"))
        self.assertEqual(tr_rec_str1,
                         "(sp2:1.0,(sp3:0.5,sp4:0.5)nd2:0.5)nd1:0.0;\n")
        self.assertEqual(tr_rec_str2, ";\n")

    def test_extract_reconstructed_tree_origin_one_root_side_dies_three_survive_sa_before_root(self):
            """
            Test extraction of rec. tree from complete tree (has origin).

            The complete tree has three surviving taxa on only one side of the tree.
            There is a sampled ancestor before the root.
            We re-root the reconstructed tree at the MRCA of those three, which is
            the sampled ancestor! (but because we represent sampled ancestors as
            leaves, we actually re-root the reconstructed tree at the dummy node
            parent to the sampled ancestor).

            To see it on icytree:
            ((sa1:0.0,(sp1:0.8,(sp2:1.0,(sp3:0.5,sp4:0.5)nd2:0.5)nd1:1.0)root:1.0)dummy1:1.0)origin:0.0;
            """

            origin_node = Node(taxon=Taxon(label="origin"),
                           label="origin",
                           edge_length=0.0)
            origin_node.state = 0
            origin_node.alive = False
            origin_node.sampled = False
            origin_node.is_sa = False
            origin_node.is_sa_dummy_parent = False
            origin_node.is_sa_lineage = False

            # only child of origin node
            dummy_node = Node(taxon=Taxon(label="dummy1"),
                            label="dummy1",
                            edge_length=1.0)
            dummy_node.state = 0
            dummy_node.alive = False
            dummy_node.sampled = False
            dummy_node.is_sa = False
            dummy_node.is_sa_dummy_parent = True
            dummy_node.is_sa_lineage = False
            
            # left child of dummy_node
            sa_node = Node(taxon=Taxon(label="sa1"),
                        label="sa1",
                        edge_length=0.0)
            sa_node.state = 0
            sa_node.alive = False
            sa_node.sampled = True
            sa_node.is_sa = True
            sa_node.is_sa_dummy_parent = False
            sa_node.is_sa_lineage = False

            # right child of dummy_node
            root_node = Node(taxon=Taxon(label="root"),
                            label="root",
                            edge_length=1.0)
            root_node.state = 0
            root_node.alive = False
            root_node.sampled = False
            root_node.is_sa = False
            root_node.is_sa_dummy_parent = False
            root_node.is_sa_lineage = False

            dummy_node.add_child(sa_node)
            dummy_node.add_child(root_node)

            origin_node.add_child(dummy_node)

            # left child of root_node
            extinct_sp1 = Node(taxon=Taxon(label="sp1"),
                            label="sp1",
                            edge_length=0.8)
            extinct_sp1.state = 0
            extinct_sp1.alive = False
            extinct_sp1.sampled = False
            extinct_sp1.is_sa = False
            extinct_sp1.is_sa_dummy_parent = False
            extinct_sp1.is_sa_lineage = False

            # right child of root_node
            internal_node1 = Node(taxon=Taxon(label="nd1"),
                                label="nd1",
                                edge_length=1.0)
            internal_node1.state = 0
            internal_node1.alive = False
            internal_node1.sampled = False
            internal_node1.is_sa = False
            internal_node1.is_sa_dummy_parent = False
            internal_node1.is_sa_lineage = False

            root_node.add_child(extinct_sp1)
            root_node.add_child(internal_node1)

            # left child of nd1
            extant_sp1 = Node(taxon=Taxon(label="sp2"),
                            label="sp2",
                            edge_length=1.0)
            extant_sp1.state = 0
            extant_sp1.alive = True
            extant_sp1.sampled = True
            extant_sp1.is_sa = False
            extant_sp1.is_sa_dummy_parent = False
            extant_sp1.is_sa_lineage = False

            # right child of nd1
            internal_node2 = Node(taxon=Taxon(label="nd2"),
                                label="nd2",
                                edge_length=0.5)
            internal_node2.state = 0
            internal_node2.alive = False
            internal_node2.sampled = False
            internal_node2.is_sa = False
            internal_node2.is_sa_dummy_parent = False
            internal_node2.is_sa_lineage = False

            internal_node1.add_child(extant_sp1)
            internal_node1.add_child(internal_node2)

            # left child of nd2
            extant_sp2 = Node(taxon=Taxon(label="sp3"),
                            label="sp3",
                            edge_length=0.5)
            extant_sp2.state = 0
            extant_sp2.alive = True
            extant_sp2.sampled = True
            extant_sp2.is_sa = False
            extant_sp2.is_sa_dummy_parent = False
            extant_sp2.is_sa_lineage = False

            # right child of nd2
            extant_sp3 = Node(taxon=Taxon(label="sp4"),
                            label="sp4",
                            edge_length=0.5)
            extant_sp3.state = 0
            extant_sp3.alive = True
            extant_sp3.sampled = True
            extant_sp3.is_sa = False
            extant_sp3.is_sa_dummy_parent = False
            extant_sp3.is_sa_lineage = False

            internal_node2.add_child(extant_sp2)
            internal_node2.add_child(extant_sp3)

            tr_complete = Tree(seed_node=origin_node)
            tr_complete.taxon_namespace.add_taxon(origin_node.taxon)
            tr_complete.taxon_namespace.add_taxon(dummy_node.taxon)
            tr_complete.taxon_namespace.add_taxon(sa_node.taxon)
            tr_complete.taxon_namespace.add_taxon(root_node.taxon)
            tr_complete.taxon_namespace.add_taxon(internal_node1)
            tr_complete.taxon_namespace.add_taxon(internal_node2)
            tr_complete.taxon_namespace.add_taxon(extinct_sp1.taxon)
            tr_complete.taxon_namespace.add_taxon(extant_sp1.taxon)
            tr_complete.taxon_namespace.add_taxon(extant_sp2.taxon)
            tr_complete.taxon_namespace.add_taxon(extant_sp3.taxon)

            sa_global_time = 1.0
            time_to_sa_lineage_node = 1.0
            sa = pjsa.SampledAncestor(
                "sa1",
                "root",
                sa_global_time,
                time_to_lineage_node=time_to_sa_lineage_node)
            sa_lineage_dict = { "root": [sa] }

            total_state_count = 1
            
            max_age = 3.0

            ann_tr1 = pjtr.AnnotatedTree(
                tr_complete,
                DiscreteTrait(total_state_count),
                start_at_origin=True,
                max_age=max_age,
                sa_lineage_dict=sa_lineage_dict,
                epsilon=1e-12)

            ann_tr2 = copy.deepcopy(ann_tr1)

            # looking at it
            # print(ann_tr1.tree.as_string(
            #     schema="newick",
            #     suppress_annotations=True,
            #     suppress_internal_node_labels=False,
            #     suppress_internal_taxon_labels=True))

            tr_rec1 = ann_tr1.extract_reconstructed_tree(
                plotting_overhead=False,
                require_obs_both_sides=False)
            tr_rec2 = ann_tr2.extract_reconstructed_tree(
                plotting_overhead=False,
                require_obs_both_sides=True)

            tr_rec_str1 = tr_rec1.as_string(
                schema="newick",
                suppress_internal_taxon_labels=True,
                suppress_internal_node_labels=False,
                suppress_rooting=True)
            # print(tr_rec_str1)
            
            tr_rec_str2 = tr_rec2.as_string(
                schema="newick",
                suppress_internal_taxon_labels=True,
                suppress_internal_node_labels=False,
                suppress_rooting=True)
            # print(tr_rec_str2)

            self.assertEqual(
            ann_tr1.tree.as_string(
                schema="newick",
                suppress_annotations=True,
                suppress_internal_taxon_labels=True),
                ("((sa1:0.0,(sp1:0.8,(sp2:1.0,(sp3:0.5,sp4:0.5)nd2:0.5)nd1:1.0"
                 ")root:1.0)dummy1:1.0)origin:0.0;\n"))
            self.assertEqual(
                tr_rec_str1,
                ("(sa1:0.0,(sp2:1.0,(sp3:0.5,sp4:0.5)nd2:0.5)nd1:2.0)dummy1:0.0"
                 ";\n"))
            self.assertEqual(tr_rec_str2, ";\n")

    def test_make_at_dict_reflect_rec_tree_origin_two_extinct(self) -> None:
        """Test method that updates 'rec_tr_at_dict' for reconstructed tree.

        See also test_tree_plotting.test_harder_geosse_rec_tr()
        """

        origin_node = Node(taxon=Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node.state = 2  # AB
        origin_node.annotations.add_bound_attribute("state")
        origin_node.index = 0
        origin_node.annotations.add_bound_attribute("index")
        origin_node.alive = False
        origin_node.sampled = False
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        root_node = Node(taxon=Taxon(label="root"), label="root", edge_length=1.0)
        root_node.state = 2  # AB
        root_node.annotations.add_bound_attribute("state")
        root_node.index = 1
        root_node.annotations.add_bound_attribute("index")
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        # left child of root_node
        internal_node1 = Node(taxon=Taxon(label="nd6"),
                              label="nd6",
                              edge_length=1.0)
        internal_node1.state = 2  # AB
        internal_node1.annotations.add_bound_attribute("state")
        internal_node1.index = 2
        internal_node1.annotations.add_bound_attribute("index")
        internal_node1.alive = False
        internal_node1.sampled = False
        internal_node1.is_sa = False
        internal_node1.is_sa_dummy_parent = False
        internal_node1.is_sa_lineage = False

        # right child of root_node
        internal_node3 = Node(taxon=Taxon(label="nd8"),
                              label="nd8",
                              edge_length=1.0)
        internal_node3.state = 2  # AB
        internal_node3.annotations.add_bound_attribute("state")
        internal_node3.index = 3
        internal_node3.annotations.add_bound_attribute("index")
        internal_node3.alive = False
        internal_node3.sampled = False
        internal_node3.is_sa = False
        internal_node3.is_sa_dummy_parent = False
        internal_node3.is_sa_lineage = False

        # left child of internal_node1
        internal_node2 = Node(taxon=Taxon(label="nd5"),
                              label="nd5",
                              edge_length=1.0)
        internal_node2.state = 2  # AB
        internal_node2.annotations.add_bound_attribute("state")
        internal_node2.index = 4
        internal_node2.annotations.add_bound_attribute("index")
        internal_node2.alive = False
        internal_node2.sampled = False
        internal_node2.is_sa = False
        internal_node2.is_sa_dummy_parent = False
        internal_node2.is_sa_lineage = False

        # right child of internal_node1
        extinct_sp3 = Node(taxon=Taxon(label="sp3"),
                          label="sp3",
                          edge_length=2.0)
        extinct_sp3.state = 0  # A
        extinct_sp3.annotations.add_bound_attribute("state")
        extinct_sp3.index = 6
        extinct_sp3.annotations.add_bound_attribute("index")
        extinct_sp3.alive = False
        extinct_sp3.sampled = False
        extinct_sp3.is_sa = False
        extinct_sp3.is_sa_dummy_parent = False
        extinct_sp3.is_sa_lineage = False

        # left child of internal_node2
        extinct_sp1 = Node(taxon=Taxon(label="sp1"),
                          label="sp1",
                          edge_length=1.0)
        extinct_sp1.state = 1  # B
        extinct_sp1.annotations.add_bound_attribute("state")
        extinct_sp1.index = 7
        extinct_sp1.annotations.add_bound_attribute("index")
        extinct_sp1.alive = False
        extinct_sp1.sampled = False
        extinct_sp1.is_sa = False
        extinct_sp1.is_sa_dummy_parent = False
        extinct_sp1.is_sa_lineage = False

        # right child of internal_node2
        extant_sp2 = Node(taxon=Taxon(label="sp2"),
                          label="sp2",
                          edge_length=2.0)
        extant_sp2.state = 2  # AB
        extant_sp2.annotations.add_bound_attribute("state")
        extant_sp2.index = 8
        extant_sp2.annotations.add_bound_attribute("index")
        extant_sp2.alive = True
        extant_sp2.sampled = True
        extant_sp2.is_sa = False
        extant_sp2.is_sa_dummy_parent = False
        extant_sp2.is_sa_lineage = False

        # left child of internal_node1
        internal_node4 = Node(taxon=Taxon(label="nd7"),
                              label="nd7",
                              edge_length=1.0)
        internal_node4.state = 2  # AB
        internal_node4.annotations.add_bound_attribute("state")
        internal_node4.index = 5
        internal_node4.annotations.add_bound_attribute("index")
        internal_node4.alive = False
        internal_node4.sampled = False
        internal_node4.is_sa = False
        internal_node4.is_sa_dummy_parent = False
        internal_node4.is_sa_lineage = False

        # right child of internal_node1
        extinct_sp6 = Node(taxon=Taxon(label="sp6"),
                           label="sp6",
                           edge_length=2.0)
        extinct_sp6.state = 0  # A
        extinct_sp6.annotations.add_bound_attribute("state")
        extinct_sp6.index = 9
        extinct_sp6.annotations.add_bound_attribute("index")
        extinct_sp6.alive = False
        extinct_sp6.sampled = False
        extinct_sp6.is_sa = False
        extinct_sp6.is_sa_dummy_parent = False
        extinct_sp6.is_sa_lineage = False

        # left child of internal_node4
        extinct_sp4 = Node(taxon=Taxon(label="sp4"),
                           label="sp4",
                           edge_length=1.0)
        extinct_sp4.state = 1  # B
        extinct_sp4.annotations.add_bound_attribute("state")
        extinct_sp4.index = 10
        extinct_sp4.annotations.add_bound_attribute("index")
        extinct_sp4.alive = False
        extinct_sp4.sampled = False
        extinct_sp4.is_sa = False
        extinct_sp4.is_sa_dummy_parent = False
        extinct_sp4.is_sa_lineage = False

        # right child of internal_node4
        extant_sp5 = Node(taxon=Taxon(label="sp5"),
                          label="sp5",
                          edge_length=2.0)
        extant_sp5.state = 2  # AB
        extant_sp5.annotations.add_bound_attribute("state")
        extant_sp5.index = 11
        extant_sp5.annotations.add_bound_attribute("index")
        extant_sp5.alive = True
        extant_sp5.sampled = True
        extant_sp5.is_sa = False
        extant_sp5.is_sa_dummy_parent = False
        extant_sp5.is_sa_lineage = False

        # building topology
        internal_node2.add_child(extinct_sp1)
        internal_node2.add_child(extant_sp2)

        internal_node1.add_child(internal_node2)  # 'nd5'
        internal_node1.add_child(extinct_sp3)

        internal_node4.add_child(extinct_sp4)
        internal_node4.add_child(extant_sp5)

        internal_node3.add_child(internal_node4)  # 'nd7'
        internal_node3.add_child(extinct_sp6)

        root_node.add_child(internal_node1)  # 'nd6'
        root_node.add_child(internal_node3)  # 'nd8'

        origin_node.add_child(root_node)

        # wrapping up tree
        tr_complete = Tree(seed_node=origin_node)
        tr_complete.taxon_namespace.add_taxon(origin_node.taxon)
        tr_complete.taxon_namespace.add_taxon(root_node.taxon)
        tr_complete.taxon_namespace.add_taxon(internal_node1.taxon)
        tr_complete.taxon_namespace.add_taxon(internal_node2.taxon)
        tr_complete.taxon_namespace.add_taxon(internal_node3.taxon)
        tr_complete.taxon_namespace.add_taxon(internal_node4.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp3.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp4.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp6.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp2.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp5.taxon)

        at1_1_1 = pjat.AttributeTransition("state", "nd5", 2.0, 2, 1)
        at1_1_2 = pjat.AttributeTransition("state", "nd5", 2.5, 1, 2)
        at1_1_3 = pjat.AttributeTransition("state", "nd7", 2.0, 2, 1)
        at1_1_4 = pjat.AttributeTransition("state", "nd7", 2.5, 1, 2)
        at1_2_1 = pjat.AttributeTransition("state", "sp3", 2.0, 2, 0)
        at1_2_2 = pjat.AttributeTransition("state", "sp6", 2.0, 2, 0)
        at2 = pjat.AttributeTransition("state", "sp2", 3.0, 2, 0)
        at3 = pjat.AttributeTransition("state", "sp2", 4.0, 0, 2)
        at4 = pjat.AttributeTransition("state", "sp5", 3.0, 2, 0)
        at5 = pjat.AttributeTransition("state", "sp5", 4.0, 0, 2)
        at_dict = {
            "nd5": [at1_1_1, at1_1_2],
            "nd7": [at1_1_3, at1_1_4],
            "sp3": [at1_2_1],
            "sp6": [at1_2_2],
            "sp2": [at2, at3],
            "sp5": [at4, at5]
        }

        # internal_node2
        clado_at1 = pjat.AttributeTransition("state",
                                             subtending_node_label="nd5",
                                             global_time=3.0,
                                             from_state=2,
                                             to_state=0,
                                             to_state2=1)
        clado_at2 = pjat.AttributeTransition("state",
                                             subtending_node_label="nd7",
                                             global_time=3.0,
                                             from_state=2,
                                             to_state=0,
                                             to_state2=1)
        clado_at3 = pjat.AttributeTransition("state",
                                             subtending_node_label="nd6",
                                             global_time=2.0,
                                             from_state=2,
                                             to_state=0,
                                             to_state2=1)
        clado_at4 = pjat.AttributeTransition("state",
                                             subtending_node_label="nd8",
                                             global_time=2.0,
                                             from_state=2,
                                             to_state=0,
                                             to_state2=1)
        clado_at_dict = {
            "nd5": clado_at1,
            "nd7": clado_at2,
            "nd6": clado_at3,
            "nd8": clado_at4
        }


        total_state_count = 3
        max_age = 5.0
        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(total_state_count),
            at_dict=at_dict,
            clado_at_dict=clado_at_dict,
            start_at_origin=True,
            max_age=max_age,
            epsilon=1e-12)

        ann_tr.populate_nd_attr_dict("state")

        # ann_tr_str = \
        #     ann_tr.tree.as_string(
        #         schema="newick",
        #         suppress_internal_taxon_labels=True,
        #         suppress_internal_node_labels=False)
        # print(ann_tr_str)

        # updates rec tree-related members in 'ann_tr'
        tr_rec = \
            ann_tr.extract_reconstructed_tree(
                plotting_overhead=True,
                require_obs_both_sides=False)

        # tr_rec_str = \
        #     tr_rec.as_string(
        #     schema="newick",
        #     suppress_internal_taxon_labels=True,
        #     suppress_internal_node_labels=False,
        #     suppress_rooting=True)
        # print(tr_rec_str)

        rec_tr_at_dict = ann_tr.rec_tr_at_dict

        # check that the right nodes are in the rec tree's at_dict
        nd_names_in_rec_tr_at_dict = set(rec_tr_at_dict.keys())
        exp_nd_names_in_rec_tr_at_dict = {"sp2", "sp5"}

        self.assertEqual(nd_names_in_rec_tr_at_dict,
                         exp_nd_names_in_rec_tr_at_dict)

        # check that the values are correct
        at_times_dict_for_testing = dict()
        for nd_name, at_list in rec_tr_at_dict.items():
            global_times = list()

            for at in at_list:
                global_times.append(at.global_time)

            at_times_dict_for_testing[nd_name] = \
                global_times

        self.assertEqual(at_times_dict_for_testing,
                         {"sp2": [1.0, 1.5, 2.0, 3.0],
                         "sp5": [1.0, 1.5, 2.0, 3.0]})

        self.assertEqual({}, ann_tr.rec_tr_clado_at_dict)


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
    # $ python3.11 tests/data/test_tree_extract_reconstructed.py
    # 
    # or
    #
    # $ python3.11 -m tests.data.test_tree_extract_reconstructed
    #
    # or 
    #
    # $ python3.11 -m unittest tests.data.test_tree_extract_reconstructed.TestExtractReconstructedTree.test_extract_reconstructed_tree_origin_one_root_side_dies_three_survive_sa_before_root

    unittest.main()