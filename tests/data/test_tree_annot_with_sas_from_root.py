from phylojunction.data.trait import DiscreteTrait
import unittest
from dendropy import Tree, Node, Taxon

# pj imports
import phylojunction.data.tree as pjtr
import phylojunction.data.sampled_ancestor as pjsa

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestAnnotateTreeWithSAsFromRoot(unittest.TestCase):

    def test_node_counting_oneSA_survives(self):
        """
        Test counting of extant, extinct and SA nodes

        Tree: root + sp1 + sp2, then sp1 undergoes ancestor sampling, which causes a "dummy node"
        to replace the left child, and the the children of the "dummy node" will in turn become sp1 and
        a SA node. Then sp2 goes extinct.

        To see it on icytree: ((sa1:0.0,sp1:1.0)dummy1:1.0,sp2:1.5)root:0.0;
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

        # right child of root node
        extant_sp2 = Node(taxon=Taxon(label="sp2"), label="sp2", edge_length=1.5)
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
        sa_node.sampled = False
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
        extant_sp1.is_sa_lineage = False
        
        dummy_node.add_child(sa_node)
        dummy_node.add_child(extant_sp1)

        tr_sa_survives = Tree(seed_node=root_node)
        tr_sa_survives.taxon_namespace.add_taxon(root_node.taxon)
        tr_sa_survives.taxon_namespace.add_taxon(dummy_node.taxon)
        tr_sa_survives.taxon_namespace.add_taxon(sa_node.taxon)
        tr_sa_survives.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_sa_survives.taxon_namespace.add_taxon(extant_sp2.taxon)

        # debugging
        # print("tr_sa_survives.seed_age = " + str(tr_sa_survives.max_distance_from_root()))
        # print(tr_sa_survives.as_string(schema="newick"))

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

        ann_tr_sa_survives_max_age = pjtr.AnnotatedTree(
            tr_sa_survives,
            DiscreteTrait(total_state_count),
            start_at_origin=False,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        self.assertEqual(
            ann_tr_sa_survives_max_age.n_extant_terminal_nodes,
            1,
            "Count of terminal extant nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_survives_max_age.n_extinct_terminal_nodes,
            1,
            "Count of terminal extinct nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_survives_max_age.n_sa_nodes,
            1,
            "Count of sampled ancestor nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_survives_max_age.origin_edge_length,
            0.0,
            "Length of origin edge should be 0.0.")

        tr_sa_survives_no_max_age = Tree(seed_node=root_node)
        tr_sa_survives_no_max_age.taxon_namespace.add_taxon(root_node.taxon)
        tr_sa_survives_no_max_age.taxon_namespace.add_taxon(dummy_node.taxon)
        tr_sa_survives_no_max_age.taxon_namespace.add_taxon(sa_node.taxon)
        tr_sa_survives_no_max_age.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_sa_survives_no_max_age.taxon_namespace.add_taxon(extant_sp2.taxon)

        ann_tr_sa_survives_no_max_age = pjtr.AnnotatedTree(
            tr_sa_survives_no_max_age,
            DiscreteTrait(total_state_count),
            start_at_origin=False,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        self.assertEqual(
            ann_tr_sa_survives_no_max_age.n_extant_terminal_nodes,
            1,
            "Count of terminal extant nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_survives_no_max_age.n_extinct_terminal_nodes,
            1,
            "Count of terminal extinct nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_survives_no_max_age.n_sa_nodes,
            1,
            "Count of sampled ancestor nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_survives_no_max_age.origin_edge_length,
            0.0,
            "Length of origin edge should be 0.0.")    
    

    def test_node_counting_oneSA_dies(self):
        """
        Test counting of extant, extinct and SA nodes

        Tree: root + sp1 + sp2, then sp1 undergoes ancestor sampling, which causes a "dummy node"
        to replace the left child, and the the children of the "dummy node" will in turn become sp1 and
        a SA node. Then both sp1 and sp2 go extinct.

        To see it on icytree: ((sa1:0.0,sp1:0.5)dummy1:1.0,sp2:1.5)root:0.0;
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

        # right child of root node
        extant_sp2 = Node(taxon=Taxon(label="sp2"), label="sp2", edge_length=1.5)
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
        sa_node.sampled = False
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

        tr_sa_dies = Tree(seed_node=root_node)
        tr_sa_dies.taxon_namespace.add_taxon(root_node.taxon)
        tr_sa_dies.taxon_namespace.add_taxon(dummy_node.taxon)
        tr_sa_dies.taxon_namespace.add_taxon(sa_node.taxon)
        tr_sa_dies.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_sa_dies.taxon_namespace.add_taxon(extant_sp2.taxon)

        # debugging
        # print("tr_sa_dies.seed_age = " + str(tr_sa_dies.max_distance_from_root()))
        # print(tr_sa_dies.as_string(schema="newick"))

        total_state_count = 1

        sa_global_time = 1.0
        time_to_sa_lineage_node = 1.0
        sa = pjsa.SampledAncestor("sa1", "sp1", sa_global_time, time_to_lineage_node=time_to_sa_lineage_node)
        sa_lineage_dict = { "sp1": [sa] }
        
        max_age = 2.0

        ann_tr_sa_dies_max_age = pjtr.AnnotatedTree(
            tr_sa_dies,
            DiscreteTrait(total_state_count),
            start_at_origin=False,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        self.assertEqual(ann_tr_sa_dies_max_age.n_extant_terminal_nodes, 0, "Count of terminal extant nodes should be 0.")
        self.assertEqual(ann_tr_sa_dies_max_age.n_extinct_terminal_nodes, 2, "Count of terminal extinct nodes should be 2.")
        self.assertEqual(ann_tr_sa_dies_max_age.n_sa_nodes, 1, "Count of sampled ancestor nodes should be 1.")
        self.assertEqual(ann_tr_sa_dies_max_age.origin_edge_length, 0.0, "Length of origin edge should be 0.0.")

        tr_sa_dies_no_max_age = Tree(seed_node=root_node)
        tr_sa_dies.taxon_namespace.add_taxon(root_node.taxon)
        tr_sa_dies.taxon_namespace.add_taxon(dummy_node.taxon)
        tr_sa_dies.taxon_namespace.add_taxon(sa_node.taxon)
        tr_sa_dies.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_sa_dies.taxon_namespace.add_taxon(extant_sp2.taxon)

        ann_tr_sa_dies_no_max_age = pjtr.AnnotatedTree(
            tr_sa_dies_no_max_age,
            DiscreteTrait(total_state_count),
            start_at_origin=False,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        self.assertEqual(ann_tr_sa_dies_no_max_age.n_extant_terminal_nodes, 0, "Count of terminal extant nodes should be 0.")
        self.assertEqual(ann_tr_sa_dies_no_max_age.n_extinct_terminal_nodes, 2, "Count of terminal extinct nodes should be 2.")
        self.assertEqual(ann_tr_sa_dies_no_max_age.n_sa_nodes, 1, "Count of sampled ancestor nodes should be 1.")
        self.assertEqual(ann_tr_sa_dies_no_max_age.origin_edge_length, 0.0, "Length of origin edge should be 0.0.")

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
    # $ python3 tests/data/test_tree_annot_with_sas_from_root.py
    # 
    # or
    #
    # $ python3 -m tests.data.test_tree_annot_with_sas_from_root
    #
    # or 
    #
    # $ python3 -m unittest tests.data.test_tree_annot_with_sas_from_root.TestAnnotateTreeWithSAsFromRoot.test_node_counting_oneSA_survives

    unittest.main()