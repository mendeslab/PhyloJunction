from phylojunction.data.trait import DiscreteTrait
import unittest
from dendropy import Tree, Node, Taxon

# pj imports
import phylojunction.data.tree as pjtr
import phylojunction.data.sampled_ancestor as pjsa

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestAnnotateTreeWithSAsFromOrigin(unittest.TestCase):

    def test_node_counting_oneSA_no_spn_survives_max_age(self):
        """
        Test counting of extant, extinct and SA nodes.

        Tree: origin + SA then survive until with and without max_age
        Note that the way "execute_sample_ancestor" is implemented, when this event happens before a
        birth event (giving rise to the root), the origin instead undergoes a "fake birth" event, where
        a dummy node is made its child, and two children are added to the dummy node:
        
        (1) SA,
        (2) A node between the origin and the root (the "brosc" node)
    
        Then (2) stays alive and potentially undergoes a proper birth event (when the root would be created)
        or another ancestor sampling event.
        Here, the "brosc" simply stays alive until the end of the process at 'max_age' or however the process
        terminates

        To see it on icytree: ((sa1:0.0,brosc:1.0)dummy1:1.0)origin:0.0;
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

        dummy_node = Node(taxon=Taxon(label="dummy1"),
                          label="dummy1",
                          edge_length=1.0)
        dummy_node.state = 0
        dummy_node.alive = False
        dummy_node.sampled = False
        dummy_node.is_sa = False
        dummy_node.is_sa_dummy_parent = True
        
        origin_node.add_child(dummy_node)
        
        # right child of dummy_node
        sa_node = Node(taxon=Taxon(label="sa1"),
                       label="sa1",
                       edge_length=0.0)
        sa_node.state = 0
        sa_node.alive = False
        sa_node.sampled = False
        sa_node.is_sa = True
        sa_node.is_sa_dummy_parent = False
        sa_node.is_sa_lineage = False

        # left child of dummy node
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

        tr_sa_no_spn_survives = Tree(seed_node=origin_node)
        tr_sa_no_spn_survives.taxon_namespace.add_taxon(origin_node.taxon)
        tr_sa_no_spn_survives.taxon_namespace.add_taxon(dummy_node.taxon)
        tr_sa_no_spn_survives.taxon_namespace.add_taxon(sa_node.taxon)
        tr_sa_no_spn_survives.taxon_namespace.add_taxon(brosc_node.taxon)       

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

        ann_tr_sa_no_spn_survives_max_age = pjtr.AnnotatedTree(
            tr_sa_no_spn_survives,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        self.assertEqual(
            ann_tr_sa_no_spn_survives_max_age.n_extant_terminal_nodes,
            1,
            "Count of terminal extant nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_no_spn_survives_max_age.n_extinct_terminal_nodes,
            0,
            "Count of terminal extinct nodes should be 0.")
        self.assertEqual(
            ann_tr_sa_no_spn_survives_max_age.n_sa_nodes,
            1,
            "Count of sampled ancestor nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_no_spn_survives_max_age.origin_edge_length,
            2.0,
            "Length of origin edge should be 2.0.")

        tr_sa_no_spn_survives_no_max_age = Tree(seed_node=origin_node)
        tr_sa_no_spn_survives_no_max_age.taxon_namespace.\
            add_taxon(origin_node.taxon)
        tr_sa_no_spn_survives_no_max_age.taxon_namespace.\
            add_taxon(dummy_node.taxon)
        tr_sa_no_spn_survives_no_max_age.taxon_namespace.\
            add_taxon(sa_node.taxon)
        tr_sa_no_spn_survives_no_max_age.taxon_namespace.\
            add_taxon(brosc_node.taxon)

        ann_tr_sa_no_spn_survives_no_max_age = pjtr.AnnotatedTree(
            tr_sa_no_spn_survives_no_max_age,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        self.assertEqual(
            ann_tr_sa_no_spn_survives_no_max_age.n_extant_terminal_nodes,
            1,
            "Count of terminal extant nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_no_spn_survives_no_max_age.n_extinct_terminal_nodes,
            0,
            "Count of terminal extinct nodes should be 0.")
        self.assertEqual(
            ann_tr_sa_no_spn_survives_no_max_age.n_sa_nodes,
            1,
            "Count of sampled ancestor nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_no_spn_survives_no_max_age.origin_edge_length,
            2.0,
            "Length of origin edge should be 2.0.")


    def test_node_counting_twoSAs_no_spn_survives(self):
        """
        Test counting of extant, extinct and SA nodes

        Tree: origin + SA + SA then survive until with max age.
        Note that the way "execute_sample_ancestor" is implemented, when this event happens before a
        birth event (giving rise to the root), the origin instead undergoes a "fake birth" event, where
        a dummy node is made its child, and two children are added to the dummy node:
        
        (1) SA,
        (2) A node between the origin and the root (the "brosc" node)
    
        Then (2) stays alive and potentially undergoes a proper birth event (when the root would be created)
        or another ancestor sampling event.
        Here, the "brosc" node undergoes another ancestor sampling event, and the resulting daughter "brosc"
        node stays alive:
        
        (i) until the end of the process at 'max_age' and no other events,
        (ii) there is no max_age

        To see it on icytree: ((sa1:0.0,(sa2:0.0,brosc:1.0)dummy2:0.5)dummy1:0.5)origin:0.0;
        """
        
        origin_node = Node(taxon=Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node.state = 0
        origin_node.alive = False
        origin_node.sampled = True
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        dummy_node1 = Node(taxon=Taxon(label="dummy1"), label="dummy1", edge_length=0.5)
        dummy_node1.state = 0
        dummy_node1.alive = False
        dummy_node1.sampled = False
        dummy_node1.is_sa = False
        dummy_node1.is_sa_dummy_parent = True
        
        origin_node.add_child(dummy_node1)
        
        # right child of dummy_node1 
        sa_node1 = Node(taxon=Taxon(label="sa1"), label="sa1", edge_length=0.0)
        sa_node1.state = 0
        sa_node1.alive = False
        sa_node1.sampled = False
        sa_node1.is_sa = True
        sa_node1.is_sa_dummy_parent = False
        sa_node1.is_sa_lineage = False

        # left child of dummy_node1
        dummy_node2 = Node(taxon=Taxon(label="dummy2"), label="dummy2", edge_length=0.5)
        dummy_node2.state = 0
        dummy_node2.alive = False
        dummy_node2.sampled = False
        dummy_node2.is_sa = False
        dummy_node2.is_sa_dummy_parent = True

        dummy_node1.add_child(sa_node1)
        dummy_node1.add_child(dummy_node2)
        
        # right child of dummy_node2
        sa_node2 = Node(taxon=Taxon(label="sa2"), label="sa2", edge_length=0.0)
        sa_node2.state = 0
        sa_node2.alive = False
        sa_node2.sampled = False
        sa_node2.is_sa = True
        sa_node2.is_sa_dummy_parent = False
        sa_node2.is_sa_lineage = False

        # left child of dummy_node2
        brosc_node = Node(taxon=Taxon(label="brosc"), label="brosc", edge_length=1.0)
        brosc_node.state = 0
        brosc_node.alive = True
        brosc_node.sampled = True
        brosc_node.is_sa = False
        brosc_node.is_sa_dummy_parent = False
        brosc_node.is_sa_lineage = True
        
        dummy_node2.add_child(sa_node2)
        dummy_node2.add_child(brosc_node)

        tr_2sas_no_spn_survives = Tree(seed_node=origin_node)
        tr_2sas_no_spn_survives.taxon_namespace.add_taxon(origin_node.taxon)
        tr_2sas_no_spn_survives.taxon_namespace.add_taxon(dummy_node1.taxon)
        tr_2sas_no_spn_survives.taxon_namespace.add_taxon(dummy_node2.taxon)
        tr_2sas_no_spn_survives.taxon_namespace.add_taxon(sa_node1.taxon)
        tr_2sas_no_spn_survives.taxon_namespace.add_taxon(sa_node2.taxon)
        tr_2sas_no_spn_survives.taxon_namespace.add_taxon(brosc_node.taxon)

        # debugging
        # print("tr_2sas_no_spn_survives.seed_age = " + str(tr_2sas_no_spn_survives.max_distance_from_root()))
        # print(tr_2sas_no_spn_survives.as_string(schema="newick"))

        total_state_count = 1

        sa1_global_time = 0.5
        time_to_sa1_lineage_node = 0.5
        sa1 = pjsa.SampledAncestor("sa1", "dummy_node2", sa1_global_time, time_to_lineage_node=time_to_sa1_lineage_node)

        sa2_global_time = 1.0
        time_to_sa2_lineage_node = 0.5
        sa2 = pjsa.SampledAncestor("sa2", "brosc", sa2_global_time, time_to_lineage_node=time_to_sa2_lineage_node)

        sa_lineage_dict = { "brosc": [sa1, sa2] }
        
        max_age = 2.0

        ann_tr_2sas_no_spn_survives_max_age = pjtr.AnnotatedTree(
            tr_2sas_no_spn_survives,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        self.assertEqual(ann_tr_2sas_no_spn_survives_max_age.n_extant_terminal_nodes, 1, "Count of terminal extant nodes should be 1.")
        self.assertEqual(ann_tr_2sas_no_spn_survives_max_age.n_extinct_terminal_nodes, 0, "Count of terminal extinct nodes should be 0.")
        self.assertEqual(ann_tr_2sas_no_spn_survives_max_age.n_sa_nodes, 2, "Count of sampled ancestor nodes should be 2.")
        self.assertEqual(ann_tr_2sas_no_spn_survives_max_age.origin_edge_length, 2.0, "Length of origin edge should be 2.0.")

        tr_2sas_no_spn_survives_no_max_age = Tree(seed_node=origin_node)
        tr_2sas_no_spn_survives_no_max_age.taxon_namespace.add_taxon(origin_node.taxon)
        tr_2sas_no_spn_survives_no_max_age.taxon_namespace.add_taxon(dummy_node1.taxon)
        tr_2sas_no_spn_survives_no_max_age.taxon_namespace.add_taxon(dummy_node2.taxon)
        tr_2sas_no_spn_survives_no_max_age.taxon_namespace.add_taxon(sa_node1.taxon)
        tr_2sas_no_spn_survives_no_max_age.taxon_namespace.add_taxon(sa_node2.taxon)
        tr_2sas_no_spn_survives_no_max_age.taxon_namespace.add_taxon(brosc_node.taxon)

        ann_tr_2sas_no_spn_survives_no_max_age = pjtr.AnnotatedTree(
            tr_2sas_no_spn_survives_no_max_age,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        self.assertEqual(ann_tr_2sas_no_spn_survives_no_max_age.n_extant_terminal_nodes, 1, "Count of terminal extant nodes should be 1.")
        self.assertEqual(ann_tr_2sas_no_spn_survives_no_max_age.n_extinct_terminal_nodes, 0, "Count of terminal extinct nodes should be 0.")
        self.assertEqual(ann_tr_2sas_no_spn_survives_no_max_age.n_sa_nodes, 2, "Count of sampled ancestor nodes should be 2.")
        self.assertEqual(ann_tr_2sas_no_spn_survives_no_max_age.origin_edge_length, 2.0, "Length of origin edge should be 2.0.")


    def test_node_counting_oneSA_no_spn_dies(self):
        """
        Test counting of extant, extinct and SA nodes

        Tree: origin + SA then dies before max_age
        Note that the way "execute_sample_ancestor" is implemented, when this event happens before a
        birth event (giving rise to the root), the origin instead undergoes a "fake birth" event, where
        a dummy node is made its child, and two children are added to the dummy node:
        
        (1) SA,
        (2) A node between the origin and the root (the "brosc" node)
    
        Then (2) stays alive and potentially undergoes a proper birth event (when the root would be created)
        or another ancestor sampling event.
        Here, the "brosc" actually dies before:
        
        (i) max_age and no other events, 
        (ii) there is no max_age

        To see it on icytree: ((sa1:0.0,brosc:0.2)dummy1:1.0)origin:0.0;
        """
        
        origin_node = Node(taxon=Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node.state = 0
        origin_node.alive = False
        origin_node.sampled = False
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        dummy_node = Node(taxon=Taxon(label="dummy1"), label="dummy1", edge_length=1.0)
        dummy_node.state = 0
        dummy_node.alive = False
        dummy_node.sampled = False
        dummy_node.is_sa = False
        dummy_node.is_sa_dummy_parent = True
        
        origin_node.add_child(dummy_node)
        
        # right child of dummy_node
        sa_node = Node(taxon=Taxon(label="sa1"), label="sa1", edge_length=0.0)
        sa_node.state = 0
        sa_node.alive = False
        sa_node.sampled = True
        sa_node.is_sa = True
        sa_node.is_sa_dummy_parent = False
        sa_node.is_sa_lineage = False

        # left child of dummy node
        brosc_node = Node(taxon=Taxon(label="brosc"), label="brosc", edge_length=0.2)
        brosc_node.state = 0
        brosc_node.alive = False
        brosc_node.sampled = False
        brosc_node.is_sa = False
        brosc_node.is_sa_dummy_parent = False
        brosc_node.is_sa_lineage = True
        
        dummy_node.add_child(sa_node)
        dummy_node.add_child(brosc_node)

        tr_sa_no_spn_dies = Tree(seed_node=origin_node)
        tr_sa_no_spn_dies.taxon_namespace.add_taxon(origin_node.taxon)
        tr_sa_no_spn_dies.taxon_namespace.add_taxon(dummy_node.taxon)
        tr_sa_no_spn_dies.taxon_namespace.add_taxon(brosc_node.taxon)
        tr_sa_no_spn_dies.taxon_namespace.add_taxon(sa_node.taxon)

        # debugging
        # print("tr_sa_no_spn_dies.seed_age = " + str(tr_sa_no_spn_dies.max_distance_from_root()))
        # print(tr_sa_no_spn_dies.as_string(schema="newick"))

        total_state_count = 1

        sa_global_time = 1.0
        time_to_sa_lineage_node = 0.2
        sa = pjsa.SampledAncestor(
            "sa1",
            "brosc",
            sa_global_time,
            time_to_lineage_node=time_to_sa_lineage_node)
        sa_lineage_dict = { "brosc": [sa] }
        
        max_age = 2.0

        ann_tr_sa_no_spn_dies_max_age = pjtr.AnnotatedTree(
            tr_sa_no_spn_dies,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        self.assertEqual(
            ann_tr_sa_no_spn_dies_max_age.n_extant_terminal_nodes,
            0,
            "Count of terminal extant nodes should be 0.")
        self.assertEqual(
            ann_tr_sa_no_spn_dies_max_age.n_extinct_terminal_nodes,
            1,
            "Count of terminal extinct nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_no_spn_dies_max_age.n_sa_nodes,
            1,
            "Count of sampled ancestor nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_no_spn_dies_max_age.origin_edge_length,
            1.2,
            "Length of origin edge should be 1.2.")

        tr_sa_no_spn_dies_no_max_age = Tree(seed_node=origin_node)
        tr_sa_no_spn_dies_no_max_age.taxon_namespace.add_taxon(origin_node.taxon)
        tr_sa_no_spn_dies_no_max_age.taxon_namespace.add_taxon(dummy_node.taxon)
        tr_sa_no_spn_dies_no_max_age.taxon_namespace.add_taxon(brosc_node.taxon)
        tr_sa_no_spn_dies_no_max_age.taxon_namespace.add_taxon(sa_node.taxon)

        ann_tr_sa_no_spn_dies_no_max_age = pjtr.AnnotatedTree(
            tr_sa_no_spn_dies_no_max_age,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        self.assertEqual(
            ann_tr_sa_no_spn_dies_no_max_age.n_extant_terminal_nodes,
            0,
            "Count of terminal extant nodes should be 0.")
        self.assertEqual(
            ann_tr_sa_no_spn_dies_no_max_age.n_extinct_terminal_nodes,
            1,
            "Count of terminal extinct nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_no_spn_dies_no_max_age.n_sa_nodes,
            1,
            "Count of sampled ancestor nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_no_spn_dies_no_max_age.origin_edge_length,
            1.2,
            "Length of origin edge should be 1.2.")


    def test_node_counting_twoSAs_no_spn_dies(self):
        """
        Test counting of extant, extinct and SA nodes

        Tree: origin + SA + SA then dies before max_age
        Note that the way "execute_sample_ancestor" is implemented, when this event happens before a
        birth event (giving rise to the root), the origin instead undergoes a "fake birth" event, where
        a dummy node is made its child, and two children are added to the dummy node:
        
        (1) SA,
        (2) A node between the origin and the root (the "brosc" node)
    
        Then (2) stays alive and potentially undergoes a proper birth event (when the root would be created)
        or another ancestor sampling event.
        Here, the "brosc" node undergoes another ancestor sampling event, and the resulting daughter "brosc"
        node actually dies before:
        
        (i) max_age and no other events, 
        (ii) there is no max_age

        To see it on icytree: ((sa1:0.0,(sa2:0.0,brosc:0.2)dummy2:0.5)dummy1:0.5)origin:0.0;
        """
        
        origin_node = Node(taxon=Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node.state = 0
        origin_node.alive = False
        origin_node.sampled = False
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        dummy_node1 = Node(taxon=Taxon(label="dummy1"), label="dummy1", edge_length=0.5)
        dummy_node1.state = 0
        dummy_node1.alive = False
        dummy_node1.sampled = False
        dummy_node1.is_sa = False
        dummy_node1.is_sa_dummy_parent = True
        
        origin_node.add_child(dummy_node1)
        
        # right child of dummy_node1 
        sa_node1 = Node(taxon=Taxon(label="sa1"), label="sa1", edge_length=0.0)
        sa_node1.state = 0
        sa_node1.alive = False
        sa_node1.sampled = True
        sa_node1.is_sa = True
        sa_node1.is_sa_dummy_parent = False
        sa_node1.is_sa_lineage = False

        # left child of dummy_node1
        dummy_node2 = Node(taxon=Taxon(label="dummy2"), label="dummy2", edge_length=0.5)
        dummy_node2.state = 0
        dummy_node2.alive = False
        dummy_node2.sampled = False
        dummy_node2.is_sa = False
        dummy_node2.is_sa_dummy_parent = True

        dummy_node1.add_child(sa_node1)
        dummy_node1.add_child(dummy_node2)
        
        # right child of dummy_node2
        sa_node2 = Node(taxon=Taxon(label="sa2"), label="sa2", edge_length=0.0)
        sa_node2.state = 0
        sa_node2.alive = False
        sa_node2.sampled = True
        sa_node2.is_sa = True
        sa_node2.is_sa_dummy_parent = False
        sa_node2.is_sa_lineage = False

        # left child of dummy_node2
        brosc_node = Node(taxon=Taxon(label="brosc"), label="brosc", edge_length=0.2)
        brosc_node.state = 0
        brosc_node.alive = False
        brosc_node.sampled = False
        brosc_node.is_sa = False
        brosc_node.is_sa_dummy_parent = False
        brosc_node.is_sa_lineage = True
        
        dummy_node2.add_child(sa_node2)
        dummy_node2.add_child(brosc_node)

        tr_2sas_no_spn_dies = Tree(seed_node=origin_node)
        tr_2sas_no_spn_dies.taxon_namespace.add_taxon(origin_node.taxon)
        tr_2sas_no_spn_dies.taxon_namespace.add_taxon(dummy_node1.taxon)
        tr_2sas_no_spn_dies.taxon_namespace.add_taxon(dummy_node2.taxon)
        tr_2sas_no_spn_dies.taxon_namespace.add_taxon(sa_node2.taxon)
        tr_2sas_no_spn_dies.taxon_namespace.add_taxon(brosc_node.taxon)
        
        # debugging
        # print("tr_2sas_no_spn_dies.seed_age = " + str(tr_2sas_no_spn_dies.max_distance_from_root()))
        # print(tr_2sas_no_spn_dies.as_string(schema="newick"))

        total_state_count = 1

        sa1_global_time = 0.5
        time_to_sa1_lineage_node = 0.5
        sa1 = pjsa.SampledAncestor(
            "sa1",
            "dummy_node2",
            sa1_global_time,
            time_to_lineage_node=time_to_sa1_lineage_node)

        sa2_global_time = 1.0
        time_to_sa2_lineage_node = 0.5
        sa2 = pjsa.SampledAncestor(
            "sa2",
            "brosc",
            sa2_global_time,
            time_to_lineage_node=time_to_sa2_lineage_node)

        sa_lineage_dict = { "brosc": [sa1, sa2] }
        
        max_age = 2.0

        ann_tr_2sas_no_spn_dies_max_age = pjtr.AnnotatedTree(
            tr_2sas_no_spn_dies,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        self.assertEqual(
            ann_tr_2sas_no_spn_dies_max_age.n_extant_terminal_nodes,
            0,
            "Count of terminal extant nodes should be 0.")
        self.assertEqual(
            ann_tr_2sas_no_spn_dies_max_age.n_extinct_terminal_nodes,
            1,
            "Count of terminal extinct nodes should be 1.")
        self.assertEqual(
            ann_tr_2sas_no_spn_dies_max_age.n_sa_nodes,
            2,
            "Count of sampled ancestor nodes should be 2.")
        self.assertEqual(
            ann_tr_2sas_no_spn_dies_max_age.origin_edge_length,
            1.2,
            "Length of origin edge should be 1.2.")

        tr_2sas_no_spn_dies_no_max_age = Tree(seed_node=origin_node)
        tr_2sas_no_spn_dies_no_max_age.taxon_namespace.add_taxon(origin_node.taxon)
        tr_2sas_no_spn_dies_no_max_age.taxon_namespace.add_taxon(dummy_node1.taxon)
        tr_2sas_no_spn_dies_no_max_age.taxon_namespace.add_taxon(dummy_node2.taxon)
        tr_2sas_no_spn_dies_no_max_age.taxon_namespace.add_taxon(sa_node2.taxon)
        tr_2sas_no_spn_dies_no_max_age.taxon_namespace.add_taxon(brosc_node.taxon)

        ann_tr_2sas_no_spn_dies_no_max_age = pjtr.AnnotatedTree(
            tr_2sas_no_spn_dies_no_max_age,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        self.assertEqual(
            ann_tr_2sas_no_spn_dies_no_max_age.n_extant_terminal_nodes,
            0,
            "Count of terminal extant nodes should be 0.")
        self.assertEqual(
            ann_tr_2sas_no_spn_dies_no_max_age.n_extinct_terminal_nodes,
            1,
            "Count of terminal extinct nodes should be 1.")
        self.assertEqual(
            ann_tr_2sas_no_spn_dies_no_max_age.n_sa_nodes,
            2,
            "Count of sampled ancestor nodes should be 2.")
        self.assertEqual(
            ann_tr_2sas_no_spn_dies_no_max_age.origin_edge_length,
            1.2,
            "Length of origin edge should be 1.2.")


    def test_node_counting_oneSA_with_root_survives(self):
        """
        Test counting of extant, extinct and SA nodes

        Tree: origin + SA + birth (root), then one lineage dies and the other survives until with max age.
        Note that the way "execute_sample_ancestor" is implemented, when this event happens before a
        birth event (giving rise to the root), the origin instead undergoes a "fake birth" event, where
        a dummy node is made its child, and two children are added to the dummy node:
        
        (1) SA,
        (2) A node between the origin and the root (the "brosc" node)
    
        Then (2) stays alive and potentially undergoes a proper birth event (when the root would be created)
        or another ancestor sampling event.
        Here, the "brosc" undergoes a birth event, and a root with two children are created.
        These two children then survive:
        
        (i) until the end of the process at max_age and no other events, 
        (ii) there is no max_age

        To see it on icytree: ((sa1:0.0,(sp1:0.25,sp2:0.5)root:0.5)dummy1:1.0)origin:0.0;
        """
        
        origin_node = Node(taxon=Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node.state = 0
        origin_node.alive = False
        origin_node.sampled = False
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        dummy_node = Node(taxon=Taxon(label="dummy1"), label="dummy1", edge_length=1.0)
        dummy_node.state = 0
        dummy_node.alive = False
        dummy_node.sampled = False
        dummy_node.is_sa = False
        dummy_node.is_sa_dummy_parent = True
        
        origin_node.add_child(dummy_node)
        
        # right child of dummy_node
        sa_node = Node(taxon=Taxon(label="sa1"), label="sa1", edge_length=0.0)
        sa_node.state = 0
        sa_node.alive = False
        sa_node.sampled = False
        sa_node.is_sa = True
        sa_node.is_sa_dummy_parent = False
        sa_node.is_sa_lineage = False

        # left child of dummy node
        root_node = Node(taxon=Taxon(label="root"), label="root", edge_length=0.5)
        root_node.state = 0
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = True
        
        dummy_node.add_child(sa_node)
        dummy_node.add_child(root_node)

        # left child of root node
        extant_sp1 = Node(taxon=Taxon(label="sp1"), label="sp1", edge_length=0.25)
        extant_sp1.state = 0
        extant_sp1.alive = False
        extant_sp1.sampled = False
        extant_sp1.is_sa = False
        extant_sp1.is_sa_dummy_parent = False
        extant_sp1.is_sa_lineage = False

        # right child of root node
        extant_sp2 = Node(taxon=Taxon(label="sp2"), label="sp2", edge_length=0.5)
        extant_sp2.state = 0
        extant_sp2.alive = True
        extant_sp2.sampled = True
        extant_sp2.is_sa = False
        extant_sp2.is_sa_dummy_parent = False
        extant_sp2.is_sa_lineage = False

        root_node.add_child(extant_sp1)
        root_node.add_child(extant_sp2)

        tr_sa_with_root_survives = Tree(seed_node=origin_node)
        tr_sa_with_root_survives.taxon_namespace.add_taxon(origin_node.taxon)
        tr_sa_with_root_survives.taxon_namespace.add_taxon(sa_node.taxon)
        tr_sa_with_root_survives.taxon_namespace.add_taxon(dummy_node.taxon)
        tr_sa_with_root_survives.taxon_namespace.add_taxon(root_node.taxon)
        tr_sa_with_root_survives.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_sa_with_root_survives.taxon_namespace.add_taxon(extant_sp2.taxon)

        # debugging
        # print("tr_sa_with_root_survives.seed_age = " + str(tr_sa_with_root_survives.max_distance_from_root()))
        # print(tr_sa_with_root_survives.as_string(schema="newick"))

        total_state_count = 1

        sa_global_time = 1.0
        time_to_sa_lineage_node = 0.5
        # TODO: when birth-event happens on brosc_node and it gets swapped with a root node,
        # we need to update the SAs older than the root to have "root" as their lineage nodes
        sa = pjsa.SampledAncestor(
            "sa1",
            "root",
            sa_global_time,
            time_to_lineage_node=time_to_sa_lineage_node)
        sa_lineage_dict = { "root": [sa] }
        
        max_age = 2.0

        ann_tr_sa_with_root_survives_max_age = pjtr.AnnotatedTree(
            tr_sa_with_root_survives,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        self.assertEqual(
            ann_tr_sa_with_root_survives_max_age.n_extant_terminal_nodes,
            1,
            "Count of terminal extant nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_with_root_survives_max_age.n_extinct_terminal_nodes,
            1,
            "Count of terminal extinct nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_with_root_survives_max_age.n_sa_nodes,
            1,
            "Count of sampled ancestor nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_with_root_survives_max_age.origin_edge_length,
            1.5,
            "Length of origin edge should be 1.5.")

        tr_sa_with_root_survives_no_max_age = Tree(seed_node=origin_node)
        tr_sa_with_root_survives_no_max_age.taxon_namespace.add_taxon(origin_node.taxon)
        tr_sa_with_root_survives_no_max_age.taxon_namespace.add_taxon(sa_node.taxon)
        tr_sa_with_root_survives_no_max_age.taxon_namespace.add_taxon(dummy_node.taxon)
        tr_sa_with_root_survives_no_max_age.taxon_namespace.add_taxon(root_node.taxon)
        tr_sa_with_root_survives_no_max_age.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_sa_with_root_survives_no_max_age.taxon_namespace.add_taxon(extant_sp2.taxon)

        ann_tr_sa_with_root_survives_no_max_age = pjtr.AnnotatedTree(
            tr_sa_with_root_survives_no_max_age,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        self.assertEqual(
            ann_tr_sa_with_root_survives_no_max_age.n_extant_terminal_nodes,
            1,
            "Count of terminal extant nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_with_root_survives_no_max_age.n_extinct_terminal_nodes,
            1,
            "Count of terminal extinct nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_with_root_survives_no_max_age.n_sa_nodes,
            1,
            "Count of sampled ancestor nodes should be 1.")
        self.assertEqual(
            ann_tr_sa_with_root_survives_no_max_age.origin_edge_length,
            1.5,
            "Length of origin edge should be 1.5.")

        

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
    # $ python3 tests/data/test_tree_annot_with_sas_from_origin.py
    # 
    # or
    #
    # $ python3 -m tests.data.test_tree_annot_with_sas_from_origin
    #
    # or 
    #
    # $ python3 -m unittest tests.data.test_tree_annot_with_sas_from_origin.TestAnnotateTreeWithSAsFromOrigin.test_node_counting_oneSA_no_spn_survives_max_age

    unittest.main()

