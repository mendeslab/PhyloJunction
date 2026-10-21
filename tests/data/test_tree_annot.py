from phylojunction.data.trait import DiscreteTrait
import unittest
import copy
from matplotlib import pyplot as plt
from dendropy import Tree, Node, Taxon

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"

# pj imports
import phylojunction.data.tree as pjtr

class TestAnnotateTree(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        # 
        # Four trees to be read as newick strings
        #
        # rootedgeless_tr_str = "((nd6:2.5,nd7:1.0)nd3:1.0,(nd4:1.0,nd5:2.5)nd2:1.0)root:0.0;"
        # rootedge_tr_str = "(((nd6:3.0,nd7:1.0)nd3:1.0,(nd4:1.0,nd5:3.0)nd2:1.0)root:2.0)origin:0.0;"
        # rootedgeless_no_spn_living_tr_str = "(brosc:1.0)origin:0.0;" # TODO: deprecate, replace with just origin node
        # rootedgeless_no_spn_dead_tr_str = "(brosc:0.5)origin:0.0;" # TODO: deprecate, replace with just origin node

        # tr_root = Tree.get(data=rootedgeless_tr_str, schema="newick")
        # tr_origin = Tree.get(data=rootedge_tr_str, schema="newick")
        # tr_no_spn = Tree.get(data=rootedgeless_no_spn_living_tr_str, schema="newick")
        # tr_no_spn_dead = Tree.get(data=rootedgeless_no_spn_dead_tr_str, schema="newick")

        #
        # Tree 1 (building by hand):
        #
        # TODO: replace this one with an origin node that is alive and no root
        origin_node1 = Node(taxon=Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node1.state = 0
        origin_node1.alive = False
        origin_node1.sampled = False
        origin_node1.is_sa = False
        brosc_node1 = Node(taxon=Taxon(label="brosc"), label="brosc", edge_length=1.0)
        brosc_node1.state = 0
        brosc_node1.alive = True
        brosc_node1.sampled = True
        brosc_node1.is_sa = False
        origin_node1.add_child(brosc_node1)
        tr_no_spn_built = Tree(seed_node=origin_node1)

        #
        # Tree 2 (building by hand):
        #
        # TODO: replace this one with an origin node that is not alive, and that
        # has an edge length < max age (need to add code that sets that edge length
        # correctly)
        origin_node2 = Node(taxon=Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node2.state = 0
        origin_node2.alive = False
        origin_node2.sampled = False
        origin_node2.is_sa = False
        brosc_node2 = Node(taxon=Taxon(label="brosc"), label="brosc", edge_length=0.5)
        brosc_node2.state = 0
        brosc_node2.alive = False
        brosc_node2.sampled = False
        brosc_node2.is_sa = False
        origin_node2.add_child(brosc_node2)
        tr_no_spn_dead_built = Tree(seed_node=origin_node2)

        total_state_count = 1

        # trees read in as newick
        # cls.tree_root = pjtr.AnnotatedTree(tr_root, total_state_count, start_at_origin=False, epsilon=1e-12)
        # cls.tree_origin = pjtr.AnnotatedTree(tr_origin, total_state_count, start_at_origin=True, epsilon=1e-12)
        # cls.tree_no_spn = pjtr.AnnotatedTree(tr_no_spn, total_state_count, start_at_origin=True, max_age=1.0, epsilon=1e-12)
        # cls.tree_no_spn_dead = pjtr.AnnotatedTree(tr_no_spn_dead, total_state_count, start_at_origin=True, max_age=1.0, epsilon=1e-12)
        # trees built by hand
        cls.tree_no_spn_built = pjtr.AnnotatedTree(tr_no_spn_built, DiscreteTrait(total_state_count), start_at_origin=True, max_age=1.0, epsilon=1e-12)
        cls.tree_no_spn_dead_built = pjtr.AnnotatedTree(tr_no_spn_dead_built, DiscreteTrait(total_state_count), start_at_origin=True, max_age=1.0, epsilon=1e-12)

        # for testing the parsing of states into AnnotatedTree member dictionary
        origin_node3 = Node(taxon=Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node3.state = 0
        origin_node3.alive = False
        origin_node3.sampled = False
        origin_node3.is_sa = False
        root_node3 = Node(taxon=Taxon(label="root"), label="root", edge_length=0.5)
        root_node3.state = 0
        root_node3.alive = False
        root_node3.sampled = False
        root_node3.is_sa = False
        origin_node3.add_child(root_node3)
        child_left = Node(taxon=Taxon(label="nd1"), label="nd1", edge_length=1.0)
        child_left.state = 0
        child_left.alive = True
        child_left.sampled = True
        child_left.is_sa = False
        root_node3.add_child(child_left)
        child_right = Node(taxon=Taxon(label="nd2"), label="nd2", edge_length=1.0)
        child_right.state = 1
        child_right.alive = True
        child_right.sampled = True
        child_right.is_sa = False
        root_node3.add_child(child_right)
        bif_tr_two_states = Tree(seed_node=origin_node3)

        total_state_count = 2

        cls.bifurcating_tree_two_states = pjtr.AnnotatedTree(bif_tr_two_states, DiscreteTrait(total_state_count), start_at_origin=True, epsilon=1e-12)

    # Renaming storage must not change tree results; ordinary state-only tests cannot detect
    # hard-coded attribute accesses. Retire only if equivalent rename coverage replaces this.
    def test_named_discrete_trait(self):
        original = copy.deepcopy(self.bifurcating_tree_two_states)
        for node in original.tree:
            node.is_sa_dummy_parent = node.is_sa_lineage = False
        raw = copy.deepcopy(original.tree)
        for node in raw:
            # AnnotatedTree removes internal taxa during initialization; restore raw input.
            raw.taxon_namespace.add_taxon(node.taxon)
            node.habitat = node.state
            del node.state
        renamed = pjtr.AnnotatedTree(raw, DiscreteTrait(2, name="habitat"), start_at_origin=True)
        self.assertEqual(renamed.state_count_dict, original.state_count_dict)
        self.assertEqual(renamed.get_stats_dict(), original.get_stats_dict())
        for nexus in (False, True):
            self.assertEqual(renamed.get_taxon_states_str(nexus=nexus),
                             original.get_taxon_states_str(nexus=nexus))
        for reconstructed in (False, True):
            fig, ax = plt.subplots()
            try:
                renamed.plot_node(ax, draw_reconstructed=reconstructed)
                colors = [line.get_colors().copy() for line in ax.collections]
                ax.clear()
                original.plot_node(ax, draw_reconstructed=reconstructed)
                self.assertEqual([c.tolist() for c in colors],
                                 [line.get_colors().tolist() for line in ax.collections])
            finally:
                plt.close(fig)
        self.assertEqual(renamed.rec_str(), original.rec_str())
        self.assertTrue(all(hasattr(nd, "habitat") and not hasattr(nd, "state")
                            for nd in renamed.tree_reconstructed))

    def test_root_origin_age_rootedge_length(self):
        """
        Test annotation of root/origin age and root edge
        """
        
        # tree starting at root
        # self.assertEqual(self.tree_root.origin_edge_length, 0.0, "No root edge: should be 0.0.")
        # self.assertIsNone(self.tree_root.origin_age, "Origin age should be 'none'.")
        # self.assertEqual(self.tree_root.root_age, 3.5, "Root age should be 3.5.")

        # tree starting at origin
        # self.assertEqual(self.tree_origin.origin_edge_length, 2.0, "Root edge should have length 2.0.")
        # self.assertEqual(self.tree_origin.origin_age, 6.0, "Origin age should be 6.0.")
        # self.assertEqual(self.tree_origin.root_age, 4.0, "Root age should be 4.0.")
        # self.assertEqual(self.tree_no_spn.origin_edge_length, 1.0, "Root edge should have length 1.0.")
        # self.assertEqual(self.tree_no_spn_dead.origin_edge_length, 0.5, "Root edge should have length 0.5.")
        self.assertEqual(self.tree_no_spn_built.origin_edge_length, 1.0, "Root edge should have length 1.0.")
        self.assertEqual(self.tree_no_spn_dead_built.origin_edge_length, 0.5, "Root edge should have length 0.5.")

        # tree starting at origin
        # self.assertFalse(self.tree_no_spn.tree_died, msg="Tree made all the way to maximum age (1.0), so it should not have been marked as dead.")
        # self.assertTrue(self.tree_no_spn_dead.tree_died, msg="Tree did not make all the way to maximum age (1.0), so it should have been marked as dead.")
        self.assertFalse(self.tree_no_spn_built.tree_died, msg="Tree made all the way to maximum age (1.0), so it should not have been marked as dead.")
        self.assertTrue(self.tree_no_spn_dead_built.tree_died, msg="Tree did not make all the way to maximum age (1.0), so it should have been marked as dead.")

    def test_node_counting(self):
        """
        Test counting of extant and extinct nodes
        """
        
        # tree starting at root
        # extant nodes
        # self.assertEqual(self.tree_root.n_extant_obs_nodes, 2, "Count of observable extant nodes should be 2.")

        # extant nodes label check
        # self.assertSequenceEqual(self.tree_root.extant_obs_nodes_labels, ("nd6", "nd5"), "Labels should be (\"nd6\", \"nd5\").")

        # extinct nodes
        # self.assertEqual(self.tree_root.n_extinct_obs_nodes, 2, "Count of observable extinct nodes should be 2.")
        
        # extinct nodes label check
        # self.assertSequenceEqual(self.tree_root.extinct_obs_nodes_labels, ("nd7", "nd4"), "Labels should be (\"nd7\", \"nd4\").")

        # ----------------------- #

        # tree starting at origin and either surviving without spn, or dying without speciation
        # extant nodes
        self.assertEqual(self.tree_no_spn_built.n_extant_terminal_nodes, 1, "Count of observable terminal nodes should be 1.")
        self.assertEqual(self.tree_no_spn_built.extant_terminal_nodes_labels, ("brosc",), "Labels should be (\"brosc\").")

        # extinct nodes
        self.assertEqual(self.tree_no_spn_dead_built.n_extinct_terminal_nodes, 1, "Count of terminal extinct nodes should be 1.")
        self.assertEqual(self.tree_no_spn_dead_built.extinct_terminal_nodes_labels, ("brosc",), "Labels should be (\"brosc\").")
        
        # extinct node label check
        # self.assertSequenceEqual(self.tree_origin.extinct_obs_nodes_labels, ("nd7", "nd4"), "Labels should be (\"nd7\", \"nd4\").")
        # self.assertSequenceEqual(self.tree_no_spn_dead.extinct_obs_nodes_labels, ("root",), "Labels should be (\"root\").") # must have comma
        
        # self.assertSequenceEqual(self.tree_no_spn_dead_built.extinct_obs_nodes_labels, ("brosc",), "Labels should be (\"brosc\").") # must have comma

    # def test_node_attr_dict_population(self):
    #     """
    #     Test counting of extant and extinct nodes
    #     """

    #     self.bifurcating_tree_two_states.populate_nd_attr_dict(["state"])
    #     self.bifurcating_tree_two_states.get_gcf(node_attr="state")

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
    # $ python3 tests/data/test_tree_annot.py
    # 
    # or
    #
    # $ python3 -m tests.data.test_tree_annot
    #
    # or 
    #
    # $ python3 -m unittest tests.data.test_tree_annot.TestAnnotateTree.test_root_origin_age_rootedge_length

    unittest.main()
