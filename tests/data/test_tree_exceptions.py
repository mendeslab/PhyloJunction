import unittest
from dendropy import Tree, Node, Taxon

# pj imports #
import phylojunction.data.tree as pjtr
import phylojunction.data.sampled_ancestor as pjsa
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"

class TestTreeExceptions(unittest.TestCase):
    
    def test_alive_annotation_exceptions(self):
        """
        To see it on icytree:
            ((sa1:0.0,sp1:0.4):1.0,(sp2:0.8,sp3:0.8):0.75):0.0;
        """

        root_node = Node(
            taxon=Taxon(label="root"), label="root", edge_length=0.0
        )
        root_node.state = 0
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        # left child of root node after ancestor sampling happens
        # on who would have been the left child ("sp1")
        dummy_node = Node(
            taxon=Taxon(label="dummy1"), label="dummy1", edge_length=1.0
        )
        dummy_node.state = 0
        dummy_node.alive = False
        dummy_node.sampled = False
        dummy_node.is_sa = False
        dummy_node.is_sa_dummy_parent = True
        dummy_node.is_sa_lineage = False

        # right child of root
        int_node = Node(
            taxon=Taxon(label="nd1"), label="nd1", edge_length=0.75
        )
        int_node.state = 0
        int_node.alive = False
        int_node.sampled = False
        int_node.is_sa = False
        int_node.is_sa_dummy_parent = False
        int_node.is_sa_lineage = False

        extinct_sp2 = Node(
            taxon=Taxon(label="sp2"), label="sp2", edge_length=0.8
            )
        extinct_sp2.state = 0
        extinct_sp2.alive = False
        extinct_sp2.sampled = False
        extinct_sp2.is_sa = False
        extinct_sp2.is_sa_dummy_parent = False
        extinct_sp2.is_sa_lineage = False

        extant_sp3 = Node(
            taxon=Taxon(label="sp3"), label="sp3", edge_length=0.8
        )
        extant_sp3.state = 0
        extant_sp3.alive = False
        extant_sp3.sampled = False
        extant_sp3.is_sa = False
        extant_sp3.is_sa_dummy_parent = False
        extant_sp3.is_sa_lineage = False

        int_node.add_child(extinct_sp2)
        int_node.add_child(extant_sp3)
        
        # right child of dummy_node
        sa_node = Node(
            taxon=Taxon(label="sa1"), label="sa1", edge_length=0.0
        )
        sa_node.state = 0
        sa_node.alive = False
        sa_node.sampled = True
        sa_node.is_sa = True
        sa_node.is_sa_dummy_parent = False
        sa_node.is_sa_lineage = False

        # left child of dummy node
        # left child of root node
        extant_sp1 = Node(
            taxon=Taxon(label="sp1"), label="sp1", edge_length=0.5
        )
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
            "sa1", "sp1", sa_global_time,
            time_to_lineage_node=time_to_sa_lineage_node
        )
        sa_lineage_dict = { "sp1": [sa] }
        
        max_age = 2.0

        with self.assertRaises(ec.AnnotatedTreeIncorrectAnnotationError) as exc:
            ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            total_state_count,
            start_at_origin=False,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)
        
        self.assertEqual(str(exc.exception),
            ("Some aspect of the tree annotation was wrong. Taxon had "
             "non-maximal age, but had '.alive == True'. This is not "
             "allowed. Exiting..."))


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
    # $ python3 tests/data/test_tree_exceptions.py
    # 
    # or
    #
    # $ python3 -m tests.data.test_tree_exceptions
    #
    # or 
    #
    # $ python3 -m unittest tests.data.test_tree_exceptions.TestTreeExceptions.test_alive_annotation_exceptions

    unittest.main()
        
    