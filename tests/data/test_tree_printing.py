from phylojunction.data.trait import DiscreteTrait
import unittest
from dendropy import Tree, Node, Taxon

# pj imports
import phylojunction.data.tree as pjtr
import phylojunction.data.sampled_ancestor as pjsa

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"

class TestTreePrinting(unittest.TestCase):

    def test_nexus_printing(self):
        """Test if node labels look OK when Nexus printing
        
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
        root_node.is_sa_lineage = False
        
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
        time_to_sa_lineage_node = 1.0
        
        sa = pjsa.SampledAncestor("sa1", "root", sa_global_time, time_to_lineage_node=time_to_sa_lineage_node)
        sa_lineage_dict = { "root": [sa] }
        
        max_age = 2.0

        ann_tr_sa_with_root_survives_max_age = pjtr.AnnotatedTree(
            tr_sa_with_root_survives,
            DiscreteTrait(total_state_count),
            start_at_origin=True,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        expected_nexus_string = \
            "#NEXUS\n\nBEGIN TAXA;\n    DIMENSIONS NTAX=3;\n    TAXLABELS\n        sa1\n        sp1\n        sp2\n  ;\nEND;\n\nBEGIN TREES;\n" + \
            "    TREE 1 = ((sa1:0.0,(sp1:0.25,sp2:0.5)root:0.5)dummy1:1.0)origin:0.0;\nEND;\n\n"
        
        # debugging
        # print(expected_nexus_string)
        # print(repr(ann_tr_sa_with_root_survives_max_age.tree.as_string(schema="nexus", suppress_internal_taxon_labels=True)))

        self.assertEqual(
            ann_tr_sa_with_root_survives_max_age.tree.as_string(schema="nexus", suppress_internal_taxon_labels=True),
            expected_nexus_string)


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
    # $ python3 tests/data/test_tree_printing.py
    # 
    # or
    #
    # $ python3 -m tests.data.test_tree_printing
    #
    # or 
    #
    # $ python3 -m unittest tests.data.test_tree_printing.TestTreePrinting.test_nexus_printing

    unittest.main()