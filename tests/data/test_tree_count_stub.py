import unittest
from dendropy import Tree, Node, Taxon

# pj imports
import phylojunction.data.tree as pjtr


class TestTreeCountStub(unittest.TestCase):

    def test_count_stubs(self):
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

        # upper extinct child of nd6 (sampled ancestor)
        sa_nd = Node(
            taxon=Taxon(label="sa"),
            label="sa",
            edge_length=0.0
        )
        sa_nd.state = 0
        sa_nd.alive = False
        sa_nd.sampled = True
        sa_nd.is_sa = True
        sa_nd.is_sa_dummy_parent = False
        sa_nd.is_sa_lineage = False

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

        nd6.add_child(sa_nd)
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
        tr_complete.taxon_namespace.add_taxon(sa_nd.taxon)
        tr_complete.taxon_namespace.add_taxon(nd7.taxon)
        tr_complete.taxon_namespace.add_taxon(nd12.taxon)
        tr_complete.taxon_namespace.add_taxon(nd8.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp4.taxon)
        tr_complete.taxon_namespace.add_taxon(nd13.taxon)

        total_state_count = 1

        max_age = 6.0

        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            total_state_count,
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

        # debugging: checking reconstructed tree is OK
        # print(tr_rec.as_string(schema="newick"))

        ann_tr.populate_n_stubs_on_branches_of_rec_tr()

        expected_stub_counts = {
            "root": 0,
            "nd10": 1,
            "nd6": 0,
            "sa": 0,
            "nd7": 1,
            "nd12": 0,
            "nd13": 1
        }

        for nd in tr_rec.preorder_node_iter():
            self.assertEqual(nd.n_stub, expected_stub_counts[nd.label])


if __name__ == "__main__":
    # If you want to run this as a standalone from PhyloJunction/
    # on the terminal, remember to add "src/phylojunction" to
    # PYTHONPATH (system variable), or to set it if it does not
    # exist -- don't forget to export it!
    #
    # Then you can do:
    # $ python3 tests/data/test_tree_count_stub.py
    #
    # or
    #
    # $ python3 -m tests.data.test_tree_count_stub
    #
    # or
    #
    # $ python3 -m unittest tests.data.test_tree_count_stub.TestTreeCountStub.test_count_stubs

    unittest.main()

