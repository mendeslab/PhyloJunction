from phylojunction.data.trait import DiscreteTrait
import unittest
import matplotlib
from dendropy import Tree, Node, Taxon

from pathlib import Path
DATA_TEST_INPUT_DIR = Path(__file__).parent.parent / "data" / "test_input"
TR_FILE_PATH_STR = str(DATA_TEST_INPUT_DIR / "tr.tre")
SMAPS_FILE_PATH_STR = str(DATA_TEST_INPUT_DIR / "smaps.tsv")
TIPDATA_FILE_PATH_STR = str(DATA_TEST_INPUT_DIR / "tip_data.tsv")

# pj imports
import phylojunction.data.tree as pjtr
import phylojunction.data.sampled_ancestor as pjsa
import phylojunction.data.attribute_transition as pjat
import phylojunction.interface.cmdbox.cmd_parse as cmdp

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestReconstructedTreePrint(unittest.TestCase):

    def test_one_lineage_left_rec_tr_plot(self) -> None:
        """Test plot objects of FBD tree with one surviving taxon.

        To see it on icytree:
        (sp2:2.0)origin_origin:0.0;
        """

        origin_node = Node(taxon=Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node.state = 0
        origin_node.annotations.add_bound_attribute("state")
        origin_node.index = 0
        origin_node.annotations.add_bound_attribute("index")
        origin_node.alive = False
        origin_node.sampled = False
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        root_node = Node(taxon=Taxon(label="root"), label="root", edge_length=1.0)
        root_node.state = 0
        root_node.annotations.add_bound_attribute("state")
        root_node.index = 1
        root_node.annotations.add_bound_attribute("index")
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        origin_node.add_child(root_node)

        # left child of root_node
        extinct_sp1 = Node(taxon=Taxon(label="sp1"),
                           label="sp1",
                           edge_length=0.5)
        extinct_sp1.state = 0
        extinct_sp1.annotations.add_bound_attribute("state")
        extinct_sp1.index = 2
        extinct_sp1.annotations.add_bound_attribute("index")
        extinct_sp1.alive = False
        extinct_sp1.sampled = False
        extinct_sp1.is_sa = False
        extinct_sp1.is_sa_dummy_parent = False
        extinct_sp1.is_sa_lineage = False

        # right child of root_node
        # left child of internal_node2
        extant_sp2 = Node(taxon=Taxon(label="sp2"),
                          label="sp2",
                          edge_length=1.0)
        extant_sp2.state = 0
        extant_sp2.annotations.add_bound_attribute("state")
        extant_sp2.index = 3
        extant_sp2.annotations.add_bound_attribute("index")
        extant_sp2.alive = True
        extant_sp2.sampled = True
        extant_sp2.is_sa = False
        extant_sp2.is_sa_dummy_parent = False
        extant_sp2.is_sa_lineage = False

        root_node.add_child(extinct_sp1)
        root_node.add_child(extant_sp2)

        # wrapping up tree
        tr_complete = Tree(seed_node=origin_node)
        tr_complete.taxon_namespace.add_taxon(root_node.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp2.taxon)

        max_age = 2.0
        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(1),
            start_at_origin=True,
            max_age=max_age,
            epsilon=1e-12)


        tr_rec = \
            ann_tr.extract_reconstructed_tree(plotting_overhead=True,
                                              require_obs_both_sides=False)

        # print(tr_rec.as_string(schema="newick"))

        fig = matplotlib.pyplot.figure()
        ax = fig.add_axes([0.25, 0.2, 0.5, 0.6])
        ax.patch.set_alpha(0.0)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        # plotting complete or rec tree
        sa_along_branches = True  # if False, cannot draw complete tree
        draw_reconstructed = True
        x_coords, y_coords = \
            pjtr.plot_ann_tree(ann_tr,
                               ax,
                               use_age=True,
                               sa_along_branches=sa_along_branches,
                               attr_of_interest="state",
                               draw_reconstructed=draw_reconstructed)
        matplotlib.pyplot.show()  # to see it (compare to baseline figs!)
        # new_svg_path = "baseline_figs/test_one_lineage_left_rec_tr_plot.svg"
        # matplotlib.pyplot.savefig(new_svg_path)

        self.assertEqual(x_coords, {'sp2': 0.0, 'origin': 2.0})
        self.assertEqual(y_coords, {'sp2': 2.0, 'origin': 2.0})

    def test_easy_fbd_rec_tr_sa_plot(self) -> None:
        """Test plot objects for easy reconstructed FBD tree.

        To see it on icytree:
        (sa1:0.0,(sp1:2.0,(sp3:1.0,sa2:0.0)dummy_node2:1.0)nd5_nd5:1.0)dummy_node1:0.0;
        """

        origin_node = Node(taxon=Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node.state = 0
        origin_node.annotations.add_bound_attribute("state")
        origin_node.index = 0
        origin_node.annotations.add_bound_attribute("index")
        origin_node.alive = False
        origin_node.sampled = False
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        root_node = Node(taxon=Taxon(label="root"), label="root", edge_length=1.0)
        root_node.state = 0
        root_node.annotations.add_bound_attribute("state")
        root_node.index = 1
        root_node.annotations.add_bound_attribute("index")
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        origin_node.add_child(root_node)

        # left child of root
        dummy_node1 = Node(taxon=Taxon(label="dummy_node1"),
                              label="dummy_node1",
                              edge_length=1.0)
        dummy_node1.state = 0
        dummy_node1.annotations.add_bound_attribute("state")
        dummy_node1.index = 2
        dummy_node1.annotations.add_bound_attribute("index")
        dummy_node1.alive = False
        dummy_node1.sampled = False
        dummy_node1.is_sa = False
        dummy_node1.is_sa_dummy_parent = True
        dummy_node1.is_sa_lineage = False

        # right child of root
        extinct_sp4 = Node(taxon=Taxon(label="sp4"),
                           label="sp4",
                           edge_length=1.0)
        extinct_sp4.state = 0
        extinct_sp4.annotations.add_bound_attribute("state")
        extinct_sp4.index = 3
        extinct_sp4.annotations.add_bound_attribute("index")
        extinct_sp4.alive = False
        extinct_sp4.sampled = False
        extinct_sp4.is_sa = False
        extinct_sp4.is_sa_dummy_parent = False
        extinct_sp4.is_sa_lineage = False

        root_node.add_child(dummy_node1)
        root_node.add_child(extinct_sp4)

        # left child of dummy_node1
        sa_node1 = Node(taxon=Taxon(label="sa1"),
                       label="sa1",
                       edge_length=0.0)
        sa_node1.state = 0
        sa_node1.annotations.add_bound_attribute("state")
        sa_node1.index = 4
        sa_node1.annotations.add_bound_attribute("index")
        sa_node1.alive = False
        sa_node1.sampled = True
        sa_node1.is_sa = True
        sa_node1.is_sa_dummy_parent = False
        sa_node1.is_sa_lineage = False

        # right child of dummy_node1
        internal_node1 = Node(taxon=Taxon(label="nd5"),
                              label="nd5",
                              edge_length=1.0)
        internal_node1.state = 0
        internal_node1.annotations.add_bound_attribute("state")
        internal_node1.index = 5
        internal_node1.annotations.add_bound_attribute("index")
        internal_node1.alive = False
        internal_node1.sampled = False
        internal_node1.is_sa = False
        internal_node1.is_sa_dummy_parent = False
        internal_node1.is_sa_lineage = True

        dummy_node1.add_child(sa_node1)
        dummy_node1.add_child(internal_node1)

        # left child of internal_node1
        internal_node2 = Node(taxon=Taxon(label="nd6"),
                              label="nd6",
                              edge_length=1.0)
        internal_node2.state = 0
        internal_node2.annotations.add_bound_attribute("state")
        internal_node2.index = 6
        internal_node2.annotations.add_bound_attribute("index")
        internal_node2.alive = False
        internal_node2.sampled = False
        internal_node2.is_sa = False
        internal_node2.is_sa_dummy_parent = False
        internal_node2.is_sa_lineage = False

        # right child of internal_node1
        dummy_node2 = Node(taxon=Taxon(label="dummy_node2"),
                           label="dummy_node2",
                           edge_length=1.0)
        dummy_node2.state = 0
        dummy_node2.annotations.add_bound_attribute("state")
        dummy_node2.index = 7
        dummy_node2.annotations.add_bound_attribute("index")
        dummy_node2.alive = False
        dummy_node2.sampled = False
        dummy_node2.is_sa = False
        dummy_node2.is_sa_dummy_parent = True
        dummy_node2.is_sa_lineage = False

        internal_node1.add_child(internal_node2)
        internal_node1.add_child(dummy_node2)

        # left child of internal_node2
        extant_sp1 = Node(taxon=Taxon(label="sp1"),
                          label="sp1",
                          edge_length=1.0)
        extant_sp1.state = 0
        extant_sp1.annotations.add_bound_attribute("state")
        extant_sp1.index = 8
        extant_sp1.annotations.add_bound_attribute("index")
        extant_sp1.alive = True
        extant_sp1.sampled = True
        extant_sp1.is_sa = False
        extant_sp1.is_sa_dummy_parent = False
        extant_sp1.is_sa_lineage = False

        # right child of internal_node2
        extinct_sp2 = Node(taxon=Taxon(label="sp2"),
                           label="sp2",
                           edge_length=0.5)
        extinct_sp2.state = 0
        extinct_sp2.annotations.add_bound_attribute("state")
        extinct_sp2.index = 9
        extinct_sp2.annotations.add_bound_attribute("index")
        extinct_sp2.alive = False
        extinct_sp2.sampled = False
        extinct_sp2.is_sa = False
        extinct_sp2.is_sa_dummy_parent = False
        extinct_sp2.is_sa_lineage = False

        internal_node2.add_child(extant_sp1)
        internal_node2.add_child(extinct_sp2)

        # left child of dummy_node2
        extant_sp3 = Node(taxon=Taxon(label="sp3"),
                          label="sp3",
                          edge_length=1.0)
        extant_sp3.state = 0
        extant_sp3.annotations.add_bound_attribute("state")
        extant_sp3.index = 10
        extant_sp3.annotations.add_bound_attribute("index")
        extant_sp3.alive = True
        extant_sp3.sampled = True
        extant_sp3.is_sa = False
        extant_sp3.is_sa_dummy_parent = False
        extant_sp3.is_sa_lineage = True

        # right child of dummy_node2
        sa_node2 = Node(taxon=Taxon(label="sa2"),
                        label="sa2",
                        edge_length=0.0)
        sa_node2.state = 0
        sa_node2.annotations.add_bound_attribute("state")
        sa_node2.index = 11
        sa_node2.annotations.add_bound_attribute("index")
        sa_node2.alive = False
        sa_node2.sampled = True
        sa_node2.is_sa = True
        sa_node2.is_sa_dummy_parent = False
        sa_node2.is_sa_lineage = False

        dummy_node2.add_child(extant_sp3)
        dummy_node2.add_child(sa_node2)

        # wrapping up tree
        tr_complete = Tree(seed_node=origin_node)
        tr_complete.taxon_namespace.add_taxon(origin_node.taxon)
        tr_complete.taxon_namespace.add_taxon(root_node.taxon)
        tr_complete.taxon_namespace.add_taxon(internal_node1.taxon)
        tr_complete.taxon_namespace.add_taxon(internal_node2.taxon)
        tr_complete.taxon_namespace.add_taxon(dummy_node1.taxon)
        tr_complete.taxon_namespace.add_taxon(dummy_node2.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp2.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp4.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp3.taxon)
        tr_complete.taxon_namespace.add_taxon(sa_node1.taxon)
        tr_complete.taxon_namespace.add_taxon(sa_node2.taxon)

        sa1 = pjsa.SampledAncestor(
            "sa1",
            "nd5",
            2.0,
            time_to_lineage_node=1.0)

        sa2 = pjsa.SampledAncestor(
            "sa2",
            "sp3",
            4.0,
            time_to_lineage_node=1.0)

        sa_lineage_dict = {"nd5": [sa1], "sp3": [sa2]}

        max_age = 5.0
        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(1),
            start_at_origin=True,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        tr_rec = \
            ann_tr.extract_reconstructed_tree(plotting_overhead=True,
                                              require_obs_both_sides=False)

        # print(ann_tr.tree_reconstructed.as_string(schema="newick"))

        fig = matplotlib.pyplot.figure()
        ax = fig.add_axes([0.25, 0.2, 0.5, 0.6])
        ax.patch.set_alpha(0.0)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        # plotting complete or rec tree
        sa_along_branches = True # if False, cannot draw complete tree
        draw_reconstructed = True
        x_coords, y_coords = \
            pjtr.plot_ann_tree(ann_tr,
                           ax,
                           use_age=True,
                           sa_along_branches=sa_along_branches,
                           attr_of_interest="state",
                           draw_reconstructed=draw_reconstructed)
        matplotlib.pyplot.show()  # to see it (compare to baseline figs!)
        # new_svg_path = "baseline_figs/test_easy_fbd_rec_tr_sa_plot.svg"
        # matplotlib.pyplot.savefig(new_svg_path)

        exp_x_coords = \
            {"sa1": 3.0, "sp1": 0.0, "sp3": 0.0, "sa2": 1.0,
             "dummy_node2": 1.0, "nd5": 2.0, "dummy_node1": 3.0}
        exp_y_coords = \
            {"sp3": 3.0, "sp1": 2.0, "dummy_node2": 3.0, "nd5": 2.5, "dummy_node1": 2.5}

        self.assertEqual(exp_x_coords, x_coords)
        self.assertEqual(exp_y_coords, y_coords)

    def test_harder_fbd_rec_tr_sa_plot(self) -> None:
        """Test plot for harder reconstructed FBD tree.

        Here, the FBD tree has two consecutive lineages that go extinct
        and we must update 'rec_tree_sa_lineage_dict'.

        To see it on icytree:
        (sa1:0.0,(sp4:1.0,sp5:1.0)nd7:3.0)dummy_node1:0.0;
        """

        # root of complete tree
        root_node = Node(taxon=Taxon(label="root"), label="root", edge_length=0.0)
        root_node.state = 0
        root_node.annotations.add_bound_attribute("state")
        root_node.index = 0
        root_node.annotations.add_bound_attribute("index")
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        # left child of root_node
        extinct_sp1 = Node(taxon=Taxon(label="sp1"),
                           label="sp1",
                           edge_length=1.0)
        extinct_sp1.state = 0
        extinct_sp1.annotations.add_bound_attribute("state")
        extinct_sp1.index = 1
        extinct_sp1.annotations.add_bound_attribute("index")
        extinct_sp1.alive = False
        extinct_sp1.sampled = False
        extinct_sp1.is_sa = False
        extinct_sp1.is_sa_dummy_parent = False
        extinct_sp1.is_sa_lineage = False

        # right child of root_node
        # root of rec tree
        dummy_node1 = Node(taxon=Taxon(label="dummy_node1"),
                           label="dummy_node1",
                           edge_length=1.0)
        dummy_node1.state = 0
        dummy_node1.annotations.add_bound_attribute("state")
        dummy_node1.index = 2
        dummy_node1.annotations.add_bound_attribute("index")
        dummy_node1.alive = False
        dummy_node1.sampled = False
        dummy_node1.is_sa = False
        dummy_node1.is_sa_dummy_parent = True
        dummy_node1.is_sa_lineage = False

        root_node.add_child(extinct_sp1)
        root_node.add_child(dummy_node1)

        # left child of dummy_node1
        sa_node1 = Node(taxon=Taxon(label="sa1"),
                        label="sa1",
                        edge_length=0.0)
        sa_node1.state = 0
        sa_node1.annotations.add_bound_attribute("state")
        sa_node1.index = 3
        sa_node1.annotations.add_bound_attribute("index")
        sa_node1.alive = False
        sa_node1.sampled = True
        sa_node1.is_sa = True
        sa_node1.is_sa_dummy_parent = False
        sa_node1.is_sa_lineage = False

        # right child of dummy_node1
        internal_node1 = Node(taxon=Taxon(label="nd5"),
                              label="nd5",
                              edge_length=1.0)
        internal_node1.state = 0
        internal_node1.annotations.add_bound_attribute("state")
        internal_node1.index = 4
        internal_node1.annotations.add_bound_attribute("index")
        internal_node1.alive = False
        internal_node1.sampled = False
        internal_node1.is_sa = False
        internal_node1.is_sa_dummy_parent = False
        internal_node1.is_sa_lineage = True

        dummy_node1.add_child(sa_node1)
        dummy_node1.add_child(internal_node1)

        # left child of internal_node1
        extinct_sp2 = Node(taxon=Taxon(label="sp2"),
                           label="sp2",
                           edge_length=1.0)
        extinct_sp2.state = 0
        extinct_sp2.annotations.add_bound_attribute("state")
        extinct_sp2.index = 5
        extinct_sp2.annotations.add_bound_attribute("index")
        extinct_sp2.alive = False
        extinct_sp2.sampled = False
        extinct_sp2.is_sa = False
        extinct_sp2.is_sa_dummy_parent = False
        extinct_sp2.is_sa_lineage = False

        # right child of internal_node1
        internal_node2 = Node(taxon=Taxon(label="nd6"),
                              label="nd6",
                              edge_length=1.0)
        internal_node2.state = 0
        internal_node2.annotations.add_bound_attribute("state")
        internal_node2.index = 6
        internal_node2.annotations.add_bound_attribute("index")
        internal_node2.alive = False
        internal_node2.sampled = False
        internal_node2.is_sa = False
        internal_node2.is_sa_dummy_parent = False
        internal_node2.is_sa_lineage = False

        internal_node1.add_child(extinct_sp2)
        internal_node1.add_child(internal_node2)

        # left child of internal_node2
        extinct_sp3 = Node(taxon=Taxon(label="sp3"),
                           label="sp3",
                           edge_length=1.0)
        extinct_sp3.state = 0
        extinct_sp3.annotations.add_bound_attribute("state")
        extinct_sp3.index = 7
        extinct_sp3.annotations.add_bound_attribute("index")
        extinct_sp3.alive = False
        extinct_sp3.sampled = False
        extinct_sp3.is_sa = False
        extinct_sp3.is_sa_dummy_parent = False
        extinct_sp3.is_sa_lineage = False

        # right child of internal_node2
        internal_node3 = Node(taxon=Taxon(label="nd7"),
                              label="nd7",
                              edge_length=1.0)
        internal_node3.state = 0
        internal_node3.annotations.add_bound_attribute("state")
        internal_node3.index = 8
        internal_node3.annotations.add_bound_attribute("index")
        internal_node3.alive = False
        internal_node3.sampled = False
        internal_node3.is_sa = False
        internal_node3.is_sa_dummy_parent = False
        internal_node3.is_sa_lineage = False

        internal_node2.add_child(extinct_sp3)
        internal_node2.add_child(internal_node3)

        # left child of internal_node3
        extant_sp4 = Node(taxon=Taxon(label="sp4"),
                          label="sp4",
                          edge_length=1.0)
        extant_sp4.state = 0
        extant_sp4.annotations.add_bound_attribute("state")
        extant_sp4.index = 9
        extant_sp4.annotations.add_bound_attribute("index")
        extant_sp4.alive = True
        extant_sp4.sampled = True
        extant_sp4.is_sa = False
        extant_sp4.is_sa_dummy_parent = False
        extant_sp4.is_sa_lineage = False

        # right child of internal_node3
        extant_sp5 = Node(taxon=Taxon(label="sp5"),
                          label="sp5",
                          edge_length=1.0)
        extant_sp5.state = 0
        extant_sp5.annotations.add_bound_attribute("state")
        extant_sp5.index = 10
        extant_sp5.annotations.add_bound_attribute("index")
        extant_sp5.alive = True
        extant_sp5.sampled = True
        extant_sp5.is_sa = False
        extant_sp5.is_sa_dummy_parent = False
        extant_sp5.is_sa_lineage = False

        internal_node3.add_child(extant_sp4)
        internal_node3.add_child(extant_sp5)

        # wrapping up tree
        tr_complete = Tree(seed_node=root_node)
        tr_complete.taxon_namespace.add_taxon(dummy_node1.taxon)
        tr_complete.taxon_namespace.add_taxon(sa_node1.taxon)
        tr_complete.taxon_namespace.add_taxon(internal_node1.taxon)
        tr_complete.taxon_namespace.add_taxon(internal_node2.taxon)
        tr_complete.taxon_namespace.add_taxon(internal_node3.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp2.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp3.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp4.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp5.taxon)

        sa1 = pjsa.SampledAncestor(
            "sa1",
            "nd5",
            1.0,
            time_to_lineage_node=1.0)

        sa_lineage_dict = {"nd5": [sa1]}

        max_age = 5.0
        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(1),
            start_at_origin=False,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        tr_rec = \
            ann_tr.extract_reconstructed_tree(plotting_overhead=True,
                                              require_obs_both_sides=False)

        # print(tr_rec.as_string(schema="newick"))

        fig = matplotlib.pyplot.figure()
        ax = fig.add_axes([0.25, 0.2, 0.5, 0.6])
        ax.patch.set_alpha(0.0)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        # plotting complete or rec tree
        sa_along_branches = True  # if False, cannot draw complete tree
        draw_reconstructed = True
        x_coords, y_coords = \
            pjtr.plot_ann_tree(ann_tr,
                               ax,
                               use_age=True,
                               sa_along_branches=sa_along_branches,
                               attr_of_interest="state",
                               draw_reconstructed=draw_reconstructed)
        matplotlib.pyplot.show()  # to see it (compare to baseline figs!)
        # new_svg_path = "baseline_figs/test_harder_fbd_rec_tr_sa_plot.svg"
        # matplotlib.pyplot.savefig(new_svg_path)


        exp_x_coords = {"sa1": 4.0, "sp4": 0.0, "sp5": 0.0, "nd7": 1.0, "dummy_node1": 4.0}
        exp_y_coords = {"sp5": 3.0, "sp4": 2.0, "nd7": 2.5, "dummy_node1": 2.5}

        self.assertEqual(exp_x_coords, x_coords)
        self.assertEqual(exp_y_coords, y_coords)

    def test_sa_tips_fbd_rec_tr_sa_plot(self) -> None:
        """Test plot for FBD tree with SA tip.

        To see it on icytree:
        (sa1:1.0,(sa2:1.0,(sp3:2.0,sp4:2.0)nd6:1.0)nd5:1.0)root:0.0;
        """

        # root of complete tree
        root_node = Node(taxon=Taxon(label="root"), label="root", edge_length=0.0)
        root_node.state = 0
        root_node.annotations.add_bound_attribute("state")
        root_node.index = 0
        root_node.annotations.add_bound_attribute("index")
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        # left child of root_node
        dummy_node1 = Node(taxon=Taxon(label="dummy_node1"),
                           label="dummy_node1",
                           edge_length=1.0)
        dummy_node1.state = 0
        dummy_node1.annotations.add_bound_attribute("state")
        dummy_node1.index = 1
        dummy_node1.annotations.add_bound_attribute("index")
        dummy_node1.alive = False
        dummy_node1.sampled = False
        dummy_node1.is_sa = False
        dummy_node1.is_sa_dummy_parent = True
        dummy_node1.is_sa_lineage = False

        # right child of root_node
        internal_node1 = Node(taxon=Taxon(label="nd5"),
                              label="nd5",
                              edge_length=1.0)
        internal_node1.state = 0
        internal_node1.annotations.add_bound_attribute("state")
        internal_node1.index = 2
        internal_node1.annotations.add_bound_attribute("index")
        internal_node1.alive = False
        internal_node1.sampled = False
        internal_node1.is_sa = False
        internal_node1.is_sa_dummy_parent = False
        internal_node1.is_sa_lineage = False

        root_node.add_child(dummy_node1)
        root_node.add_child(internal_node1)

        # left child of dummy_node1
        sa_node1 = Node(taxon=Taxon(label="sa1"),
                        label="sa1",
                        edge_length=0.0)
        sa_node1.state = 0
        sa_node1.annotations.add_bound_attribute("state")
        sa_node1.index = 3
        sa_node1.annotations.add_bound_attribute("index")
        sa_node1.alive = False
        sa_node1.sampled = True
        sa_node1.is_sa = True
        sa_node1.is_sa_dummy_parent = False
        sa_node1.is_sa_lineage = False

        # right child of dummy_node1
        extinct_sp1 = Node(taxon=Taxon(label="sp1"),
                           label="sp1",
                           edge_length=1.0)
        extinct_sp1.state = 0
        extinct_sp1.annotations.add_bound_attribute("state")
        extinct_sp1.index = 4
        extinct_sp1.annotations.add_bound_attribute("index")
        extinct_sp1.alive = False
        extinct_sp1.sampled = False
        extinct_sp1.is_sa = False
        extinct_sp1.is_sa_dummy_parent = False
        extinct_sp1.is_sa_lineage = True

        dummy_node1.add_child(sa_node1)
        dummy_node1.add_child(extinct_sp1)

        # left child of internal_node1
        dummy_node2 = Node(taxon=Taxon(label="dummy_node2"),
                          label="dummy_node2",
                          edge_length=1.0)
        dummy_node2.state = 0
        dummy_node2.annotations.add_bound_attribute("state")
        dummy_node2.index = 5
        dummy_node2.annotations.add_bound_attribute("index")
        dummy_node2.alive = False
        dummy_node2.sampled = False
        dummy_node2.is_sa = False
        dummy_node2.is_sa_dummy_parent = True
        dummy_node2.is_sa_lineage = False

        # right child of internal_node1
        internal_node2 = Node(taxon=Taxon(label="nd6"),
                              label="nd6",
                              edge_length=1.0)
        internal_node2.state = 0
        internal_node2.annotations.add_bound_attribute("state")
        internal_node2.index = 6
        internal_node2.annotations.add_bound_attribute("index")
        internal_node2.alive = False
        internal_node2.sampled = False
        internal_node2.is_sa = False
        internal_node2.is_sa_dummy_parent = False
        internal_node2.is_sa_lineage = False

        internal_node1.add_child(dummy_node2)
        internal_node1.add_child(internal_node2)

        # left child of dummy_node2
        sa_node2 = Node(taxon=Taxon(label="sa2"),
                        label="sa2",
                        edge_length=0.0)
        sa_node2.state = 0
        sa_node2.annotations.add_bound_attribute("state")
        sa_node2.index = 7
        sa_node2.annotations.add_bound_attribute("index")
        sa_node2.alive = False
        sa_node2.sampled = True
        sa_node2.is_sa = True
        sa_node2.is_sa_dummy_parent = False
        sa_node2.is_sa_lineage = False

        # right child of dummy_node2
        extinct_sp2 = Node(taxon=Taxon(label="sp2"),
                           label="sp2",
                           edge_length=1.0)
        extinct_sp2.state = 0
        extinct_sp2.annotations.add_bound_attribute("state")
        extinct_sp2.index = 8
        extinct_sp2.annotations.add_bound_attribute("index")
        extinct_sp2.alive = False
        extinct_sp2.sampled = False
        extinct_sp2.is_sa = False
        extinct_sp2.is_sa_dummy_parent = False
        extinct_sp2.is_sa_lineage = True

        dummy_node2.add_child(sa_node2)
        dummy_node2.add_child(extinct_sp2)

        # left child of internal_node2
        extant_sp3 = Node(taxon=Taxon(label="sp3"),
                          label="sp3",
                          edge_length=2.0)
        extant_sp3.state = 0
        extant_sp3.annotations.add_bound_attribute("state")
        extant_sp3.index = 9
        extant_sp3.annotations.add_bound_attribute("index")
        extant_sp3.alive = True
        extant_sp3.sampled = True
        extant_sp3.is_sa = False
        extant_sp3.is_sa_dummy_parent = False
        extant_sp3.is_sa_lineage = False

        # right child of internal_node2
        extant_sp4 = Node(taxon=Taxon(label="sp4"),
                          label="sp4",
                          edge_length=2.0)
        extant_sp4.state = 0
        extant_sp4.annotations.add_bound_attribute("state")
        extant_sp4.index = 10
        extant_sp4.annotations.add_bound_attribute("index")
        extant_sp4.alive = True
        extant_sp4.sampled = True
        extant_sp4.is_sa = False
        extant_sp4.is_sa_dummy_parent = False
        extant_sp4.is_sa_lineage = False

        internal_node2.add_child(extant_sp3)
        internal_node2.add_child(extant_sp4)

        # wrapping up tree
        tr_complete = Tree(seed_node=root_node)
        tr_complete.taxon_namespace.add_taxon(dummy_node1.taxon)
        tr_complete.taxon_namespace.add_taxon(dummy_node2.taxon)
        tr_complete.taxon_namespace.add_taxon(sa_node1.taxon)
        tr_complete.taxon_namespace.add_taxon(sa_node2.taxon)
        tr_complete.taxon_namespace.add_taxon(internal_node1.taxon)
        tr_complete.taxon_namespace.add_taxon(internal_node2.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp2.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp3.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp4.taxon)

        sa1 = pjsa.SampledAncestor(
            "sa1",
            "sp1",
            1.0,
            time_to_lineage_node=1.0)

        sa2 = pjsa.SampledAncestor(
            "sa2",
            "sp2",
            2.0,
            time_to_lineage_node=1.0)

        sa_lineage_dict = {"sp1": [sa1], "sp2": [sa2]}

        max_age = 4.0
        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(1),
            start_at_origin=False,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        # print(ann_tr.tree.as_string(schema="newick"))

        tr_rec = \
            ann_tr.extract_reconstructed_tree(plotting_overhead=True,
                                              require_obs_both_sides=False)

        # print(ann_tr.tree_reconstructed.as_string(schema="newick"))

        fig = matplotlib.pyplot.figure()
        ax = fig.add_axes([0.25, 0.2, 0.5, 0.6])
        ax.patch.set_alpha(0.0)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        # plotting complete or rec tree
        sa_along_branches = True  # if False, cannot draw complete tree
        draw_reconstructed = True
        x_coords, y_coords = \
            pjtr.plot_ann_tree(ann_tr,
                               ax,
                               use_age=True,
                               sa_along_branches=sa_along_branches,
                               attr_of_interest="state",
                               draw_reconstructed=draw_reconstructed)
        matplotlib.pyplot.show()  # to see it (compare to baseline figs!)
        # new_svg_path = "baseline_figs/test_sa_tips_fbd_rec_tr_sa_plot.svg"
        # matplotlib.pyplot.savefig(new_svg_path)

        exp_x_coords = {"sa1": 3.0, "sa2": 2.0, "sp3": 0.0, "sp4": 0.0, "nd6": 2.0, "nd5": 3.0, "root": 4.0}
        exp_y_coords = {"sp4": 5.0, "sp3": 4.0, "sa2": 3.0, "sa1": 2.0, "nd6": 4.5, "nd5": 3.75, "root": 2.875}

        self.assertEqual(exp_x_coords, x_coords)
        self.assertEqual(exp_y_coords, y_coords)

    def test_sa_followed_by_sa_tip_fbd_rec_tr_plot(self) -> None:
        """Test plot for FBD tree with SA followed by tip.

        To see it on icytree:
        ((sa1:0.0,sa2:1.0)dummy_node1:1.0,(sa3:0.0,sp2:3.0)dummy_node3:1.0)root:0.0;
        """

        # root of complete tree
        root_node = Node(taxon=Taxon(label="root"), label="root", edge_length=0.0)
        root_node.state = 0
        root_node.annotations.add_bound_attribute("state")
        root_node.index = 0
        root_node.annotations.add_bound_attribute("index")
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        # left child of root_node
        dummy_node1 = Node(taxon=Taxon(label="dummy_node1"),
                           label="dummy_node1",
                           edge_length=1.0)
        dummy_node1.state = 0
        dummy_node1.annotations.add_bound_attribute("state")
        dummy_node1.index = 1
        dummy_node1.annotations.add_bound_attribute("index")
        dummy_node1.alive = False
        dummy_node1.sampled = False
        dummy_node1.is_sa = False
        dummy_node1.is_sa_dummy_parent = True
        dummy_node1.is_sa_lineage = False

        # right child of root_node
        dummy_node3 = Node(taxon=Taxon(label="dummy_node3"),
                           label="dummy_node3",
                           edge_length=1.0)
        dummy_node3.state = 0
        dummy_node3.annotations.add_bound_attribute("state")
        dummy_node3.index = 2
        dummy_node3.annotations.add_bound_attribute("index")
        dummy_node3.alive = False
        dummy_node3.sampled = False
        dummy_node3.is_sa = False
        dummy_node3.is_sa_dummy_parent = True
        dummy_node3.is_sa_lineage = False

        root_node.add_child(dummy_node1)
        root_node.add_child(dummy_node3)

        # left child of dummy_node1
        sa_node1 = Node(taxon=Taxon(label="sa1"),
                        label="sa1",
                        edge_length=0.0)
        sa_node1.state = 0
        sa_node1.annotations.add_bound_attribute("state")
        sa_node1.index = 3
        sa_node1.annotations.add_bound_attribute("index")
        sa_node1.alive = False
        sa_node1.sampled = True
        sa_node1.is_sa = True
        sa_node1.is_sa_dummy_parent = False
        sa_node1.is_sa_lineage = False

        # right child of dummy_node1
        dummy_node2 = Node(taxon=Taxon(label="dummy_node2"),
                           label="dummy_node2",
                           edge_length=1.0)
        dummy_node2.state = 0
        dummy_node2.annotations.add_bound_attribute("state")
        dummy_node2.index = 4
        dummy_node2.annotations.add_bound_attribute("index")
        dummy_node2.alive = False
        dummy_node2.sampled = False
        dummy_node2.is_sa = False
        dummy_node2.is_sa_dummy_parent = True
        dummy_node2.is_sa_lineage = False

        dummy_node1.add_child(sa_node1)
        dummy_node1.add_child(dummy_node2)

        # left child of dummy_node2
        sa_node2 = Node(taxon=Taxon(label="sa2"),
                        label="sa2",
                        edge_length=0.0)
        sa_node2.state = 0
        sa_node2.annotations.add_bound_attribute("state")
        sa_node2.index = 7
        sa_node2.annotations.add_bound_attribute("index")
        sa_node2.alive = False
        sa_node2.sampled = True
        sa_node2.is_sa = True
        sa_node2.is_sa_dummy_parent = False
        sa_node2.is_sa_lineage = False

        # right child of dummy_node2
        extinct_sp1 = Node(taxon=Taxon(label="sp1"),
                           label="sp1",
                           edge_length=1.0)
        extinct_sp1.state = 0
        extinct_sp1.annotations.add_bound_attribute("state")
        extinct_sp1.index = 8
        extinct_sp1.annotations.add_bound_attribute("index")
        extinct_sp1.alive = False
        extinct_sp1.sampled = False
        extinct_sp1.is_sa = False
        extinct_sp1.is_sa_dummy_parent = False
        extinct_sp1.is_sa_lineage = True

        dummy_node2.add_child(sa_node2)
        dummy_node2.add_child(extinct_sp1)

        # left child of dummy_node3
        sa_node3 = Node(taxon=Taxon(label="sa3"),
                        label="sa3",
                        edge_length=0.0)
        sa_node3.state = 0
        sa_node3.annotations.add_bound_attribute("state")
        sa_node3.index = 5
        sa_node3.annotations.add_bound_attribute("index")
        sa_node3.alive = False
        sa_node3.sampled = True
        sa_node3.is_sa = True
        sa_node3.is_sa_dummy_parent = False
        sa_node3.is_sa_lineage = False

        # right child of dummy_node3
        extant_sp2 = Node(taxon=Taxon(label="sp2"),
                          label="sp2",
                          edge_length=3.0)
        extant_sp2.state = 0
        extant_sp2.annotations.add_bound_attribute("state")
        extant_sp2.index = 6
        extant_sp2.annotations.add_bound_attribute("index")
        extant_sp2.alive = True
        extant_sp2.sampled = True
        extant_sp2.is_sa = False
        extant_sp2.is_sa_dummy_parent = False
        extant_sp2.is_sa_lineage = False

        dummy_node3.add_child(sa_node3)
        dummy_node3.add_child(extant_sp2)

        # wrapping up tree
        tr_complete = Tree(seed_node=root_node)
        tr_complete.taxon_namespace.add_taxon(dummy_node1.taxon)
        tr_complete.taxon_namespace.add_taxon(dummy_node2.taxon)
        tr_complete.taxon_namespace.add_taxon(dummy_node3.taxon)
        tr_complete.taxon_namespace.add_taxon(sa_node1.taxon)
        tr_complete.taxon_namespace.add_taxon(sa_node2.taxon)
        tr_complete.taxon_namespace.add_taxon(sa_node3.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp2.taxon)

        sa1 = pjsa.SampledAncestor(
            "sa1",
            "sp1",
            1.0,
            time_to_lineage_node=2.0)

        sa2 = pjsa.SampledAncestor(
            "sa2",
            "sp1",
            2.0,
            time_to_lineage_node=1.0)

        sa3 = pjsa.SampledAncestor(
            "sa3",
            "sp2",
            1.0,
            time_to_lineage_node=1.0)

        sa_lineage_dict = {"sp1": [sa1, sa2], "sp2": [sa3]}

        max_age = 4.0
        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(1),
            start_at_origin=False,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        # print(ann_tr.tree.as_string(schema="newick"))

        tr_rec = \
            ann_tr.extract_reconstructed_tree(plotting_overhead=True,
                                              require_obs_both_sides=False)

        # print(ann_tr.tree_reconstructed.as_string(schema="newick"))

        fig = matplotlib.pyplot.figure()
        ax = fig.add_axes([0.25, 0.2, 0.5, 0.6])
        ax.patch.set_alpha(0.0)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        # plotting complete or rec tree
        sa_along_branches = True  # if False, cannot draw complete tree
        draw_reconstructed = True
        x_coords, y_coords = \
            pjtr.plot_ann_tree(ann_tr,
                               ax,
                               use_age=True,
                               sa_along_branches=sa_along_branches,
                               attr_of_interest="state",
                               draw_reconstructed=draw_reconstructed)
        matplotlib.pyplot.show()  # to see it (compare to baseline figs!)
        # new_svg_path = "baseline_figs/test_sa_followed_by_sa_tip_fbd_rec_tr_plot.svg"
        # matplotlib.pyplot.savefig(new_svg_path)

        exp_x_coords = {"sa1": 3.0, "sa2": 2.0, "dummy_node1": 3.0, "sa3": 3.0, "sp2": 0.0, "dummy_node3": 3.0, "root": 4.0}
        exp_y_coords = {"sp2": 2.0, "sa2": 1.0, "dummy_node1": 1.0, "dummy_node3": 2.0, "root": 1.5}

        self.assertEqual(exp_x_coords, x_coords)
        self.assertEqual(exp_y_coords, y_coords)

    def test_easy_geosse_rec_tr(self) -> None:
        """Test plot for GeoSSE tree.

        Only transition is on root edge, which should not have
        any bearing on reconstructed tree's plot (should be all
        the same state 0.
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
        root_node.state = 0  # A
        root_node.annotations.add_bound_attribute("state")
        root_node.index = 1
        root_node.annotations.add_bound_attribute("index")
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        origin_node.add_child(root_node)

        # left child of root_node
        extant_sp1 = Node(taxon=Taxon(label="sp1"),
                          label="sp1",
                          edge_length=2.0)
        extant_sp1.state = 0  # A
        extant_sp1.annotations.add_bound_attribute("state")
        extant_sp1.index = 2
        extant_sp1.annotations.add_bound_attribute("index")
        extant_sp1.alive = True
        extant_sp1.sampled = True
        extant_sp1.is_sa = False
        extant_sp1.is_sa_dummy_parent = False
        extant_sp1.is_sa_lineage = False

        # right child of root_node
        internal_node1 = Node(taxon=Taxon(label="nd3"),
                              label="nd3",
                              edge_length=1.0)
        internal_node1.state = 0  # A
        internal_node1.annotations.add_bound_attribute("state")
        internal_node1.index = 3
        internal_node1.annotations.add_bound_attribute("index")
        internal_node1.alive = False
        internal_node1.sampled = False
        internal_node1.is_sa = False
        internal_node1.is_sa_dummy_parent = False
        internal_node1.is_sa_lineage = False

        root_node.add_child(extant_sp1)
        root_node.add_child(internal_node1)

        # left child of internal_node1
        extinct_sp2 = Node(taxon=Taxon(label="sp2"),
                           label="sp2",
                           edge_length=0.5)
        extinct_sp2.state = 0  # A
        extinct_sp2.annotations.add_bound_attribute("state")
        extinct_sp2.index = 4
        extinct_sp2.annotations.add_bound_attribute("index")
        extinct_sp2.alive = False
        extinct_sp2.sampled = False
        extinct_sp2.is_sa = False
        extinct_sp2.is_sa_dummy_parent = False
        extinct_sp2.is_sa_lineage = False

        # right child of internal_node1
        extant_sp3 = Node(taxon=Taxon(label="sp3"),
                          label="sp3",
                          edge_length=1.0)
        extant_sp3.state = 0  # A
        extant_sp3.annotations.add_bound_attribute("state")
        extant_sp3.index = 5
        extant_sp3.annotations.add_bound_attribute("index")
        extant_sp3.alive = True
        extant_sp3.sampled = True
        extant_sp3.is_sa = False
        extant_sp3.is_sa_dummy_parent = False
        extant_sp3.is_sa_lineage = False

        internal_node1.add_child(extinct_sp2)
        internal_node1.add_child(extant_sp3)

        # wrapping up tree
        tr_complete = Tree(seed_node=origin_node)
        tr_complete.taxon_namespace.add_taxon(origin_node.taxon)
        tr_complete.taxon_namespace.add_taxon(root_node.taxon)
        tr_complete.taxon_namespace.add_taxon(internal_node1.taxon)
        tr_complete.taxon_namespace.add_taxon(extinct_sp2.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp3.taxon)

        at1 = pjat.AttributeTransition("state", "root", 0.5, 2, 0)

        at_dict = { "root": [at1] }

        total_state_count = 3
        max_age = 2.0
        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(total_state_count),
            at_dict=at_dict,
            start_at_origin=True,
            max_age=max_age,
            epsilon=1e-12)

        ann_tr.populate_nd_attr_dict("state")

        tr_rec = \
            ann_tr.extract_reconstructed_tree(plotting_overhead=True,
                                              require_obs_both_sides=False)

        # print(ann_tr.tree_reconstructed.as_string(schema="newick"))

        fig = matplotlib.pyplot.figure()
        ax = fig.add_axes([0.25, 0.2, 0.5, 0.6])
        ax.patch.set_alpha(0.0)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        # plotting complete or rec tree
        draw_reconstructed = True
        x_coords, y_coords = \
            pjtr.plot_ann_tree(ann_tr,
                               ax,
                               use_age=False,
                               sa_along_branches=False,
                               attr_of_interest="state",
                               draw_reconstructed=draw_reconstructed)
        matplotlib.pyplot.show()  # to see it (compare to baseline figs!)
        # new_svg_path = "baseline_figs/test_easy_geosse_rec_tr.svg"
        # matplotlib.pyplot.savefig(new_svg_path)

        self.assertEqual({'sp1': 2.0, 'sp3': 2.0, 'root': 0.0}, x_coords)
        self.assertEqual({'sp3': 3.0, 'sp1': 2.0, 'root': 2.5}, y_coords)

    def test_harder_geosse_rec_tr(self) -> None:
        """Test plot for GeoSSE tree with 4 extinct tips."""

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

        total_state_count = 3
        max_age = 5.0
        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            DiscreteTrait(total_state_count),
            at_dict=at_dict,
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

        fig = matplotlib.pyplot.figure()
        ax = fig.add_axes([0.25, 0.2, 0.5, 0.6])
        ax.patch.set_alpha(0.0)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        # plotting complete or rec tree
        draw_reconstructed = True
        x_coords, y_coords = \
            pjtr.plot_ann_tree(ann_tr,
                           ax,
                           use_age=False,
                           sa_along_branches=False,
                           attr_of_interest="state",
                           draw_reconstructed=draw_reconstructed)
        matplotlib.pyplot.show()  # to see it (compare to baseline figs!)
        # new_svg_path = "baseline_figs/test_harder_geosse_rec_tr.svg"
        # matplotlib.pyplot.savefig(new_svg_path)

        self.assertEqual({'sp2': 4.0, 'sp5': 4.0, 'root': 0.0},
                         x_coords)
        self.assertEqual({'sp5': 3.0, 'sp2': 2.0, 'root': 2.5},
                         y_coords)

    def test_read_smaps_rec_tr(self) -> None:
        """Test stochastic-map reading function and tree plotting.

        The tree here is much larger (10 taxa), assumed to be
        an observed tree (i.e., effectively a 'reconstructed' tree)
        with multiple anagenetic and cladogenetic attribute transitions.
        """

        read_smap_str = \
            ('trs <- read_tree(n=1, file_path=\"' +
             TR_FILE_PATH_STR + '\",'
             'node_name_attr="index")\n'
             'mapped_trs := map_attr(tree=trs, fun="smap", '
             'maps_file_path=\"' +
             SMAPS_FILE_PATH_STR + '\",'
             'tip_attr_file_path=\"' +
             TIPDATA_FILE_PATH_STR + '\",'
             'attr_name="state", n_regions=2, geo="true")')

        dag_obj = cmdp.script2dag(read_smap_str, in_pj_file=False)

        fig = matplotlib.pyplot.figure()

        ax = fig.add_axes([0.25, 0.2, 0.5, 0.6])
        ax.patch.set_alpha(0.0)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        for node_name, node_dag in dag_obj.name_node_dict.items():
            if node_name == "mapped_trs":
                i = 0  # 1
                ann_tr = node_dag.value[i]

                x_coords, y_coords = \
                    pjtr.plot_ann_tree(ann_tr,
                                       ax,
                                       use_age=False,
                                       start_at_origin=False,
                                       sa_along_branches=False,
                                       attr_of_interest="state")

        matplotlib.pyplot.show()

        self.assertEqual(
            x_coords,
            { 'nd19': 1.0, 'nd20': 1.0, 'nd17': 0.9000279732601156,
              'nd18': 1.0, 'nd7': 0.8511350440846258, 'nd13': 1.0,
              'nd14': 1.0, 'nd8': 0.7890570597123967, 'nd1': 0.2996978898812708,
              'nd15': 1.0, 'nd16': 1.0,
              'nd9': 0.8082580225200156, 'nd21': 1.0, 'nd22': 1.0,
              'nd10': 0.9451562520235127, 'nd6': 0.3666853234617953,
              'nd12': 1.0, 'nd2': 0.20093516010407542, 'nd0': 0.0 }
        )

        self.assertEqual(
            y_coords,
            { 'nd12': 11.0, 'nd22': 10.0, 'nd21': 9.0, 'nd16': 8.0,
              'nd15': 7.0, 'nd14': 6.0, 'nd13': 5.0, 'nd18': 4.0,
              'nd20': 3.0, 'nd19': 2.0, 'nd17': 2.5, 'nd7': 3.25, 'nd8': 5.5,
              'nd1': 4.375, 'nd9': 7.5, 'nd10': 9.5,'nd6': 8.5, 'nd2': 9.75,
              'nd0': 7.0625}
        )


if __name__ == "__main__":
    # $ python3.11 tests/data/test_tree_plotting.py
    #
    # or
    #
    # $ python3.11 -m tests.data.test_tree_plotting
    #
    # or
    #
    # $ python3.11 -m unittest tests.data.test_tree_plotting.TestReconstructedTreePrint.test_one_lineage_left_rec_tr_plot

    unittest.main()

class TestColorbarLayout(unittest.TestCase):
    # Colorbar cleanup must preserve automatic layout without changing manual CLI positions.
    # Existing tree images do not cover repeated resizing or axes participation in layout.
    def test_colorbar_cleanup_preserves_layout(self):
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        import numpy as np

        for managed in (False, True):
            with self.subTest(managed=managed):
                fig = Figure(figsize=(8, 5), constrained_layout=managed)
                FigureCanvasAgg(fig)
                ax = fig.add_subplot(111) if managed else fig.add_axes([.075, .25, .6, .7])
                fig.canvas.draw()
                original = ax.get_position().bounds
                for width in (8, 12, 8):
                    fig.set_size_inches(width, 5)
                    for _ in range(3):
                        marker = ax.scatter([0, 1], [0, 1], c=[0, 1])
                        ax._pj_trait_position = ax.get_position(original=True).frozen()
                        ax._pj_trait_colorbar = fig.colorbar(marker, ax=ax, label='trait')
                        fig.canvas.draw()
                        self.assertEqual(len(fig.axes), 2)
                        pjtr.clear_trait_colorbar(ax)
                        marker.remove()
                        fig.canvas.draw()
                        self.assertEqual(len(fig.axes), 1)
                        if managed:
                            self.assertTrue(ax.get_in_layout())
                            self.assertGreater(ax.get_position().width, .8)
                        else:
                            np.testing.assert_allclose(ax.get_position().bounds, original)
