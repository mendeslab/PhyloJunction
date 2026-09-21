from phylojunction.data.trait import DiscreteTrait
import io
import unittest
import matplotlib
from tabulate import tabulate  # type: ignore
from dendropy import Node, Taxon, Tree

# pj imports
import phylojunction.data.tree as pjtr
import phylojunction.data.attribute_transition as pjat
import phylojunction.pgm.pgm as pgm
import phylojunction.readwrite.pj_write as pjw

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestSmapWrite(unittest.TestCase):

    def test_make_smap_str(self) -> None:
        """Test stoch mappings are correctly produced for two trees."""

        # to see plot, uncomment plotting command inside function
        # (sp2:4.0[&state=2,index=8],sp5:4.0[&state=2,index=11])root:0.0[&state=2,index=1];
        def build_tree1() -> pjtr.AnnotatedTree:
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
                                                 to_state2=1,
                                                 at_speciation=True)
            clado_at2 = pjat.AttributeTransition("state",
                                                 subtending_node_label="nd7",
                                                 global_time=3.0,
                                                 from_state=2,
                                                 to_state=0,
                                                 to_state2=1,
                                                 at_speciation=True)
            clado_at3 = pjat.AttributeTransition("state",
                                                 subtending_node_label="nd6",
                                                 global_time=2.0,
                                                 from_state=2,
                                                 to_state=0,
                                                 to_state2=1,
                                                 at_speciation=True)
            clado_at4 = pjat.AttributeTransition("state",
                                                 subtending_node_label="nd8",
                                                 global_time=2.0,
                                                 from_state=2,
                                                 to_state=0,
                                                 to_state2=1,
                                                 at_speciation=True)
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
            #         suppress_annotations=False,
            #         suppress_internal_taxon_labels=True,
            #         suppress_internal_node_labels=False)
            # print(ann_tr_str)

            # updates rec tree-related members in 'ann_tr'
            tr_rec = \
                ann_tr.extract_reconstructed_tree(
                    plotting_overhead=True,
                    require_obs_both_sides=False)

            rec_tr_str = \
                ann_tr.tree_reconstructed.as_string(
                    schema="newick",
                    suppress_annotations=False,
                    suppress_internal_taxon_labels=True,
                    suppress_internal_node_labels=False)
            # print(rec_tr_str)

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
            _, _ = \
                pjtr.plot_ann_tree(ann_tr,
                                   ax,
                                   use_age=False,
                                   sa_along_branches=False,
                                   attr_of_interest="state",
                                   draw_reconstructed=draw_reconstructed)
            # matplotlib.pyplot.show()  # to see it (compare to baseline figs!)
            # new_png_path = "tree1.png"
            # matplotlib.pyplot.savefig(new_png_path)

            return ann_tr

        # to see plot, uncomment plotting command inside function
        # (sp2:4.0[&state=2,index=8],((sp4:2.0[&state=1,index=10],sp5:2.0[&state=0,index=11])nd7:1.0[&state=2,index=5],sp6:3.0[&state=0,index=9])nd8:1.0[&state=2,index=3])root:0.0[&state=2,index=1];
        def build_tree2() -> pjtr.AnnotatedTree:
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
            extant_sp6 = Node(taxon=Taxon(label="sp6"),
                               label="sp6",
                               edge_length=3.0)
            extant_sp6.state = 0  # A
            extant_sp6.annotations.add_bound_attribute("state")
            extant_sp6.index = 9
            extant_sp6.annotations.add_bound_attribute("index")
            extant_sp6.alive = True
            extant_sp6.sampled = True
            extant_sp6.is_sa = False
            extant_sp6.is_sa_dummy_parent = False
            extant_sp6.is_sa_lineage = False

            # left child of internal_node4
            extant_sp4 = Node(taxon=Taxon(label="sp4"),
                               label="sp4",
                               edge_length=2.0)
            extant_sp4.state = 1  # B
            extant_sp4.annotations.add_bound_attribute("state")
            extant_sp4.index = 10
            extant_sp4.annotations.add_bound_attribute("index")
            extant_sp4.alive = True
            extant_sp4.sampled = True
            extant_sp4.is_sa = False
            extant_sp4.is_sa_dummy_parent = False
            extant_sp4.is_sa_lineage = False

            # right child of internal_node4
            extant_sp5 = Node(taxon=Taxon(label="sp5"),
                              label="sp5",
                              edge_length=2.0)
            extant_sp5.state = 0  # A
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

            internal_node4.add_child(extant_sp4)
            internal_node4.add_child(extant_sp5)

            internal_node3.add_child(internal_node4)  # 'nd7'
            internal_node3.add_child(extant_sp6)

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
            tr_complete.taxon_namespace.add_taxon(extant_sp4.taxon)
            tr_complete.taxon_namespace.add_taxon(extant_sp6.taxon)
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
            at4 = pjat.AttributeTransition("state", "sp4", 3.0, 2, 1)
            at5 = pjat.AttributeTransition("state", "sp5", 4.0, 2, 0)
            at_dict = {
                "nd5": [at1_1_1, at1_1_2],
                "nd7": [at1_1_3, at1_1_4],
                "sp3": [at1_2_1],
                "sp6": [at1_2_2],
                "sp2": [at2, at3],
                "sp4": [at4],
                "sp5": [at5]
            }

            # internal_node2
            clado_at1 = pjat.AttributeTransition("state",
                                                 subtending_node_label="nd5",
                                                 global_time=3.0,
                                                 from_state=2,
                                                 to_state=0,
                                                 to_state2=1,
                                                 at_speciation=True)
            clado_at2 = pjat.AttributeTransition("state",
                                                 subtending_node_label="nd7",
                                                 global_time=3.0,
                                                 from_state=2,
                                                 to_state=2,
                                                 to_state2=1,
                                                 at_speciation=True)
            clado_at3 = pjat.AttributeTransition("state",
                                                 subtending_node_label="nd6",
                                                 global_time=2.0,
                                                 from_state=2,
                                                 to_state=0,
                                                 to_state2=1,
                                                 at_speciation=True)
            clado_at4 = pjat.AttributeTransition("state",
                                                 subtending_node_label="nd8",
                                                 global_time=2.0,
                                                 from_state=2,
                                                 to_state=0,
                                                 to_state2=1,
                                                 at_speciation=True)
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
            #         suppress_annotations=False,
            #         suppress_internal_taxon_labels=True,
            #         suppress_internal_node_labels=False)
            # print(ann_tr_str)

            # updates rec tree-related members in 'ann_tr'
            # tr_rec = \
            #     ann_tr.extract_reconstructed_tree(
            #         plotting_overhead=True,
            #         require_obs_both_sides=False)

            # rec_tr_str = \
            #     ann_tr.tree_reconstructed.as_string(
            #         schema="newick",
            #         suppress_annotations=False,
            #         suppress_internal_taxon_labels=True,
            #         suppress_internal_node_labels=False)
            # print(rec_tr_str)

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
            # draw_reconstructed = True
            # _, _ = \
            #     pjtr.plot_ann_tree(ann_tr,
            #                        ax,
            #                        use_age=False,
            #                        sa_along_branches=False,
            #                        attr_of_interest="state",
            #                        draw_reconstructed=draw_reconstructed)
            # matplotlib.pyplot.show()  # to see it (compare to baseline figs!)
            # new_png_path = "tree2.png"
            # matplotlib.pyplot.savefig(new_png_path)

            return ann_tr

        ann_tr1 = build_tree1()
        # print("Done building tree 1\n")

        ann_tr2 = build_tree2()
        # print("Done building tree 2")

        node_val1 = [ann_tr1, ann_tr2]
        node_val2 = [ann_tr1, ann_tr2]

        dag_obj = pgm.DirectedAcyclicGraph()

        # one sample, two replicates
        stoch_node1 = pgm.StochasticNodeDAG(
            "tr1",
            sample_size=1,
            replicate_size=2,
            value=node_val1,
            clamped=True)

        # two samples, each with one replicate
        stoch_node2 = pgm.StochasticNodeDAG(
            "tr2",
            sample_size=2,
            replicate_size=1,
            value=node_val2,
            clamped=True)

        dag_obj.add_node(stoch_node1)
        dag_obj.add_node(stoch_node2)

        all_nds_df_list_dict, all_nds_df_content_str_list_dict = \
            pjw.prep_trees_rb_smap_dfs(dag_obj, ["tr1", "tr2"], "state")

        # writing to file handle
        rb_smap_dfs_outfile_tr1 = io.StringIO()
        rb_smap_dfs_outfile_tr2_1 = io.StringIO()
        rb_smap_dfs_outfile_tr2_2 = io.StringIO()

        # call the smap parsing method inside
        pjw.write_data_df(rb_smap_dfs_outfile_tr1, all_nds_df_list_dict["tr1"][0], format="tsv")
        pjw.write_data_df(rb_smap_dfs_outfile_tr2_1, all_nds_df_list_dict["tr2"][0], format="tsv")
        pjw.write_data_df(rb_smap_dfs_outfile_tr2_2, all_nds_df_list_dict["tr2"][1], format="tsv")

        # checking content
        rb_smap_dfs_outfile_tr1.seek(0)  # move to the start of the file handle
        smaps_str_tr1 = rb_smap_dfs_outfile_tr1.read()

        # debugging
        # print(smaps_str_tr1)

        rb_smap_dfs_outfile_tr2_1.seek(0)
        smaps_str_tr2_1 = rb_smap_dfs_outfile_tr2_1.read()
        rb_smap_dfs_outfile_tr2_2.seek(0)
        smaps_str_tr2_2 = rb_smap_dfs_outfile_tr2_2.read()

        # looking at file names
        # pjw.dump_trees_rb_smap_dfs("./", dag_obj, ["tr1"], "state")
        # pjw.dump_trees_rb_smap_dfs("./", dag_obj, ["tr1", "tr2"], "state")

        # printing the dataframes
        # for dag_node_name, sample_df_list in all_nds_df_list_dict.items():
        #     for sample_df in sample_df_list:
        #         print(tabulate(sample_df,
        #                        headers=sample_df.head(),
        #                        tablefmt="pretty",
        #                        showindex=False))
        #         print("\n\n")

        # printing dataframes' contents
        # for dag_node_name, sample_df_content_list in all_nds_df_content_str_list_dict.items():
        #     print("size of sample_df_content_list", len(sample_df_content_list))
        #     for sample_df_content in sample_df_content_list:
        #         print("\n".join("\t".join(l) for l in sample_df_content))
        #         print("\n\n\n")

        exp_content_tr1 = \
            ("iteration\tnode_index\tbranch_start_time\tbranch_end_time\t"
             "start_state\tend_state\ttransition_time\ttransition_type\t"
             "parent_index\tchild1_index\tchild2_index\n") + \
            "1\t1\t4.0\t4.0\t2\t2\tNA\tno_change\tNA\t8\t11\n" + \
            "1\t8\t4.0\t0.0\t2\t1\t3.0\tanagenetic\t1\tNA\tNA\n" + \
            "1\t8\t4.0\t0.0\t1\t2\t2.5\tanagenetic\t1\tNA\tNA\n" + \
            "1\t8\t4.0\t0.0\t2\t0\t2.0\tanagenetic\t1\tNA\tNA\n" + \
            "1\t8\t4.0\t0.0\t0\t2\t1.0\tanagenetic\t1\tNA\tNA\n" + \
            "1\t11\t4.0\t0.0\t2\t1\t3.0\tanagenetic\t1\tNA\tNA\n" + \
            "1\t11\t4.0\t0.0\t1\t2\t2.5\tanagenetic\t1\tNA\tNA\n" + \
            "1\t11\t4.0\t0.0\t2\t0\t2.0\tanagenetic\t1\tNA\tNA\n" + \
            "1\t11\t4.0\t0.0\t0\t2\t1.0\tanagenetic\t1\tNA\tNA\n" + \
            "2\t1\t4.0\t4.0\t2\t2\tNA\tno_change\tNA\t8\t3\n" + \
            "2\t8\t4.0\t0.0\t2\t1\t3.0\tanagenetic\t1\tNA\tNA\n" + \
            "2\t8\t4.0\t0.0\t1\t2\t2.5\tanagenetic\t1\tNA\tNA\n" + \
            "2\t8\t4.0\t0.0\t2\t0\t2.0\tanagenetic\t1\tNA\tNA\n" + \
            "2\t8\t4.0\t0.0\t0\t2\t1.0\tanagenetic\t1\tNA\tNA\n" + \
            "2\t3\t4.0\t3.0\t2\t0\t3.0\tcladogenetic\t1\t5\t9\n" + \
            "2\t3\t4.0\t3.0\t2\t1\t3.0\tcladogenetic\t1\t5\t9\n" + \
            "2\t3\t4.0\t3.0\t2\t2\tNA\tno_change\t1\t5\t9\n" + \
            "2\t5\t3.0\t2.0\t2\t1\t2.0\tcladogenetic\t3\t10\t11\n" + \
            "2\t5\t3.0\t2.0\t1\t2\t2.5\tanagenetic\t3\t10\t11\n" + \
            "2\t10\t2.0\t0.0\t1\t1\tNA\tno_change\t5\tNA\tNA\n" + \
            "2\t11\t2.0\t0.0\t2\t0\t1.0\tanagenetic\t5\tNA\tNA\n" + \
            "2\t9\t3.0\t0.0\t0\t0\tNA\tno_change\t3\tNA\tNA\n"

        self.assertEqual(smaps_str_tr1, exp_content_tr1)

        exp_content_tr2_1 = \
            ("iteration\tnode_index\tbranch_start_time\tbranch_end_time\t"
             "start_state\tend_state\ttransition_time\ttransition_type\t"
             "parent_index\tchild1_index\tchild2_index\n") + \
            "1\t1\t4.0\t4.0\t2\t2\tNA\tno_change\tNA\t8\t11\n" + \
            "1\t8\t4.0\t0.0\t2\t1\t3.0\tanagenetic\t1\tNA\tNA\n" + \
            "1\t8\t4.0\t0.0\t1\t2\t2.5\tanagenetic\t1\tNA\tNA\n" + \
            "1\t8\t4.0\t0.0\t2\t0\t2.0\tanagenetic\t1\tNA\tNA\n" + \
            "1\t8\t4.0\t0.0\t0\t2\t1.0\tanagenetic\t1\tNA\tNA\n" + \
            "1\t11\t4.0\t0.0\t2\t1\t3.0\tanagenetic\t1\tNA\tNA\n" + \
            "1\t11\t4.0\t0.0\t1\t2\t2.5\tanagenetic\t1\tNA\tNA\n" + \
            "1\t11\t4.0\t0.0\t2\t0\t2.0\tanagenetic\t1\tNA\tNA\n" + \
            "1\t11\t4.0\t0.0\t0\t2\t1.0\tanagenetic\t1\tNA\tNA\n"

        self.assertEqual(smaps_str_tr2_1, exp_content_tr2_1)

        exp_content_tr2_2 = \
            ("iteration\tnode_index\tbranch_start_time\tbranch_end_time\t"
             "start_state\tend_state\ttransition_time\ttransition_type\t"
             "parent_index\tchild1_index\tchild2_index\n") + \
            "1\t1\t4.0\t4.0\t2\t2\tNA\tno_change\tNA\t8\t3\n" + \
            "1\t8\t4.0\t0.0\t2\t1\t3.0\tanagenetic\t1\tNA\tNA\n" + \
            "1\t8\t4.0\t0.0\t1\t2\t2.5\tanagenetic\t1\tNA\tNA\n" + \
            "1\t8\t4.0\t0.0\t2\t0\t2.0\tanagenetic\t1\tNA\tNA\n" + \
            "1\t8\t4.0\t0.0\t0\t2\t1.0\tanagenetic\t1\tNA\tNA\n" + \
            "1\t3\t4.0\t3.0\t2\t0\t3.0\tcladogenetic\t1\t5\t9\n" + \
            "1\t3\t4.0\t3.0\t2\t1\t3.0\tcladogenetic\t1\t5\t9\n" + \
            "1\t3\t4.0\t3.0\t2\t2\tNA\tno_change\t1\t5\t9\n" + \
            "1\t5\t3.0\t2.0\t2\t1\t2.0\tcladogenetic\t3\t10\t11\n" + \
            "1\t5\t3.0\t2.0\t1\t2\t2.5\tanagenetic\t3\t10\t11\n" + \
            "1\t10\t2.0\t0.0\t1\t1\tNA\tno_change\t5\tNA\tNA\n" + \
            "1\t11\t2.0\t0.0\t2\t0\t1.0\tanagenetic\t5\tNA\tNA\n" + \
            "1\t9\t3.0\t0.0\t0\t0\tNA\tno_change\t3\tNA\tNA\n"

        self.assertEqual(smaps_str_tr2_2, exp_content_tr2_2)


if __name__ == "__main__":
    # $ python3.11 tests/readwrite/test_write_smap.py
    #
    # or
    #
    # $ python3.11 -m tests.readwrite.test_write_smap
    #
    # or
    #
    # $ python3.11 -m unittest tests.readwrite.test_write_smap.TestSmapWrite.test_make_smap_str

    unittest.main()


