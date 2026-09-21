import unittest
from phylojunction.data.trait import DiscreteTrait
from phylojunction.pgm.pgm import DeterministicNodeDAG
from phylojunction.interface.grammar.det_fn_map_attribute import make_mapped_ann_tree

# pj imports
import phylojunction.readwrite.pj_read as pjr
import phylojunction.functionality.biogeo as pjbio
import phylojunction.functionality.stoch_map as pjsmap

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class TestStochMaps(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:

        tree_file_path = \
            "examples/trees_maps_files/geosse_dummy_tree2.tre"
        n_states = 3

        ann_tr_list = [pjr.read_nwk_tree_str(tree_file_path,
                                             "read_tree",
                                             node_names_attribute="index",
                                             n_states=n_states,
                                             in_file=True)]

        n_chars = 2
        state2bit_lookup = pjbio.State2BitLookup(n_chars,
                                                 2,
                                                 geosse=True)

        maps_file_path = \
            "examples/trees_maps_files/geosse_dummy_tree2_maps.tsv"
        node_states_file_path = \
            "examples/trees_maps_files/geosse_dummy_tree2_tip_states.tsv"

        cls.stoch_mapcoll = \
            pjsmap.StochMapsOnTreeCollection(maps_file_path,
                                             ann_tr_list,
                                             state2bit_lookup,
                                             node_states_file_path=node_states_file_path,
                                             stoch_map_attr_name="state")

    # Mapping must preserve custom storage names when replacing the state count.
    # Default-name tests miss this; retire if a shared mapping test covers renamed attributes.
    def test_named_trait_mapping(self):
        tree = pjr.read_nwk_tree_str("examples/trees_maps_files/geosse_dummy_tree2.tre",
                                    "read_tree", node_names_attribute="index",
                                    n_states=1, in_file=True)
        tree.trait = DiscreteTrait(1, name="habitat")
        for nd in tree.tree:
            nd.habitat = nd.state
            del nd.state
        mapped = make_mapped_ann_tree("map_attr", {
            "tree": [DeterministicNodeDAG("trees", value=[tree])],
            "fun": ['"smap"'], "n_regions": ["2"], "geo": ['"true"'],
            "maps_file_path": ['"examples/trees_maps_files/geosse_dummy_tree2_maps.tsv"'],
            "tip_attr_file_path": ['"examples/trees_maps_files/geosse_dummy_tree2_tip_states.tsv"'],
        })
        self.assertEqual(tree.trait, DiscreteTrait(1, name="habitat"))
        for result in mapped:
            self.assertEqual(result.trait, DiscreteTrait(3, name="habitat"))
            self.assertTrue(all(hasattr(nd, "habitat") and not hasattr(nd, "state")
                                for nd in result.tree))
            transitions = [event for events in result.at_dict.values() for event in events]
            transitions += list(result.clado_at_dict.values())
            self.assertTrue(transitions)
            self.assertTrue(all(event.attr_label == "habitat" for event in transitions))

    def test_read_stoch_maps(self) -> None:
        """
        Test that manually constructed stochastic maps (for iteration
        1) on five-taxon tree are read in correctly.
        """

        expected_str_representation = \
            ("Stochastic maps in MCMC iteration 1:\n"
             "    Number of anagenetic changes = 3 (50.0%)\n"
             "        Range contractions = 0\n"
             "        Range expansions = 3\n"
             "    Number of higher-order anagenetic changes = 0 (0.0%)\n"
             "    Number of cladogenetic changes = 3 (50.0%)\n"
             "    Number of identical cladogenetic changes = 0 (0.0%)")

        self.assertEqual(
            str(self.stoch_mapcoll.stoch_maps_tree_dict[1]),
            expected_str_representation)


if __name__ == '__main__':
    # From PhyloJunction/
    #
    # $ python3 tests/functionality/test_stoch_maps.py
    #
    # or
    #
    # $ python3 -m tests.functionality.test_stoch_maps
    #
    # or
    #
    # $ python3 -m unittest tests.functionality.test_stoch_maps.TestStochMaps.test_read_stoch_maps

    unittest.main()
