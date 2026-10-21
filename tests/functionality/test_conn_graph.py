import unittest
from itertools import permutations

# pj imports
import phylojunction.functionality.feature_io as pjgeo

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestConnGraph(unittest.TestCase):

    # Protect symmetric adjacency and existing neighbors under repeated insertion;
    # directed-only tests miss this. Remove if adjacency storage is replaced.
    def test_undirected_adjacency(self):
        g = pjgeo.GeoGraph(3)
        for a, b in [(0, 1), (0, 2), (0, 1), (2, 0)]:
            g.add_edge(a, b)
        self.assertEqual(g.edge_dict, {0: {1, 2}, 1: {0}, 2: {0}})
        self.assertEqual(g.edge_set, {(0, 1), (1, 0), (0, 2), (2, 0)})
        g.populate_comm_class_members()
        self.assertTrue(g.are_connected(1, 2))

    # Connected components must partition nodes independently of insertion order and
    # rebuild cleanly; fixed-order examples miss overlapping merges. Retire if graphs are replaced.
    def test_component_partition(self):
        edges = [(0, 1), (2, 3), (4, 5), (5, 3)]
        for order in permutations(edges):
            with self.subTest(order=order):
                g = pjgeo.GeoGraph(7)
                for source, target in order:
                    g.add_edge(source, target, is_directed=True)
                for _ in range(2):
                    g.populate_comm_class_members()
                    self.assertEqual(g.comm_class_set_list, [{0, 1}, {2, 3, 4, 5}, {6}])
                    self.assertEqual(g.n_comm_classes, 3)
                    self.assertEqual(sum(map(len, g.comm_class_set_list)), 7)
                    self.assertEqual(set.union(*g.comm_class_set_list), set(range(7)))
                    for source, target in edges:
                        self.assertTrue(g.are_connected(source, target))
                    self.assertTrue(g.are_connected(2, 4))
                    self.assertFalse(g.are_connected(0, 6))

    def test_graph1(self) -> None:
        """Test multiple comm. classes are correctly built."""

        n_regions = 7

        g = pjgeo.GeoGraph(n_regions)

        # comm classes: {0}, {1}, {2}, {3}, {4}, {5}, {6}
        # (no edges!)

        g.populate_comm_class_members()

        # debugging
        # print(g.comm_class_set_list)

        self.assertEqual(g.comm_class_set_list,
                         [{0}, {1}, {2}, {3}, {4}, {5}, {6}])
        self.assertFalse(g.are_connected(0, 1))
        self.assertFalse(g.are_connected(2, 3))
        self.assertFalse(g.are_connected(4, 5))
        self.assertFalse(g.are_connected(0, 2))
        self.assertFalse(g.are_connected(0, 4))
        self.assertFalse(g.are_connected(2, 4))
        self.assertFalse(g.are_connected(6, 0))
        self.assertFalse(g.are_connected(6, 3))
        self.assertFalse(g.are_connected(6, 4))

    def test_graph2(self) -> None:
        """Test multiple comm. classes are correctly built."""

        n_regions = 7

        g = pjgeo.GeoGraph(n_regions)

        # comm classes: {0,1}, {2,3}, {4,5}, {6}
        g.add_edge(0, 1, is_directed=True)
        g.add_edge(2, 3, is_directed=True)
        g.add_edge(4, 5, is_directed=True)

        g.populate_comm_class_members()

        # debugging
        # print(g.comm_class_set_list)

        self.assertEqual(g.comm_class_set_list,
                         [{0, 1}, {2, 3}, {4, 5}, {6}])
        self.assertTrue(g.are_connected(0, 1))
        self.assertTrue(g.are_connected(2, 3))
        self.assertTrue(g.are_connected(4, 5))
        self.assertFalse(g.are_connected(0, 2))
        self.assertFalse(g.are_connected(0, 4))
        self.assertFalse(g.are_connected(2, 4))
        self.assertFalse(g.are_connected(6, 0))
        self.assertFalse(g.are_connected(6, 3))
        self.assertFalse(g.are_connected(6, 4))

    def test_graph3(self) -> None:
        """Test multiple single comm. classe is correctly built."""

        n_regions = 7

        g = pjgeo.GeoGraph(n_regions)

        # comm classes: {0,1,2,3,4,5,6}
        g.add_edge(0, 1, is_directed=True)
        g.add_edge(2, 3, is_directed=True)
        g.add_edge(4, 5, is_directed=True)
        g.add_edge(6, 0, is_directed=True)
        g.add_edge(6, 2, is_directed=True)
        g.add_edge(6, 4, is_directed=True)

        g.populate_comm_class_members()

        # debugging
        # print(g.comm_class_set_list)

        self.assertEqual(g.comm_class_set_list,
                         [{0, 1, 2, 3, 4, 5, 6}])
        self.assertTrue(g.are_connected(0, 1))
        self.assertTrue(g.are_connected(2, 3))
        self.assertTrue(g.are_connected(4, 5))
        self.assertTrue(g.are_connected(0, 2))
        self.assertTrue(g.are_connected(0, 4))
        self.assertTrue(g.are_connected(2, 4))
        self.assertTrue(g.are_connected(6, 0))
        self.assertTrue(g.are_connected(6, 3))
        self.assertTrue(g.are_connected(6, 4))


if __name__ == '__main__':
    # From PhyloJunction/
    #
    # $ python3 tests/functionality/test_conn_graph.py
    #
    # or
    #
    # $ python3 -m tests.functionality.test_conn_graph
    #
    # or
    #
    # $ python3 -m unittest tests.functionality.test_conn_graph.TestConnGraph.test_graph1

    unittest.main()