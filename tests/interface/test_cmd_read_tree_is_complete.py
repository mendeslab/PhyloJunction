import unittest

# pj imports
import phylojunction.interface.cmdbox.cmd_parse as cmd
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestReadTreeIsComplete(unittest.TestCase):
    """Test 'read_tree's 'is_complete'/'keep_unsampled_alive' grammar.

    Complements 'tests.readwrite.test_pj_read_annotated_newick',
    which exercises 'readwrite.pj_read.read_nwk_tree_str' directly --
    these tests instead go through the full '.pj' grammar/DAG-
    building pipeline ('cmd.script2dag'), confirming the two new
    arguments are correctly parsed and threaded all the way down to
    it.
    """

    def test_is_complete_true_with_fossil_tip(self):
        """'is_complete="true"' is parsed and honored end to end."""

        script = (
            'trs <- read_tree(n=1, string="(tipa[&alive=True]:2.0,'
            'fossiltip[&alive=False,sampled=True]:1.0):0.0;", '
            'is_complete="true")\n')

        dag_obj = cmd.script2dag(script, in_pj_file=False)
        trs_node = dag_obj.get_node_dag_by_name("trs")
        ann_tr = trs_node.value[0]

        fossiltip = ann_tr.tree.find_node_with_taxon_label("fossiltip")
        self.assertFalse(fossiltip.alive)
        self.assertTrue(fossiltip.sampled)

        tr_rec = ann_tr.extract_reconstructed_tree(
            plotting_overhead=False, require_obs_both_sides=False)
        rec_leaf_labels = {
            nd.taxon.label for nd in tr_rec.leaf_node_iter()}
        self.assertIn("fossiltip", rec_leaf_labels)

    def test_is_complete_defaults_to_false(self):
        """Omitting 'is_complete' keeps the safe, annotation-free default."""

        script = (
            'trs <- read_tree(n=1, string="(nd1:2.0,nd2:2.0):0.0;")\n')

        # must not raise -- no annotations needed
        dag_obj = cmd.script2dag(script, in_pj_file=False)
        trs_node = dag_obj.get_node_dag_by_name("trs")
        ann_tr = trs_node.value[0]

        tipa = ann_tr.tree.find_node_with_taxon_label("nd1")
        self.assertTrue(tipa.sampled)

    def test_is_complete_true_missing_annotation_raises(self):
        """A missing required annotation surfaces through the DAG build.

        'AnnotatedTreeMissingNodeAnnotationError' is raised from deep
        inside 'read_nwk_tree_str', while the tree's value is actually
        generated -- confirm it is not swallowed or rewrapped into a
        generic parsing error anywhere between there and 'script2dag'.
        """

        script = (
            'trs <- read_tree(n=1, string="(nd1[&alive=True]:2.0,'
            'deadtip:1.0):0.0;", is_complete="true")\n')

        with self.assertRaises(
                ec.AnnotatedTreeMissingNodeAnnotationError) as exc:
            cmd.script2dag(script, in_pj_file=False)

        self.assertIn("deadtip", str(exc.exception))

    def test_is_complete_invalid_value_raises(self):
        """'is_complete' must be a quoted "true"/"false" string."""

        script = (
            'trs <- read_tree(n=1, string="(nd1:2.0,nd2:2.0):0.0;", '
            'is_complete="yes")\n')

        with self.assertRaises(ec.ParseInvalidArgumentError) as exc:
            cmd.script2dag(script, in_pj_file=False)

        self.assertIn("is_complete", str(exc.exception))


if __name__ == "__main__":
    unittest.main()
