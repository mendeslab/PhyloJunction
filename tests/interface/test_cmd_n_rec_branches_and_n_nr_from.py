import unittest

# pj imports
import phylojunction.interface.cmdbox.cmd_parse as cmd
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestNRecBranchesAndNNrFrom(unittest.TestCase):
    """Test the 'n_rec_branches' det-fn and lognormal's 'n_nr_from' argument.

    These let a per-branch relative rate be drawn outside of
    'gamma_spike_clock' (see 'calculation.gamma_spike_clock_model')
    and fed back into it via its own 'br_specific_rel_rate' argument
    -- these tests cover the grammar-level plumbing of that whole
    chain, the same way 'test_cmd_det_fn_clock.py' covers
    'gamma_spike_clock' itself.
    """

    preamble = (
        "birth_rate <- 1.0\n"
        "death_rate <- 0.2\n"
        "det_birth_rate := sse_rate(name=\"birth_rate\", value=birth_rate, "
        "event=\"speciation\")\n"
        "det_death_rate := sse_rate(name=\"death_rate\", value=death_rate, "
        "event=\"extinction\")\n"
        "stash := sse_stash(flat_rate_mat=[det_birth_rate, det_death_rate], "
        "seed_age=3.0)\n"
        "ori <- 3.0\n"
        "trs ~ discrete_sse(n=2, nr=2, stash=stash, start_state=[0,0], "
        "stop=\"age\", stop_value=ori, origin=\"false\", cond_spn=\"true\", "
        "cond_surv=\"true\")\n"
        "n_br := n_rec_branches(tree=trs)\n"
    )

    def test_n_rec_branches_shape_and_values(self):
        """'n_rec_branches' mirrors its tree's shape; values match directly."""

        dag_obj = cmd.script2dag(self.preamble, in_pj_file=False, random_seed=7)
        n_br_node = dag_obj.get_node_dag_by_name("n_br")
        trs_node = dag_obj.get_node_dag_by_name("trs")

        self.assertEqual(n_br_node.sample_size, trs_node.sample_size)
        self.assertEqual(n_br_node.repl_size, trs_node.repl_size)
        self.assertEqual(len(n_br_node.value), len(trs_node.value))

        for i, ann_tr in enumerate(trs_node.value):
            self.assertEqual(ann_tr.n_rec_branches(), n_br_node.value[i])

    def test_n_rec_branches_requires_a_tree_dag_node(self):
        """'n_rec_branches(tree=...)' with a non-DAG-node argument errors."""

        with self.assertRaises(ec.ScriptSyntaxError) as exc:
            cmd.script2dag(
                self.preamble +
                "n_br2 := n_rec_branches(tree=\"not_a_tree\")\n",
                in_pj_file=False)

        self.assertIn(
            "Argument must be a DAG node holding a phylogenetic tree.",
            str(exc.exception))

    def test_n_nr_from_happy_path(self):
        """'n_nr_from' infers shape and draws one vector per tree."""

        script = self.preamble + (
            "br_specific_rel_rate ~ lognormal(n_nr_from=n_br, "
            "meanlog=-0.25, sdlog=0.4, log_space=\"true\")\n")

        dag_obj = cmd.script2dag(script, in_pj_file=False, random_seed=7)
        n_br_node = dag_obj.get_node_dag_by_name("n_br")
        rate_node = dag_obj.get_node_dag_by_name("br_specific_rel_rate")

        self.assertEqual(rate_node.sample_size, n_br_node.sample_size)
        self.assertEqual(rate_node.repl_size, n_br_node.repl_size)
        self.assertEqual(
            [len(v) for v in rate_node.value], n_br_node.value)

    def test_n_nr_from_rejects_literal_value(self):
        """'n_nr_from' must reference a DAG node, not a literal."""

        script = self.preamble + (
            "br_specific_rel_rate ~ lognormal(n_nr_from=5, meanlog=-0.25, "
            "sdlog=0.4)\n")

        with self.assertRaises(ec.ScriptSyntaxError) as exc:
            cmd.script2dag(script, in_pj_file=False)

        self.assertIn(
            "Argument must be a DAG node reference", str(exc.exception))

    def test_n_redirects_to_n_nr_from_when_node_has_real_shape(self):
        """'n=' given a node with real sample/replicate structure errors.

        Rather than the generic "too many values" error, the user
        should be pointed at 'n_nr_from' instead -- see 'nr=' below
        for the symmetric case.
        """

        script = self.preamble + (
            "br_specific_rel_rate ~ lognormal(n=n_br, nr=1, meanlog=-0.25, "
            "sdlog=0.4)\n")

        with self.assertRaises(ec.ScriptSyntaxError) as exc:
            cmd.script2dag(script, in_pj_file=False)

        self.assertIn("use 'n_nr_from' instead of 'n'", str(exc.exception))

    def test_nr_redirects_to_n_nr_from_when_node_has_real_shape(self):
        """Symmetric check: 'nr=' given such a node errors the same way."""

        script = self.preamble + (
            "br_specific_rel_rate ~ lognormal(n=1, nr=n_br, meanlog=-0.25, "
            "sdlog=0.4)\n")

        with self.assertRaises(ec.ScriptSyntaxError) as exc:
            cmd.script2dag(script, in_pj_file=False)

        self.assertIn("use 'n_nr_from' instead of 'nr'", str(exc.exception))

    def test_n_nr_from_mutually_exclusive_with_n(self):
        """'n_nr_from' combined with 'n' (even a harmless literal) errors."""

        script = self.preamble + (
            "br_specific_rel_rate ~ lognormal(n_nr_from=n_br, n=4, "
            "meanlog=-0.25, sdlog=0.4)\n")

        with self.assertRaises(ec.ScriptSyntaxError) as exc:
            cmd.script2dag(script, in_pj_file=False)

        self.assertIn("Cannot be combined with 'n_nr_from'", str(exc.exception))

    def test_gamma_spike_clock_br_specific_rel_rate_end_to_end(self):
        """'br_specific_rel_rate' values end up exactly as each branch's raw rate."""

        script = self.preamble + (
            "br_specific_rel_rate ~ lognormal(n_nr_from=n_br, "
            "meanlog=-0.25, sdlog=0.4, log_space=\"true\")\n"
            "phylo := gamma_spike_clock(tree=trs, "
            "attr_name=\"mutation_rate\", clock_rate=0.01, spike_mean=0.1, "
            "spike_shape=2.0, "
            "br_specific_rel_rate=br_specific_rel_rate, "
            "strict_clock_when_no_stubs=\"false\")\n")

        dag_obj = cmd.script2dag(script, in_pj_file=False, random_seed=7)
        rate_node = dag_obj.get_node_dag_by_name("br_specific_rel_rate")
        phylo_node = dag_obj.get_node_dag_by_name("phylo")

        for i, ann_tr in enumerate(phylo_node.value):
            expected = list(rate_node.value[i])
            idx = 0

            for nd in ann_tr.tree_reconstructed.preorder_node_iter():
                if nd.branch_raw_rate is not None:
                    self.assertEqual(nd.branch_raw_rate, expected[idx])
                    idx += 1

            self.assertEqual(idx, len(expected))


if __name__ == "__main__":
    unittest.main()
