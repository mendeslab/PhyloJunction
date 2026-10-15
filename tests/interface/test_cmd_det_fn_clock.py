import unittest

# pj imports
import phylojunction.interface.cmdbox.cmd_parse as cmd
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestDetFnClockErrors(unittest.TestCase):
    """Test 'gamma_spike_clock' through actual .pj script parsing.

    Unlike 'tests/calculation/test_gamma_spike_clock_model.py' (which
    calls 'GammaSpikeClockModel'/its free function directly), these
    tests go through the real ':=' grammar pathway
    ('interface.grammar.det_fn_clock_makers.make_gamma_spike_clock_tree'
    and its wrapping in 'interface.grammar.det_fn_grammar' and
    'interface.cmdbox.cmd_parse'), the same way
    'test_cmd_det_fn_discrete_sse.py' does for 'sse_rate'/'sse_stash'.
    """

    # shared preamble: a minimal birth-death chronogram, plus the
    # 'br_specific_rel_rate' pipeline 'gamma_spike_clock' now always
    # requires (see 'test_cmd_n_rec_branches_and_n_nr_from.py' for
    # dedicated coverage of that pipeline itself)
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
        "trs ~ discrete_sse(n=1, nr=1, stash=stash, start_state=[0], "
        "stop=\"age\", stop_value=ori, origin=\"false\", cond_spn=\"true\", "
        "cond_surv=\"true\")\n"
        "n_br := n_rec_branches(tree=trs)\n"
        "br_specific_rel_rate ~ lognormal(n_nr_from=n_br, meanlog=-0.125, "
        "sdlog=0.4, log_space=\"true\")\n"
    )

    def test_missing_attr_name_error(self):
        """Test that omitting a required parameter ('attr_name') errors.

        Only one required parameter is tested missing here (rather
        than one test per parameter) -- 'attr_name' is picked because
        it is specific to this det-fn (unlike the shared float
        hyperparameters), so this also doubles as a check that
        'attr_name' really is treated as required.
        """

        missing_attr_name_line = (
            "phylo := gamma_spike_clock(tree=trs, clock_rate=0.01, "
            "br_specific_rel_rate=br_specific_rel_rate, spike_mean=0.1, "
            "spike_shape=2.0)")

        with self.assertRaises(ec.ScriptSyntaxError) as exc:
            cmd.script2dag(
                self.preamble + missing_attr_name_line, in_pj_file=False)

        self.assertEqual(
            str(exc.exception),
            ("ERROR: " + missing_attr_name_line + "\n\n"
             "The line above had a syntax problem and could not be "
             "tokenized. Deterministic output from 'gamma_spike_clock' "
             "could not be instantiated. Parameter 'attr_name' is "
             "missing."))

    def test_missing_br_specific_rel_rate_error(self):
        """Test that omitting 'br_specific_rel_rate' errors.

        This model never draws its gradual, relaxed-clock component
        internally -- 'br_specific_rel_rate' must always be given.
        """

        missing_rate_line = (
            "phylo := gamma_spike_clock(tree=trs, "
            "attr_name=\"mutation_rate\", clock_rate=0.01, "
            "spike_mean=0.1, spike_shape=2.0)")

        with self.assertRaises(ec.ScriptSyntaxError) as exc:
            cmd.script2dag(self.preamble + missing_rate_line, in_pj_file=False)

        self.assertIn(
            "Parameter 'br_specific_rel_rate' is missing.",
            str(exc.exception))

    def test_invalid_strict_clock_when_no_stubs_error(self):
        """Test that an invalid 'strict_clock_when_no_stubs' value errors.

        Only one parameter is tested with an invalid argument here
        (rather than one test per parameter) -- this one is picked
        because its "true"/"false" string validation is new code
        added this session, with no prior coverage at all.
        """

        invalid_line = (
            "phylo := gamma_spike_clock(tree=trs, "
            "attr_name=\"mutation_rate\", clock_rate=0.01, "
            "br_specific_rel_rate=br_specific_rel_rate, "
            "spike_mean=0.1, spike_shape=2.0, "
            "strict_clock_when_no_stubs=\"maybe\")")

        with self.assertRaises(ec.ScriptSyntaxError) as exc:
            cmd.script2dag(self.preamble + invalid_line, in_pj_file=False)

        self.assertEqual(
            str(exc.exception),
            ("ERROR: " + invalid_line + "\n\n"
             "The line above had a syntax problem and could not be "
             "tokenized. Deterministic output from 'gamma_spike_clock' "
             "could not be instantiated. Argument '\"maybe\"' is invalid "
             "for parameter 'strict_clock_when_no_stubs'. Argument must "
             "be \"true\" or \"false\"."))

    def test_strict_clock_when_no_stubs_defaults_to_false(self):
        """Test 'strict_clock_when_no_stubs' defaults to 'False' if omitted.

        There is no way to inspect the ephemeral 'GammaSpikeClockModel'
        instance 'make_gamma_spike_clock_tree' builds internally (only
        its already-'.generate()'-d trees survive, as the det-fn
        node's value), so this is a behavioral check instead: under
        the real model ('strict_clock_when_no_stubs=False'), EVERY
        branch (stub-free or not) gets its own independent log-normal
        'branch_raw_rate' draw, so a stub-free branch's draw is (with
        probability 1, since the log-normal is a continuous
        distribution) never exactly 1.0. Under
        'strict_clock_when_no_stubs=True', by contrast, a stub-free
        branch's 'branch_raw_rate' is *always* pinned to exactly 1.0
        (see 'populate_gamma_spike_mol_rates_on_rec_tr_branches'). So
        finding a stub-free branch whose 'branch_raw_rate' differs
        from 1.0 rules out the argument having silently defaulted to
        'True'.
        """

        # 'strict_clock_when_no_stubs' is deliberately left out
        no_flag_line = (
            "phylo := gamma_spike_clock(tree=trs, "
            "attr_name=\"mutation_rate\", clock_rate=0.01, "
            "br_specific_rel_rate=br_specific_rel_rate, "
            "spike_mean=0.1, spike_shape=2.0)")

        dag_obj = cmd.script2dag(
            self.preamble + no_flag_line, in_pj_file=False, random_seed=42)

        ann_tr = dag_obj.get_node_dag_by_name("phylo").value[0]

        checked_a_stub_free_branch = False
        for rec_nd in ann_tr.tree_reconstructed.preorder_node_iter():
            if rec_nd.n_stub == 0 and rec_nd.branch_raw_rate is not None:
                self.assertNotEqual(rec_nd.branch_raw_rate, 1.0)
                checked_a_stub_free_branch = True
                break

        # sanity check on the fixture/seed itself: if this ever stops
        # finding a stub-free branch, the assertion above is never
        # actually exercised
        self.assertTrue(checked_a_stub_free_branch)


if __name__ == "__main__":
    unittest.main()
