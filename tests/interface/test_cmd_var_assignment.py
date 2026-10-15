import unittest
import re
import tempfile
import os
from pathlib import Path

EXAMPLES_DIR = Path(__file__).parent.parent.parent / "examples" / "trees_maps_files"
FILE_PATH1_STR = str(EXAMPLES_DIR / "tree_to_read.tre")
FILE_PATH2_STR = str(EXAMPLES_DIR / "trees_to_read.tre")

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.interface.cmdbox.cmd_parse as cmd
import phylojunction.interface.cmdbox.cmd_parse_utils as cmdu
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestVarAssignment(unittest.TestCase):

    def test_var_assignment(self):
        """
        Test variable assignment.
        
        Test a series of different variable assignments produce the
        right DAG.
        """

        dag_obj = pgm.DirectedAcyclicGraph()
        
        cmd_line1 = "a <- 1"

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.assign_regex, cmd_line1)
        cmd.parse_variable_assignment(dag_obj,
                                      stoch_node_name,
                                      stoch_node_spec,
                                      cmd_line1)
        a_node_dag = dag_obj.get_node_dag_by_name("a")
        
        self.assertEqual(1, dag_obj.n_nodes)
        self.assertEqual(type(a_node_dag.value), list)
        self.assertEqual(len(a_node_dag.value), 1)
        self.assertEqual(a_node_dag.value, ["1"])
        
        # --- #

        cmd_line2 = "b <- [1, 2, 3]"

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.assign_regex, cmd_line2)
        cmd.parse_variable_assignment(dag_obj,
                                      stoch_node_name,
                                      stoch_node_spec,
                                      cmd_line2)
        b_node_dag = dag_obj.get_node_dag_by_name("b")

        self.assertEqual(2, dag_obj.n_nodes)
        self.assertEqual(type(b_node_dag.value), list)
        self.assertEqual(len(b_node_dag.value), 3)
        self.assertEqual(b_node_dag.value, ["1", "2", "3"])
        
        # --- #
        
        cmd_line3 = "a <- [1]"

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.assign_regex, cmd_line3)
        cmd.parse_variable_assignment(dag_obj,
                                      stoch_node_name,
                                      stoch_node_spec,
                                      cmd_line3)
        a_node_dag = dag_obj.get_node_dag_by_name("a")

        self.assertEqual(2, dag_obj.n_nodes)
        self.assertEqual(type(a_node_dag.value), list)
        self.assertEqual(len(a_node_dag.value), 1)
        self.assertEqual(a_node_dag.value, ["1"])

        # --- #
        
        cmd_line4 = "c <- b"

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.assign_regex, cmd_line4)
        cmd.parse_variable_assignment(dag_obj,
                                      stoch_node_name,
                                      stoch_node_spec,
                                      cmd_line4)
        c_node_dag = dag_obj.get_node_dag_by_name("c")

        self.assertEqual(3, dag_obj.n_nodes)
        self.assertEqual(type(c_node_dag.value), list)
        self.assertEqual(len(c_node_dag.value), 3)
        self.assertEqual(c_node_dag.value, ["1", "2", "3"])

        # --- # 

        cmd_line5 = "d <- [c, 4, 5, 6]"

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.assign_regex, cmd_line5)
        cmd.parse_variable_assignment(dag_obj,
                                      stoch_node_name,
                                      stoch_node_spec,
                                      cmd_line5)
        d_node_dag = dag_obj.get_node_dag_by_name("d")

        self.assertEqual(4, dag_obj.n_nodes)
        self.assertEqual(type(d_node_dag.value), list)
        self.assertEqual(len(d_node_dag.value), 6)
        self.assertEqual(d_node_dag.value, ["1", "2", "3", "4", "5", "6"])

    def test_var_assignment_read_tree_string(self):
        """Test read_tree() from string.
        
        Test read_tree() calls using Newick strings directly produce
        the right DAG.
        """

        dag_obj = pgm.DirectedAcyclicGraph()

        cmd_line1 = ('tr <- read_tree(string="((sp1[&index=1]:1.0,sp2'
                     '[&index=2]:1.0)[&index=4]:1.0,sp3[&index=3]:2.0)'
                     '[&index=5];", node_name_attr="index")')

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.assign_regex, cmd_line1)
        cmd.parse_variable_assignment(dag_obj,
                                      stoch_node_name,
                                      stoch_node_spec,
                                      cmd_line1)
        a_node_dag = dag_obj.get_node_dag_by_name("tr")

        self.assertEqual(1, dag_obj.n_nodes)
        self.assertEqual(type(a_node_dag.value), list)
        self.assertEqual(len(a_node_dag.value), 1)
        self.assertEqual(a_node_dag.value[0].__str__().rstrip(),
                         ('((nd1:1.0[&index=1],nd2:1.0[&index=2])nd4:1.0'
                          '[&index=4],nd3:2.0[&index=3])nd5[&index=5];'))
        
        cmd_line2 = 'tr <- read_tree(string="((sp1:1.0,sp2:1.0):1.0,sp3:2.0);")'

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.assign_regex, cmd_line2)
        cmd.parse_variable_assignment(dag_obj, stoch_node_name, stoch_node_spec, cmd_line2)
        a_node_dag = dag_obj.get_node_dag_by_name("tr")

        self.assertEqual(1, dag_obj.n_nodes)
        self.assertEqual(type(a_node_dag.value), list)
        self.assertEqual(len(a_node_dag.value), 1)
        self.assertEqual(a_node_dag.value[0].__str__().rstrip(),
                         ('((sp1:1.0,sp2:1.0)nd1:1.0,sp3:2.0)root;'))
        
    def test_var_assignment_read_tree_sample_size(self):
        """'read_tree()' (via '<-') correctly sets 'sample_size'.

        'read_tree' is a constant function, not a sampling
        distribution -- its value only comes into existence lazily,
        inside 'StochasticNodeDAG.generate_value()', well after
        'parse_variable_assignment' has to decide this node's
        'sample_size'. That decision used to always come out as 0,
        regardless of how many trees 'read_tree' actually produced,
        because the code path only populated an (always-empty, for
        this kind of call) 'values_list' and read 'sample_size' off of
        that -- rather than off of the constant function's own
        'n_samples' (an attribute every concrete 'ConstantFn' is
        required to expose). A 'sample_size' of 0 does not matter for
        'read_tree' on its own (nothing else inspects it), but it does
        matter the moment this tree is then fed into something that
        checks shapes against another DAG node's own correct
        'sample_size' -- e.g. 'gamma_spike_clock's 'br_specific_rel_rate'
        argument, built from a 'lognormal(n_nr_from=...)' node whose
        own 'sample_size' is correctly inferred.
        """

        dag_obj = pgm.DirectedAcyclicGraph()

        cmd_line = ('tr <- read_tree(string="((sp1:1.0,sp2:1.0):1.0,'
                    'sp3:2.0);")')

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.assign_regex, cmd_line)
        cmd.parse_variable_assignment(
            dag_obj, stoch_node_name, stoch_node_spec, cmd_line)
        a_node_dag = dag_obj.get_node_dag_by_name("tr")

        self.assertEqual(a_node_dag.sample_size, 1)
        self.assertEqual(a_node_dag.repl_size, 1)
        self.assertEqual(len(a_node_dag.value), 1)

    def test_var_assignment_read_tree_from_file(self):
        """Test read_tree() from file.
        
        Test read_tree()calls using tree files containing Newick
        strings directly produce the right DAG.
        """

        dag_obj = pgm.DirectedAcyclicGraph()

        cmd_line1 = ('tr <- read_tree(file_path=\"' +
                     FILE_PATH1_STR + '\",'
                     'node_name_attr="index")')

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.assign_regex, cmd_line1)
        cmd.parse_variable_assignment(dag_obj,
                                      stoch_node_name,
                                      stoch_node_spec,
                                      cmd_line1)
        a_node_dag = dag_obj.get_node_dag_by_name("tr")

        self.assertEqual(1, dag_obj.n_nodes)
        self.assertEqual(type(a_node_dag.value), list)
        self.assertEqual(len(a_node_dag.value), 1)
        self.assertEqual(a_node_dag.value[0].__str__().rstrip(),
                         ('((nd1:1.0[&index=1],nd2:1.0[&index=2])nd4:1.0'
                          '[&index=4],nd3:2.0[&index=3])nd5[&index=5];'))
        
        cmd_line2 = ('tr <- read_tree(file_path=\"' +
                     FILE_PATH2_STR + '\",'
                     'node_name_attr="index")')

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.assign_regex, cmd_line2)
        cmd.parse_variable_assignment(dag_obj, stoch_node_name, stoch_node_spec, cmd_line2)
        a_node_dag = dag_obj.get_node_dag_by_name("tr")

        self.assertEqual(1, dag_obj.n_nodes)
        self.assertEqual(type(a_node_dag.value), list)
        self.assertEqual(len(a_node_dag.value), 2)
        self.assertEqual(a_node_dag.value[1].__str__().rstrip(),
                         ('((nd1:0.1[&index=1],nd2:0.1[&index=2])nd4:1.9'
                          '[&index=4],nd3:2.0[&index=3])nd5[&index=5];'))

    def test_var_assignment_read_tree_string_exceptions(self):
        """Test read_tree() from string raised exceptions."""

        dag_obj = pgm.DirectedAcyclicGraph()

        cmd_line1 = ('tr <- read_tree(n=2, string="((sp1[&index=1]:1.0,'
                     'sp2[&index=2]:1.0)[&index=4]:1.0,sp3[&index=3]:2.0)'
                     '[&index=5];", node_name_attr="index")')

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.assign_regex, cmd_line1)
        
        with self.assertRaises(ec.ObjInitIncorrectDimensionError) as exc:
            cmd.parse_variable_assignment(dag_obj,
                                          stoch_node_name,
                                          stoch_node_spec,
                                          cmd_line1)

        expected_exception_message1 = \
            ("Could not initialize the object of CtFnTreeReader. "
             "Incorrect dimension of container string, which was of "
             "size 1. The expected dimension was 2.")

        self.assertEqual(str(exc.exception),
                         expected_exception_message1)

        cmd_line2 = ('tr <- read_tree(file_path="examples/trees_maps_files/'
                     'tree_to_read.tre", string="((sp1[&index=1]:1.0,'
                     'sp2[&index=2]:1.0)[&index=4]:1.0,sp3[&index=3]:2.0)'
                     '[&index=5];", node_name_attr="index")')

        stoch_node_name, _, stoch_node_spec = \
            re.split(cmdu.assign_regex, cmd_line2)

        expected_exception_message2 = \
            ("Argument file_path cannot be specified at the same time as 'string'.")

        with self.assertRaises(ec.ParseMutuallyExclusiveParametersError) as exc:
            cmd.parse_variable_assignment(dag_obj,
                                          stoch_node_name,
                                          stoch_node_spec,
                                          cmd_line2)

        self.assertEqual(str(exc.exception).rstrip(),
                         expected_exception_message2)
        
    def test_var_assignment_read_tree_file_exceptions(self):
        """Test read_tree() from file raised exceptions."""

        dag_obj = pgm.DirectedAcyclicGraph()

        cmd_line1 = ('tr <- read_tree(nr=13, file_path=\"' +
                     FILE_PATH1_STR + '\", '
                     'node_name_attr="index")')

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.assign_regex, cmd_line1)

        expected_exception_message1 = \
            ("Could not initialize the object of CtFnTreeReader. "
             "Incorrect dimension of container file_path, which was of "
             "size 1. The expected dimension was 13.")

        with self.assertRaises(ec.ObjInitIncorrectDimensionError) as exc:
            cmd.parse_variable_assignment(dag_obj,
                                          stoch_node_name,
                                          stoch_node_spec,
                                          cmd_line1)

        self.assertEqual(str(exc.exception),
                         expected_exception_message1)
        
        cmd_line2 = ('tr <- read_tree(file_path=\"' +
                     FILE_PATH2_STR + '\", '
                     'node_name_attr="name")')

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.assign_regex, cmd_line2)

        expected_exception_message2 = \
            ("Could not initialize the object of \'read_tree\'. "
             "The argument provided to \'node_name_attribute\' was invalid. "
             "The attribute for node names was not in the newick strings.")
        
        with self.assertRaises(ec.ObjInitInvalidArgError) as exc:
            cmd.parse_variable_assignment(dag_obj,
                                          stoch_node_name,
                                          stoch_node_spec,
                                          cmd_line2)

        self.assertEqual(str(exc.exception),
                         expected_exception_message2)

    def test_var_assignment_read_tree_neither_n_nor_nr_infers_from_file(self):
        """Neither 'n' nor 'nr' given: every tree in the file goes into 'n'.

        A file holding more than one tree string (one per line) is a
        valid, intentional way to provide several samples at once --
        as long as neither 'n' nor 'nr' is given at all, in which case
        every tree string found goes straight into 'n_samples' (never
        'n_repl').
        """

        dag_obj = pgm.DirectedAcyclicGraph()

        with tempfile.NamedTemporaryFile(
                mode="w", suffix=".tre", delete=False) as tmp_file:
            tmp_file.write("(tipa:1.0,tipb:1.0):0.0;\n")
            tmp_file.write("(tipc:2.0,tipd:2.0):0.0;\n")
            tmp_file_path = tmp_file.name

        try:
            cmd_line = 'tr <- read_tree(file_path="' + tmp_file_path + '")'

            stoch_node_name, _, stoch_node_spec = \
                re.split(cmdu.assign_regex, cmd_line)
            cmd.parse_variable_assignment(
                dag_obj, stoch_node_name, stoch_node_spec, cmd_line)
            a_node_dag = dag_obj.get_node_dag_by_name("tr")

            self.assertEqual(a_node_dag.sample_size, 2)
            self.assertEqual(a_node_dag.repl_size, 1)
            self.assertEqual(len(a_node_dag.value), 2)

        finally:
            os.remove(tmp_file_path)

    def test_var_assignment_read_tree_explicit_n_one_rejects_multi_tree_file(self):
        """Explicit 'n=1' is taken literally, even with a multi-tree file.

        Unlike when 'n'/'nr' are both omitted (see
        'test_var_assignment_read_tree_neither_n_nor_nr_infers_from_file'
        above), an explicitly-given 'n=1' is never silently
        reinterpreted just because the file happens to hold more than
        one tree string -- it is checked against the file's actual
        tree count, and an exception is raised on a mismatch, exactly
        as it already was for any other explicit 'n'/'nr' value (see
        'test_var_assignment_read_tree_string_exceptions'/
        'test_var_assignment_read_tree_file_exceptions' above).
        """

        dag_obj = pgm.DirectedAcyclicGraph()

        with tempfile.NamedTemporaryFile(
                mode="w", suffix=".tre", delete=False) as tmp_file:
            tmp_file.write("(tipa:1.0,tipb:1.0):0.0;\n")
            tmp_file.write("(tipc:2.0,tipd:2.0):0.0;\n")
            tmp_file_path = tmp_file.name

        try:
            cmd_line = \
                'tr <- read_tree(n=1, file_path="' + tmp_file_path + '")'

            stoch_node_name, _, stoch_node_spec = \
                re.split(cmdu.assign_regex, cmd_line)

            with self.assertRaises(
                    ec.ObjInitIncorrectDimensionError) as exc:
                cmd.parse_variable_assignment(
                    dag_obj, stoch_node_name, stoch_node_spec, cmd_line)

            expected_exception_message = \
                ("Could not initialize the object of CtFnTreeReader. "
                 "Incorrect dimension of container file_path, which was "
                 "of size 2. The expected dimension was 1.")

            self.assertEqual(str(exc.exception), expected_exception_message)

        finally:
            os.remove(tmp_file_path)

    def test_var_assignment_read_tree_explicit_nr_matches_multi_tree_file(self):
        """Explicit 'nr' (not 'n') is also taken literally and checked.

        'nr=2' alone (no 'n' given) against a 2-tree file satisfies
        'n_samples * n_repl == total_tr_count' (1 * 2 == 2), so this
        does not raise -- either of 'n'/'nr' being explicitly given is
        enough to opt out of the "neither given" auto-inference case.
        """

        dag_obj = pgm.DirectedAcyclicGraph()

        with tempfile.NamedTemporaryFile(
                mode="w", suffix=".tre", delete=False) as tmp_file:
            tmp_file.write("(tipa:1.0,tipb:1.0):0.0;\n")
            tmp_file.write("(tipc:2.0,tipd:2.0):0.0;\n")
            tmp_file_path = tmp_file.name

        try:
            cmd_line = \
                'tr <- read_tree(nr=2, file_path="' + tmp_file_path + '")'

            stoch_node_name, _, stoch_node_spec = \
                re.split(cmdu.assign_regex, cmd_line)
            cmd.parse_variable_assignment(
                dag_obj, stoch_node_name, stoch_node_spec, cmd_line)
            a_node_dag = dag_obj.get_node_dag_by_name("tr")

            self.assertEqual(a_node_dag.sample_size, 1)
            self.assertEqual(a_node_dag.repl_size, 2)
            self.assertEqual(len(a_node_dag.value), 2)

        finally:
            os.remove(tmp_file_path)

    def test_var_assignment_read_tree_sample_size_conflicts_with_later_stochastic_node(self):
        """'read_tree's own sample size is enforced against later '~' nodes.

        'read_tree' (assigned via '<-') is conceptually equivalent to
        a tree-valued stochastic node (e.g. 'discrete_sse') here, not
        to a plain literal constant (e.g. 'birth_rate <- 1.0') -- it
        has its own meaningful, explicit 'n'/'nr', so a later '~' node
        asking for a genuinely different (> 1) sample size must raise
        'DAGCannotAddNodeError', exactly as it already would between
        two ordinary stochastic nodes. A plain literal constant has
        no such function (and so no sample/replicate shape of its
        own) and remains exempt -- see
        'test_var_assignment_read_tree_sample_size_one_still_lenient'
        below for why 'n=1' itself does not trigger this (it does not
        for ordinary stochastic nodes either).
        """

        with tempfile.NamedTemporaryFile(
                mode="w", suffix=".tre", delete=False) as tmp_file:
            tmp_file.write("(tipa:1.0,tipb:1.0):0.0;\n")
            tmp_file.write("(tipc:2.0,tipd:2.0):0.0;\n")
            tmp_file_path = tmp_file.name

        try:
            script = (
                'trs <- read_tree(n=2, file_path="' + tmp_file_path + '")\n'
                'scalar_node ~ exponential(n=5, rate=1.0)\n')

            with self.assertRaises(ec.ScriptSyntaxError) as exc:
                cmd.script2dag(script, in_pj_file=False)

            self.assertIn(
                "Specified number of simulations was both > 1 and "
                "different from that in a previous command.",
                str(exc.exception))

        finally:
            os.remove(tmp_file_path)

    def test_var_assignment_read_tree_sample_size_one_still_lenient(self):
        """'read_tree's own 'n=1' does not lock out a later '~' node.

        Mirrors ordinary stochastic-node behavior: a sample size of
        exactly 1 is never treated as a hard commitment (only two
        *different*, both-greater-than-1 sizes conflict -- see
        'test_var_assignment_read_tree_sample_size_conflicts_with_later_stochastic_node'
        above) -- so a later '~' node is free to set its own, larger
        'n' without raising.
        """

        with tempfile.NamedTemporaryFile(
                mode="w", suffix=".tre", delete=False) as tmp_file:
            tmp_file.write("(tipa:1.0,tipb:1.0):0.0;\n")
            tmp_file_path = tmp_file.name

        try:
            script = (
                'trs <- read_tree(n=1, file_path="' + tmp_file_path + '")\n'
                'scalar_node ~ exponential(n=5, rate=1.0)\n')

            # must not raise
            dag_obj = cmd.script2dag(script, in_pj_file=False)
            scalar_node = dag_obj.get_node_dag_by_name("scalar_node")
            self.assertEqual(scalar_node.sample_size, 5)

        finally:
            os.remove(tmp_file_path)

    def test_var_assignment_read_tree_literal_constant_still_exempt(self):
        """A plain literal constant ('<-') stays exempt, unlike 'read_tree'.

        'birth_rate <- 1.0' has no constant *function* behind it (see
        'pgm.StochasticNodeDAG.constant_fn'), so it never participates
        in the DAG's sample-size bookkeeping at all -- a later '~'
        node is free to pick any 'n', same as before this behavior
        was added for 'read_tree'.
        """

        script = (
            'birth_rate <- 1.0\n'
            'scalar_node ~ exponential(n=5, rate=1.0)\n')

        # must not raise
        dag_obj = cmd.script2dag(script, in_pj_file=False)
        scalar_node = dag_obj.get_node_dag_by_name("scalar_node")
        self.assertEqual(scalar_node.sample_size, 5)


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
    # on the terminal, remember to add "src/" to
    # PYTHONPATH (system variable), or to set it if it does not
    # exist -- don't forget to export it!
    # 
    # Then you can do:
    # $ python3 tests/interface/test_cmd_var_assignment
    # 
    # or
    #
    # $ python3 -m tests.interface.test_cmd_var_assignment
    #
    # or 
    #
    # $ python3 -m unittest tests.interface.test_cmd_var_assignment.TestVarAssignment.test_var_assignment

    unittest.main()