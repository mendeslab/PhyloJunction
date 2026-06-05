import unittest
import re
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