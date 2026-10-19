from phylojunction.utility.exception_classes import ScriptEvaluationError
import unittest
import re
from pathlib import Path
import tempfile

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.interface.cmdbox.cmd_parse as cmd
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestVarAssignment(unittest.TestCase):
    # Keep the existing file-reader assertions independent of absent example downloads.
    # The inputs retain the annotations and branch lengths asserted by those tests.
    @classmethod
    def setUpClass(cls):
        cls.tree_directory = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tree_directory.cleanup)
        cls.first_tree = str(Path(cls.tree_directory.name) / "one.tre")
        cls.two_trees = str(Path(cls.tree_directory.name) / "two.tre")
        first = "((sp1[&index=1]:1.0,sp2[&index=2]:1.0)[&index=4]:1.0,sp3[&index=3]:2.0)[&index=5];"
        second = "((sp1[&index=1]:0.1,sp2[&index=2]:0.1)[&index=4]:1.9,sp3[&index=3]:2.0)[&index=5];"
        Path(cls.first_tree).write_text(first + "\n")
        Path(cls.two_trees).write_text(first + "\n" + second + "\n")


    def test_var_assignment(self):
        """Assignments preserve scalar/vector shape and sequential name bindings."""
        result = cmd.compile_model('a <- 1\nb <- [1,2,3]').execute()
        self.assertEqual(result.value('a'), 1)
        self.assertEqual(result.value('b'), [1,2,3])
        self.assertEqual(result.model.n_nodes, 2)
        result = result.append('a <- [1]\nc <- b\nd <- [c,[4],[5],[6]]')
        self.assertEqual(result.value('a'), [1])
        self.assertEqual(result.value('c'), [1,2,3])
        self.assertEqual(result.value('d'), [[1,2,3],[4],[5],[6]])
        self.assertEqual(result.model.n_nodes, 5)
        self.assertEqual(len(result.model.name_node_dict), 4)

    def test_var_assignment_read_tree_string(self):
        """Test read_tree() from string.
        
        Test read_tree() calls using Newick strings directly produce
        the right DAG.
        """

        dag_obj = pgm.DirectedAcyclicGraph().execute()

        cmd_line1 = ('tr <- read_tree(string="((sp1[&index=1]:1.0,sp2'
                     '[&index=2]:1.0)[&index=4]:1.0,sp3[&index=3]:2.0)'
                     '[&index=5];", node_name_attr="index")')

        
        dag_obj = dag_obj.append(cmd_line1)
        a_node_dag = dag_obj.value("tr")

        self.assertEqual(1, len(dag_obj.model.name_node_dict))
        self.assertEqual(type(a_node_dag), list)
        self.assertEqual(len(a_node_dag), 1)
        self.assertEqual(a_node_dag[0].__str__().rstrip(),
                         ('((nd1:1.0[&index=1],nd2:1.0[&index=2])nd4:1.0'
                          '[&index=4],nd3:2.0[&index=3])nd5[&index=5];'))
        
        cmd_line2 = 'tr <- read_tree(string="((sp1:1.0,sp2:1.0):1.0,sp3:2.0);")'

        
        dag_obj = dag_obj.append(cmd_line2)
        a_node_dag = dag_obj.value("tr")

        self.assertEqual(1, len(dag_obj.model.name_node_dict))
        self.assertEqual(type(a_node_dag), list)
        self.assertEqual(len(a_node_dag), 1)
        self.assertEqual(a_node_dag[0].__str__().rstrip(),
                         ('((sp1:1.0,sp2:1.0)nd1:1.0,sp3:2.0)root;'))
        
    def test_var_assignment_read_tree_from_file(self):
        """Test read_tree() from file.
        
        Test read_tree()calls using tree files containing Newick
        strings directly produce the right DAG.
        """

        dag_obj = pgm.DirectedAcyclicGraph().execute()

        cmd_line1 = ('tr <- read_tree(file_path=\"' +
                     self.first_tree + '\",'
                     'node_name_attr="index")')

        
        dag_obj = dag_obj.append(cmd_line1)
        a_node_dag = dag_obj.value("tr")

        self.assertEqual(1, len(dag_obj.model.name_node_dict))
        self.assertEqual(type(a_node_dag), list)
        self.assertEqual(len(a_node_dag), 1)
        self.assertEqual(a_node_dag[0].__str__().rstrip(),
                         ('((nd1:1.0[&index=1],nd2:1.0[&index=2])nd4:1.0'
                          '[&index=4],nd3:2.0[&index=3])nd5[&index=5];'))
        
        cmd_line2 = ('tr <- read_tree(file_path=\"' +
                     self.two_trees + '\",'
                     'node_name_attr="index")')

        
        dag_obj = dag_obj.append(cmd_line2)
        a_node_dag = dag_obj.value("tr")

        self.assertEqual(1, len(dag_obj.model.name_node_dict))
        self.assertEqual(type(a_node_dag), list)
        self.assertEqual(len(a_node_dag), 2)
        self.assertEqual(a_node_dag[1].__str__().rstrip(),
                         ('((nd1:0.1[&index=1],nd2:0.1[&index=2])nd4:1.9'
                          '[&index=4],nd3:2.0[&index=3])nd5[&index=5];'))

        # The execution API preserves a tree collection within each run, rather than
        # interpreting file length as the number of model executions.
        for source, count in ((cmd_line1, 1), (cmd_line2, 2)):
            results = cmd.compile_model(source + '\nfirst := tr[1]').execute(2)
            self.assertEqual([len(v) for v in results.values('tr')], [count, count])
            for i in range(2):
                self.assertIs(results.value('first', i), results.value('tr', i)[0])

    def test_var_assignment_read_tree_string_exceptions(self):
        """Test read_tree() from string raised exceptions."""

        dag_obj = pgm.DirectedAcyclicGraph().execute()

        cmd_line1 = ('tr <- read_tree(n=2, string="((sp1[&index=1]:1.0,'
                     'sp2[&index=2]:1.0)[&index=4]:1.0,sp3[&index=3]:2.0)'
                     '[&index=5];", node_name_attr="index")')

        
        
        with self.assertRaises(ScriptEvaluationError) as exc:
            dag_obj = dag_obj.append(cmd_line1)

        self.assertIn('iid', str(exc.exception))

        cmd_line2 = ('tr <- read_tree(file_path="examples/trees_maps_files/'
                     'tree_to_read.tre", string="((sp1[&index=1]:1.0,'
                     'sp2[&index=2]:1.0)[&index=4]:1.0,sp3[&index=3]:2.0)'
                     '[&index=5];", node_name_attr="index")')

        

        with self.assertRaisesRegex(ScriptEvaluationError, 'exactly one'):
            dag_obj.append(cmd_line2)

    def test_var_assignment_read_tree_file_exceptions(self):
        """Test read_tree() from file raised exceptions."""

        dag_obj = pgm.DirectedAcyclicGraph().execute()

        cmd_line1 = ('tr <- read_tree(nr=13, file_path=\"' +
                     self.first_tree + '\", '
                     'node_name_attr="index")')

        

        with self.assertRaisesRegex(ScriptEvaluationError, 'iid'):
            dag_obj.append(cmd_line1)

        cmd_line2 = ('tr <- read_tree(file_path=\"' +
                     self.two_trees + '\", '
                     'node_name_attr="name")')

        

        expected_exception_message2 = \
            ("Could not initialize the object of \'read_tree\'. "
             "The argument provided to \'node_name_attribute\' was invalid. "
             "The attribute for node names was not in the newick strings.")
        
        with self.assertRaises(ScriptEvaluationError) as exc:
            dag_obj = dag_obj.append(cmd_line2)

        self.assertIsInstance(exc.exception.__cause__, ec.ObjInitInvalidArgError)
        self.assertEqual(str(exc.exception.__cause__),
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