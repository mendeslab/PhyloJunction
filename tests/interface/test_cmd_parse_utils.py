import unittest
import ast
import random
import numpy as np
from phylojunction.interface.cmdbox.cmd_parse import compile_model
from phylojunction.interface.cmdbox import expressions

# pj imports

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"

class TestCommandParseUtils(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tokenizer_exp_dict = { "name": "\"lambda0\"", "value": "1.0", "event":"\"w_speciation\"", "states":"[0,0,0]" }

    def test_tokenizer(self):
        """
        Test if tokenizer produces the correct dictionary outputs
        """
        fn_spec1 = "name=\"lambda0\", value=1.0, event=\"w_speciation\", states=[0, 0, 0]"
        fn_spec2 = "name=\"lambda0\", states=[0, 0, 0], value=1.0, event=\"w_speciation\""
        fn_spec3 = "states=[0, 0, 0], name=\"lambda0\", value=1.0, event=\"w_speciation\""
        
        cmd_line1 = "a := sse_rate(" + fn_spec1 + ")"    
        cmd_line2 = "a := sse_rate(" + fn_spec2 + ")" 
        cmd_line3 = "a := sse_rate(" + fn_spec3 + ")"     
        
        token_dict1 = {kw.arg: ast.literal_eval(kw.value) for kw in expressions.parse_expression('sse_rate(' + fn_spec1 + ')').keywords}
        token_dict2 = {kw.arg: ast.literal_eval(kw.value) for kw in expressions.parse_expression('sse_rate(' + fn_spec2 + ')').keywords}
        token_dict3 = {kw.arg: ast.literal_eval(kw.value) for kw in expressions.parse_expression('sse_rate(' + fn_spec3 + ')').keywords}

        for k, v in self.tokenizer_exp_dict.items():
            self.assertEqual(ast.literal_eval(v), token_dict1[k])

        for k, v in self.tokenizer_exp_dict.items():
            self.assertEqual(ast.literal_eval(v), token_dict2[k])

        for k, v in self.tokenizer_exp_dict.items():
            self.assertEqual(ast.literal_eval(v), token_dict3[k])

    # The old flat-token tests cannot protect nested syntax or dependency discovery.
    # Keep this while expressions are shared by compilation and constructor dispatch.
    def test_nested_expressions(self):
        _, _, expression = expressions.parse_statement('x ~ iid(n,normal(matrix[i][j],1))')
        self.assertEqual(list(expressions.variable_names(expression)), ['n', 'matrix', 'i', 'j'])
        for source in ('x.real', '__import__("os").system("pwd")', '[x for x in y if true]', 'x[1:2]'):
            with self.assertRaises(ValueError):
                expressions.parse_expression(source)
        call = expressions.parse_expression('f(name="a b", value=[1, 2])')
        self.assertEqual({kw.arg: ast.literal_eval(kw.value) for kw in call.keywords},
                         {'name': 'a b', 'value': [1, 2]})

    # Shared shape dispatch must distinguish list arguments from parameter mapping.
    # Scalar constructor tests cannot protect this boundary; retain while mapping is implicit.
    def test_vectorized_calls(self):
        from phylojunction.interface.grammar.constructors import construct
        ds = construct('normal', [0, 2], 1)
        self.assertEqual([d.mean for d in ds], [0, 2])
        self.assertEqual(construct('normal', [], -1), [])
        self.assertEqual(construct('categorical', [1, 2], [1, 1]).cats, [1, 2])
        self.assertEqual(len(construct('sse_rate', ['a', 'b'], [1, 2], 'speciation', [0, 0, 0])), 2)
        for args in (([0], [1, 2]), ([[0]], 1)):
            with self.assertRaisesRegex(ValueError, 'lengths|shapes'):
                construct('normal', *args)
        with self.assertRaises(TypeError):
            construct('normal', [], bogus=1)

    # End-to-end comparisons protect local scope and recursive sampling, which scalar
    # constructor tests cannot check. Retain while both implicit and explicit mapping exist.
    def test_comprehensions_and_joint(self):
        prefix = 'mu <- 99\nmus <- [0,2]\nsds <- [1,3]\n'
        a = compile_model(prefix + 'xs ~ normal(mus,sds)').execute(2, 42)
        b = compile_model(prefix + 'xs ~ [normal(mu,sd) for (mu,sd) in zip_strict(mus,sds)]').execute(2,42)
        self.assertEqual(a.values('xs'), b.values('xs'))
        self.assertEqual(b.model.name_node_dict['xs'].bindings, {'mus': 1, 'sds': 2})
        nested = compile_model('rows <- [[0],[1,2]]\nx ~ [[normal(mu,1) for mu in row] for row in rows]').execute(1,42)
        self.assertEqual([len(row) for row in nested.value('x')], [1,2])
        self.assertEqual(compile_model('xs ~ normal([], -1)').execute().value('xs'), [])
        self.assertEqual(compile_model('xs <- [1,2.5]').execute().value('xs'), [1,2.5])
        for source in ('x <- [1,"a"]', 'x <- [1,[2]]', 'x ~ [normal(0,1),2]',
                       'x := [a for (a,b) in [[1]]]', 'x := zip_strict([1],[1,2])',
                       'x ~ [[normal(0,1,clamp=true)]]'):
            with self.assertRaises((TypeError, ValueError)):
                compile_model(source).execute()

    # The collection mapper can expand histories without the PJ wrapper returning them.
    # Check the public expression result as well as shared-input immutability.
    def test_history_values(self):
        from pathlib import Path
        source = Path('examples/see_stoch_maps.pj').read_text()
        result = compile_model(source).execute()
        self.assertEqual(len(result.value('histories')), 91)
        self.assertEqual(len(result.value('mapped_trs')), 91)
        self.assertIsNot(result.value('mapped_trs')[0], result.value('mapped_trs')[1])
        import tempfile
        from phylojunction.readwrite.execution_export import write_maps
        with tempfile.TemporaryDirectory() as directory:
            write_maps(directory, result, ['mapped_trs'], 'state')
            self.assertEqual(len(list(Path(directory).glob('*_smaps.tsv'))), 91)

        original = result.value('trs')[0]
        self.assertNotEqual(original.at_dict, result.value('mapped_trs')[0].at_dict)
        paired = result.append('two := map_attr(tree=[trs[1],trs[1]], history=[histories[1],histories[2]])')
        self.assertEqual(len(paired.value('two')), 2)
        with self.assertRaisesRegex(ValueError, 'equal lengths'):
            result.append('bad := map_attr(tree=trs, history=histories)')

    # Legacy batch tests do not exercise shared definitions or transactional additions.
    # These invariants remain necessary while executions and IID values are distinct.
    def test_model_executions(self):
        model = compile_model('d := normal(1,2)\nx ~ iid(3,iid(2,d))\ni <- 2\ny := x[i][1]')
        result = model.execute(3, 42)
        self.assertEqual(result.n_executions, 3)
        self.assertEqual(len(result.model.definitions), 4)
        self.assertEqual(result.model.name_node_dict['y'].bindings, {'x': 1, 'i': 2})
        for execution in range(3):
            self.assertEqual(result.value('y', execution), result.value('x', execution)[1][0])
        self.assertEqual(result.values('x'), model.execute(3, 42).values('x'))
        other = result.append('x <- 5\nz := x')
        self.assertEqual(other.values('z'), [5, 5, 5])
        self.assertEqual(other.values('y'), result.values('y'))
        self.assertNotIn('z', result.model.name_node_dict)
        self.assertFalse(hasattr(result.model.name_node_dict['x'], 'value'))
        model.add_definition('extra <- 0')
        self.assertNotIn('extra', result.model.name_node_dict)
        for source in ('v := x[0]', 'v := x[4]', 'v := x[1.0]', 'v := x[true]', 'v ~ 2'):
            with self.assertRaises((IndexError, TypeError)):
                result.append(source)
        before_py, before_np = random.getstate(), np.random.get_state()
        with self.assertRaises(IndexError):
            result.append('new ~ normal(0,1)\nv := x[4]')
        self.assertEqual(before_py, random.getstate())
        after_np = np.random.get_state()
        np.testing.assert_array_equal(before_np[1], after_np[1])
        self.assertEqual(before_np[2:], after_np[2:])
        with self.assertRaisesRegex(ValueError, 'iid'):
            compile_model('x ~ normal(n=2, mean=0, sd=1)').execute()


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
    # $ python3 tests/interface/test_cmd_parse_utils.py
    # 
    # or
    #
    # $ python3 -m tests.interface.test_cmd_parse_utils
    #
    # or 
    #
    # $ python3 -m unittest tests.interface.test_cmd_parse_utils.TestCommandParseUtils.test_tokenizer

    unittest.main()