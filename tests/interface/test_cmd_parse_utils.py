import unittest

# pj imports
import phylojunction.interface.cmdbox.cmd_parse_utils as cmdu

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
        
        token_dict1 = cmdu.tokenize_fn_spec(fn_spec1, cmd_line1)
        token_dict2 = cmdu.tokenize_fn_spec(fn_spec2, cmd_line2)
        token_dict3 = cmdu.tokenize_fn_spec(fn_spec3, cmd_line3)

        for k, v in self.tokenizer_exp_dict.items():
            self.assertEqual(v, token_dict1[k])

        for k, v in self.tokenizer_exp_dict.items():
            self.assertEqual(v, token_dict2[k])

        for k, v in self.tokenizer_exp_dict.items():
            self.assertEqual(v, token_dict3[k])

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