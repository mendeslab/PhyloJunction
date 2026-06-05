import unittest

# pj imports
import phylojunction.utility.helper_functions as pjh
import phylojunction.calculation.math_utils as pjmath

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"

class TestExpectations(unittest.TestCase):

    def test_yule_height(self):
        """
        Test if expectation of Yule tree height
        """

        exp_h = pjmath.exp_root_height_yule_ntaxa(0.3, 7)
        self.assertAlmostEqual(exp_h, 5.309524, msg="Expected height should be 5.309524.", delta=1e-6)

    def test_bd_extant_count(self):
        """
        Test count of extant species of birth-death tree
        """

        exp_n = pjmath.exp_extant_count_bd(0.9, 0.7, 7.5)
        self.assertAlmostEqual(exp_n, 4.481689, msg="Expected count of extant taxa should be 4.481689.", delta=1e-6)

if __name__ == '__main__':
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
    # on the terminal, remember to add "src/phylojunction" to
    # PYTHONPATH (system variable), or to set it if it does not
    # exist -- don't forget to export it!
    # 
    # Then you can do:
    # $ python3 tests/distribution/test_math.py
    # 
    # or
    #
    # $ python3 -m tests.utility.test_math
    #
    # or 
    # 
    # $ python3 -m unittest tests.utility.test_math.TestExpectations.test_yule_height

    unittest.main()