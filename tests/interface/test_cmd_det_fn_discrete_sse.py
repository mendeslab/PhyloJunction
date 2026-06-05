import unittest
import re

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.interface.cmdbox.cmd_parse as cmd
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestDetFnDiscreteSSEErrors(unittest.TestCase):
    def test_rate_epoch_index_too_large_error(self):
        """
        Test if rate annotated with epoch index that is too large
        (e.g., 2 epochs, index provided is 3)
        """
        
        yule_timehet_invalid_epoch = \
            ("l0ratet1 := sse_rate(name=\"lambdat1\", value=1.0, event=\"w_speciation\", epoch=1)\n"
             "l0ratet2 := sse_rate(name=\"lambdat2\", value=0.25, event=\"w_speciation\", epoch=2)\n"
             "l0ratet3 := sse_rate(name=\"lambdat3\", value=3.0, event=\"w_speciation\", epoch=3)\n"
             "stash := sse_stash(flat_rate_mat=[l0ratet1, l0ratet2, l0ratet3], n_states=1, n_epochs=2, "
             "seed_age=3.0, epoch_age_ends=[1.2, 0.7])")

        with self.assertRaises(ec.ParseDetFnInitFailError) as exc:
            dag_obj = cmd.script2dag(yule_timehet_invalid_epoch, in_pj_file=False)

        self.assertEqual(
            str(exc.exception),
            ("Deterministic output from 'sse_stash' could not be instantiated. "
             "Incorrect dimension of container epoch_age_ends, which was of size 2. "
             "The expected dimension was 1."))
        
    def test_rate_epoch_index_zero(self):
        """
        Test if rate annotated with epoch index that is zero
        (e.g., user thinks first epoch has index 0) causes error
        """
        
        yule_zero_index = \
            ("l0ratet1 := sse_rate(name=\"lambdat1\", value=1.0, event=\"w_speciation\", epoch=0)\n"
             "l0ratet2 := sse_rate(name=\"lambdat2\", value=0.25, event=\"w_speciation\", epoch=1)\n"
             "stash := sse_stash(flat_rate_mat=[l0ratet1, l0ratet2], n_states=1, n_epochs=2, "
             "seed_age=3.0, epoch_age_ends=[0.7])")

        with self.assertRaises(ec.ParseRequirePositiveIntegerError) as exc:
            dag_obj = cmd.script2dag(yule_zero_index, in_pj_file=False)

        self.assertEqual(
            str(exc.exception),
            ("When specifying object sse_rate's parameter 'epoch', something other than a positive "
             "integer was provided. A positive integer is required."))

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
    # $ python3 tests/interface/test_cmd_det_fn_discrete_sse
    # 
    # or
    #
    # $ python3 -m tests.interface.test_cmd_det_fn_discrete_sse
    #
    # or 
    #
    # $ python3 -m unittest tests.interface.test_cmd_det_fn_discrete_sse.TestDetFnDiscreteSSEErrors.test_rate_epoch_annotation_error

    unittest.main()