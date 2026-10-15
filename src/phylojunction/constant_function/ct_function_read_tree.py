import os
import typing as ty

# dendropy imports for error checking
from dendropy.dataio.newickreader import NewickReader as dpnr
from dendropy.dataio.tokenizer import Tokenizer as dptk

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.utility.exception_classes as ec
import phylojunction.readwrite.pj_read as pjr
from phylojunction.data.tree import AnnotatedTree

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class CtFnTreeReader(pgm.ConstantFn):
    """Build one or more AnnotatedTree(s) from Newick string(s).

    'n'/'nr' ("number of samples"/"number of replicates") behave in
    one of exactly two ways, depending on whether the user wrote
    either of them in the script at all:

    (i) neither 'n' nor 'nr' is given: every tree string found (a
        'string' argument is always exactly one tree; a 'file_path'
        argument may hold any number of them, one per line) goes
        straight into 'n_samples' -- never into 'n_repl' -- so, e.g.,
        a 5-line file becomes 5 samples, 1 replicate each, with no
        further input needed from the user; or

    (ii) 'n' and/or 'nr' is given explicitly (even if the value
         written happens to be '1'): that value is taken completely
         literally, and checked -- via 'init_check_vectorize_sample_size'
         -- against how many tree strings were actually found
         ('n_samples * n_repl' must equal that count exactly).
         Mismatches raise 'exception_classes.ObjInitIncorrectDimensionError'
         rather than being silently reconciled one way or another.

    There is deliberately no third, in-between behavior (e.g., "read
    only the first tree of a multi-tree file") -- the whole file is
    always consumed, and 'n'/'nr' only ever either describe its
    exact shape (case ii) or get inferred from it outright (case i).
    """

    CT_FN_NAME = "read_tree"

    # both defaults to 1 in make_tree_reader()
    n_samples: int
    n_repl: int
    node_name_attr: str
    slice_age_ends: ty.List[float]
    slice_t_ends: ty.List[float]
    
    in_file: bool
    total_tr_count: int
    tr_str_list: ty.List[str]
    is_complete: bool
    keep_unsampled_alive: bool
    n_provided: bool
    nr_provided: bool

    # validation of pars happens in phylojunction_grammar
    def __init__(self,
                 n_samples: int,
                 n_repl: int,
                 tr_fp_or_str: str,
                 node_name_attr: str,
                 in_file: bool,
                 list_time_slice_age_ends: ty.Optional[ty.List[float]] = None,
                 is_complete: bool = False,
                 keep_unsampled_alive: bool = False,
                 n_provided: bool = False,
                 nr_provided: bool = False) -> None:

        self.n_samples = n_samples
        self.n_repl = n_repl
        self.node_name_attr = node_name_attr
        self.in_file = in_file
        self.is_complete = is_complete
        self.keep_unsampled_alive = keep_unsampled_alive
        self.n_provided = n_provided
        self.nr_provided = nr_provided

        self.slice_t_ends = None
        self.slice_age_ends = list_time_slice_age_ends

        if len(list_time_slice_age_ends) in (0, 1):
            # if 1 epoch is provided throw away last epoch end
            self.slice_age_ends = None

        # if there is 2 or more entries in time slice age ends,
        # throw away last epoch end
        elif list_time_slice_age_ends[-1] == 0.0:
            self.slice_age_ends = list_time_slice_age_ends[:-1]

        # will hold one or more tree Newick strings
        self.tr_str_list = list()
        if in_file:
            with open(tr_fp_or_str, "r") as tr_in:
                self.tr_str_list = tr_in.read().splitlines()
        
        else:
             self.tr_str_list = [tr_fp_or_str]

        self.total_tr_count = len(self.tr_str_list)

        # 'n'/'nr' are either both respected literally (whatever the
        # user wrote, even "n=1"/"nr=1"), or neither is given at all,
        # in which case every tree string found goes straight into
        # 'n_samples' (never 'n_repl') -- there is no silent middle
        # ground where a user-written "n=1" gets reinterpreted as
        # something else just because the file happens to hold more
        # than one tree. 'n_provided'/'nr_provided' (set in
        # 'interface.grammar.ct_fn_treereader_makers.make_tree_reader',
        # the only place that can actually tell "the user wrote this"
        # apart from "this is just the default") are what make that
        # distinction possible here.
        if not self.n_provided and not self.nr_provided:
            self.n_samples = self.total_tr_count

        self.init_check_vectorize_sample_size()

    def init_check_vectorize_sample_size(self) -> None:
        """Check provided tree strings and count members match.

        Always enforced -- when neither 'n' nor 'nr' was provided (see
        '__init__' above), 'n_samples' was already set to match
        'total_tr_count' exactly, so this can never raise in that
        case; it only ever does when the user gave an explicit 'n'
        and/or 'nr' that does not actually match how many tree strings
        were found.
        """

        arg: str = "string"
        if self.in_file:
            arg = "file_path"

        if (self.n_samples * self.n_repl) != self.total_tr_count:
            raise ec.ObjInitIncorrectDimensionError(
                    "CtFnTreeReader",
                    arg,
                    self.total_tr_count,
                    exp_len=self.n_samples * self.n_repl
                )

    def _update_ann_trs_time_slice_ends(self,
                                        generated_values:
                                        ty.List[AnnotatedTree]) -> None:

        for ann_tr in generated_values:
            seed_age_for_time_slicing = ann_tr.seed_age

            # convert into time ends
            # (0.0 at origin, seed_age = stop_val for "age" stop
            # condition in the present)
            if seed_age_for_time_slicing:
                # we make sure that the end of the present age is here
                # even if user already provided it (we forcefully removed
                # it at initialization)
                self.slice_age_ends.append(0.0)

                # initializes/resets slice_t_ends member
                self.slice_t_ends = list()

                # old first, young last
                n_slices_to_ignore = 0
                for age_end in self.slice_age_ends:
                    # we ignore user-specified age ends that for some reason
                    # older than the seed (origin/root) age; must have been
                    # a user oversight...
                    if age_end > seed_age_for_time_slicing:
                        n_slices_to_ignore += 1
                        continue

                    else:
                        self.slice_t_ends.append(seed_age_for_time_slicing - age_end)

            ann_tr.slice_age_ends = self.slice_age_ends
            ann_tr.slice_t_ends = self.slice_t_ends

    def get_rev_inference_spec_info(self) -> ty.List[str]:
        pass

    def generate(self) -> ty.List[AnnotatedTree]:
        generated_values: ty.List[AnnotatedTree] = list()
        
        try:
            for tr_str in self.tr_str_list:
                ann_tr = pjr.read_nwk_tree_str(
                    tr_str,
                    in_file=False,
                    node_names_attribute=self.node_name_attr,
                    is_complete=self.is_complete,
                    keep_unsampled_alive=self.keep_unsampled_alive)
                
                generated_values.append(ann_tr)

        except (TypeError,
                AttributeError,
                dptk.UnexpectedEndOfStreamError,
                dpnr.NewickReaderMalformedStatementError,
                dpnr.NewickReaderDuplicateTaxonError) as e:
            if not self.in_file:
                raise ec.ParseInvalidNewickStringError("string")
            
            else:
                raise ec.ParseInvalidNewickStringError("file_path")

        # if user provided time slice ends (in forward time)
        if self.slice_age_ends is not None:
            self._update_ann_trs_time_slice_ends(generated_values)
        
        return generated_values
            

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
    # $ python3 src/phylojunction/constant_function/ct_function_read_tree.py
    # 
    # or
    #
    # $ python3 -m phylojunction.constant_function.ct_function_read_tree

    n_samples = 1
    n_repl = 1

    tr1 = CtFnTreeReader(1, 1, "((A:1.0,B:1.0):1.0,C:2.0);", False)
    
    try:
        tr_error1 = CtFnTreeReader(1, 2, "((A:1.0,B:1.0):1.0,C:2.0);", False)
    except Exception as e:
        print(e)

    try:
        tr_error2 = CtFnTreeReader(2, 1, "((A:1.0,B:1.0):1.0,C:2.0);", False)
    except Exception as e:
        print(e)
    
    tr2 = CtFnTreeReader(1, 1, "examples/trees_maps_files/tree_to_read.tre", True)
    tr3 = CtFnTreeReader(1, 1, "examples/trees_maps_files/trees_to_read.tre", True)
    
    try:
        tr_error3 = CtFnTreeReader(1, 3, "examples/trees_maps_files/trees_to_read.tre", True)
    except Exception as e:
        print(e)

    try:
        tr_error4 = CtFnTreeReader(3, 1, "examples/trees_maps_files/trees_to_read.tre", True)
    except Exception as e:
        print(e)

    print(tr1.generate())
    print(tr2.generate())
    print(tr3.generate())