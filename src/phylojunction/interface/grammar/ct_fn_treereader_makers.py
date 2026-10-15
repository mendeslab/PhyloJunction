import os
import typing as ty

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.utility.exception_classes as ec
from phylojunction.constant_function.ct_function_read_tree import CtFnTreeReader

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


def make_tree_reader(ct_fn_name: str,
                     ct_fn_param_dict: \
                        ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
        -> pgm.ConstantFn:
    """Parse 'read_tree's arguments and build a 'CtFnTreeReader'.

    This is the only place that can tell "the user wrote 'n=1'" apart
    from "the user did not write 'n' at all, so it defaulted to 1" --
    'ct_fn_param_dict' either has an 'n'/'nr' key or it does not,
    whereas by the time a plain 'int' value reaches 'CtFnTreeReader'
    itself, that distinction is already lost. 'n_provided'/
    'nr_provided' (below) capture it here and thread it through, so
    'CtFnTreeReader' can tell its own two cases apart -- see its own
    docstring for what each one does.

    Args:
        ct_fn_name (str): Name of the constant function being called
            (always 'read_tree').
        ct_fn_param_dict (dict): Dictionary containing constant
            function parameter names (str) as keys and lists (of
            either strings or NodeDAGs) as values.

    Returns:
        (ConstantFn): A 'CtFnTreeReader' instance, ready to
            '.generate()' one or more AnnotatedTree(s).
    """

    #############################
    # IMPORTANT: Default values #
    #############################
    n_samples: int = 1
    n_repl: int = 1
    # tracked separately from 'n_samples'/'n_repl' themselves, so that
    # 'CtFnTreeReader' can tell "the user wrote 'n=1'" apart from "the
    # user did not write 'n' at all, so it defaulted to 1" -- the two
    # cases are handled differently (see there)
    n_provided: bool = False
    nr_provided: bool = False
    tr_fp_or_str: str = str()
    node_name_attr: str = str()
    time_slice_age_ends: ty.List[float] = list()
    in_file: bool = False
    is_complete: bool = False
    keep_unsampled_alive: bool = False

    # args were already grammar-checked in constant_fn_grammar
    for arg, val in ct_fn_param_dict.items():
        # if element in val is string, it remains unchanged,
        # if NodeDAG, we get its string-fied value
        extracted_val: ty.List[str] = \
            pgm.extract_vals_as_str_from_node_dag(val)

        # if True:
        if arg == "n":
            if len(extracted_val) > 1:
                raise ec.ParseRequireSingleValueError("read_tree", arg)

            # only one element always
            try:
                n_samples = int(extracted_val[0])
                n_provided = True

            except ValueError:
                raise ec.ParseRequireIntegerError("read_tree", arg)

        elif arg == "nr":
            if len(extracted_val) > 1:
                raise ec.ParseRequireSingleValueError("read_tree", arg)

            # only one element always
            try:
                n_repl = int(extracted_val[0])
                nr_provided = True

            except ValueError:
                raise ec.ParseRequireIntegerError("read_tree", arg)

        elif arg == "file_path":
            if len(extracted_val) > 1:
                raise ec.ParseRequireSingleValueError("read_tree", arg)

            tr_fp_or_str = extracted_val[0].replace('"', '')

            if "string" in ct_fn_param_dict:
                raise ec.ParseMutuallyExclusiveParametersError(arg, "string")

            if not os.path.isfile(tr_fp_or_str):
                raise ec.ParsePathDoesNotExistError(
                    arg,
                    tr_fp_or_str)

            in_file = True

        elif arg == "string":
            if len(extracted_val) > 1:
                raise ec.ParseRequireSingleValueError("read_tree", arg)

            if "file_path" in ct_fn_param_dict:
                raise ec.ParseMutuallyExclusiveParametersError(arg, "file_path")

            tr_fp_or_str = extracted_val[0].replace('"', '')

        elif arg == "node_name_attr":
            node_name_attr = extracted_val[0].replace('"', '')

        elif arg == "epoch_age_ends":
            try:
                time_slice_age_ends = \
                    [float(v) for v in val if isinstance(v, str)]

            except ValueError:
                raise ec.ParseRequireNumericError(CtFnTreeReader.CT_FN_NAME,
                                                  arg)

        elif arg == "is_complete":
            if len(extracted_val) > 1:
                raise ec.ParseRequireSingleValueError("read_tree", arg)

            val_str = extracted_val[0]

            if val_str not in ('"true"', '"false"'):
                raise ec.ParseInvalidArgumentError(
                    "is_complete",
                    "'" + val_str + "'",
                    "Argument must be \"true\" or \"false\".")

            is_complete = val_str == '"true"'

        elif arg == "keep_unsampled_alive":
            if len(extracted_val) > 1:
                raise ec.ParseRequireSingleValueError("read_tree", arg)

            val_str = extracted_val[0]

            if val_str not in ('"true"', '"false"'):
                raise ec.ParseInvalidArgumentError(
                    "keep_unsampled_alive",
                    "'" + val_str + "'",
                    "Argument must be \"true\" or \"false\".")

            keep_unsampled_alive = val_str == '"true"'

    return CtFnTreeReader(n_samples,
                          n_repl,
                          tr_fp_or_str,
                          node_name_attr,
                          in_file,
                          time_slice_age_ends,
                          is_complete=is_complete,
                          keep_unsampled_alive=keep_unsampled_alive,
                          n_provided=n_provided,
                          nr_provided=nr_provided)