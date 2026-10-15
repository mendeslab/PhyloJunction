import typing as ty

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.data.tree as pjtr
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


def make_n_rec_branches(
        det_fn_name: str,
        det_fn_param_dict:
            ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
        -> ty.List[int]:
    """Build and return the number of branches of each tree in 'tree'.

    Called when the ':=' operator is used together with the
    'n_rec_branches' deterministic function, e.g.:

        n_br := n_rec_branches(tree=trs)

    One int per tree held by the 'tree' argument's DAG node, in the
    same (sample, replicate) order, counting only real branches of
    each tree's *reconstructed* topology -- see 'tree.is_real_branch'
    for exactly what that excludes (the root, and any zero-length
    branch), and 'AnnotatedTree.n_rec_branches' for the counting itself.
    Because this is pure tree-topology bookkeeping (no randomness
    involved), the result can in principle differ for every
    (sample, replicate) tree, even within the same sample.

    This is meant to be usable on its own (e.g., for summary
    statistics), but is also the building block for feeding
    per-branch base rates into a clock model -- see
    'interface.grammar.dn_grammar's 'n_nr_from' argument, which can
    point a parametric distribution (e.g. 'lognormal') at this det-fn's
    output so it draws one value per branch, per tree.

    Args:
        det_fn_name (str): Name of the deterministic function being
            called (always 'n_rec_branches').
        det_fn_param_dict (dict): Dictionary containing deterministic
            function parameter names (str) as keys and lists (of
            either strings or NodeDAGs) as values.

    Returns:
        (list of int): One branch count per tree in 'tree', in the
            same order.
    """

    tree_dag_node: ty.Optional[pgm.NodeDAG] = None

    for arg, val in det_fn_param_dict.items():
        if not val:
            raise ec.ParseMissingArgumentError(arg)

        if arg == "tree":
            if not isinstance(val[0], pgm.NodeDAG):
                raise ec.ParseInvalidArgumentError(
                    "tree",
                    "'" + str(val[0]) + "'",
                    "Argument must be a DAG node holding a phylogenetic "
                    "tree.")

            tree_dag_node = val[0]

    if tree_dag_node is None:
        raise ec.ParseMissingParameterError("tree")

    return [
        ann_tr.n_rec_branches() for ann_tr in
        ty.cast(ty.List[pjtr.AnnotatedTree], tree_dag_node.value)
    ]
