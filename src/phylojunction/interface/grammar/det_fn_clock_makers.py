import typing as ty

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.data.tree as pjtr
import phylojunction.utility.exception_classes as ec
import phylojunction.calculation.gamma_spike_clock_model as pjgsc

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


def make_gamma_spike_clock_tree(
        det_fn_name: str,
        det_fn_param_dict:
            ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
        -> ty.List[pjtr.AnnotatedTree]:
    """Build and return trees annotated by a Gamma-spike clock model.

    Called when the ':=' operator is used together with the
    'gamma_spike_clock' deterministic function, e.g.:

        n_br := n_rec_branches(tree=trs)
        br_specific_rel_rate ~ lognormal(n_nr_from=n_br, meanlog=-0.08, sdlog=0.4, log_space="true")
        phylo := gamma_spike_clock(tree=trs, attr_name="mutation_rate", clock_rate=0.01, br_specific_rel_rate=br_specific_rel_rate, spike_mean=0.01, spike_shape=2.0, strict_clock_when_no_stubs="true")

    There is no 'seed' parameter: reproducibility is controlled once,
    globally, for the whole DAG (see
    'calculation.gamma_spike_clock_model.populate_gamma_spike_mol_rates_on_rec_tr_branches'
    for why).

    'br_specific_rel_rate' is required: a DAG node already holding
    one pre-generated list of per-branch *relative* rates per tree in
    'tree' -- dimensionless multipliers, conventionally centered at
    1.0 in real space (e.g., a 'lognormal(meanlog=-sdlog**2/2, ...)'
    node, so that 'clock_rate' alone remains the mean absolute rate
    -- the same identifiability convention BEAST2's UCLN relaxed
    clock uses), e.g. built with 'lognormal(n_nr_from=n_br, ...)',
    where 'n_br := n_rec_branches(tree=trs)'. This model never draws
    the gradual, relaxed-clock component internally -- its own
    hyperparameters (and the distribution used) are always chosen
    independently, in the script. See
    'calculation.gamma_spike_clock_model.populate_gamma_spike_mol_rates_on_rec_tr_branches's
    'br_specific_rel_rates' parameter for exactly how these are
    consumed (the spike component is always still drawn internally
    here).

    'strict_clock_when_no_stubs' is optional (a quoted '"true"' or
    '"false"' string; defaults to '"false"' if not given) -- see
    'calculation.gamma_spike_clock_model.populate_gamma_spike_mol_rates_on_rec_tr_branches'
    for what it does. NOTE: setting it to '"true"' is NOT a realistic
    model of molecular evolution; it is a debugging/diagnostic device
    only (again, see that function's docstring), not something to use
    for actual analyses.

    Args:
        det_fn_name (str): Name of the deterministic function being
            called (always 'gamma_spike_clock').
        det_fn_param_dict (dict): Dictionary containing deterministic
            function parameter names (str) as keys and lists (of
            either strings or NodeDAGs) as values.

    Returns:
        (list of AnnotatedTree): One new tree per tree held by the
            'tree' argument's DAG node, each annotated by a
            'GammaSpikeClockModel'.
    """

    tree_dag_node: ty.Optional[pgm.NodeDAG] = None
    attr_name: ty.Optional[str] = None
    clock_rate: ty.Optional[float] = None
    br_specific_rel_rate_dag_node: ty.Optional[pgm.NodeDAG] = None
    spike_mean: ty.Optional[float] = None
    spike_shape: ty.Optional[float] = None
    strict_clock_when_no_stubs: bool = False

    def _extract_single_float(
            arg: str,
            val: ty.List[ty.Union[str, pgm.NodeDAG]]) -> float:
        """Extract a single float out of a (possibly DAG-node) argument.

        Vectorized clock hyperparameters (more than one value for
        'arg') are not yet supported -- see
        'GammaSpikeClockModel.init_check_vectorize_sample_size'.
        """

        if len(val) > 1:
            raise ec.ParseRequireSingleValueError(det_fn_name, arg)

        if isinstance(val[0], pgm.NodeDAG):
            cast_val = ty.cast(ty.List[pgm.NodeDAG], val)
            extracted = pgm.extract_value_from_dagnodes(cast_val)

            if len(extracted) != 1:
                raise ec.ParseRequireSingleValueError(det_fn_name, arg)

            try:
                return float(extracted[0])

            except (TypeError, ValueError):
                raise ec.ParseRequireNumericError(det_fn_name, arg)

        try:
            return float(ty.cast(str, val[0]))

        except ValueError:
            raise ec.ParseRequireNumericError(det_fn_name, arg)

    # val is a list of strings or nodes
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

        elif arg == "attr_name":
            if isinstance(val[0], pgm.NodeDAG):
                raise ec.ParseInvalidArgumentError(
                    "attr_name",
                    "'" + val[0].node_name + "'",
                    "Argument must be a quoted string naming the node "
                    "attribute to populate (e.g., \"mutation_rate\").")

            attr_name = ty.cast(str, val[0]).replace('"', "")

        elif arg == "clock_rate":
            clock_rate = _extract_single_float(arg, val)

        elif arg == "br_specific_rel_rate":
            if not isinstance(val[0], pgm.NodeDAG):
                raise ec.ParseInvalidArgumentError(
                    "br_specific_rel_rate",
                    "'" + str(val[0]) + "'",
                    "Argument must be a DAG node reference (e.g., a "
                    "'lognormal(n_nr_from=...)'-sampled one), not a "
                    "literal value.")

            br_specific_rel_rate_dag_node = val[0]

        elif arg == "spike_mean":
            spike_mean = _extract_single_float(arg, val)

        elif arg == "spike_shape":
            spike_shape = _extract_single_float(arg, val)

        elif arg == "strict_clock_when_no_stubs":
            if isinstance(val[0], pgm.NodeDAG):
                raise ec.ParseInvalidArgumentError(
                    "strict_clock_when_no_stubs",
                    "'" + val[0].node_name + "'",
                    "Argument must be \"true\" or \"false\".")

            val_str = ty.cast(str, val[0])

            if val_str not in ('"true"', '"false"'):
                raise ec.ParseInvalidArgumentError(
                    "strict_clock_when_no_stubs",
                    "'" + val_str + "'",
                    "Argument must be \"true\" or \"false\".")

            strict_clock_when_no_stubs = val_str == '"true"'

    ##############################
    # Required-parameter checks #
    ##############################

    if tree_dag_node is None:
        raise ec.ParseMissingParameterError("tree")

    if attr_name is None:
        raise ec.ParseMissingParameterError("attr_name")

    if clock_rate is None:
        raise ec.ParseMissingParameterError("clock_rate")

    if br_specific_rel_rate_dag_node is None:
        raise ec.ParseMissingParameterError("br_specific_rel_rate")

    if spike_mean is None:
        raise ec.ParseMissingParameterError("spike_mean")

    if spike_shape is None:
        raise ec.ParseMissingParameterError("spike_shape")

    clock_model = pjgsc.GammaSpikeClockModel(
        tree_dag_node=tree_dag_node,
        attr_name=attr_name,
        clock_rate=clock_rate,
        br_specific_rel_rate_dag_node=br_specific_rel_rate_dag_node,
        spike_mean=spike_mean,
        spike_shape=spike_shape,
        strict_clock_when_no_stubs=strict_clock_when_no_stubs)

    return clock_model.generate()
