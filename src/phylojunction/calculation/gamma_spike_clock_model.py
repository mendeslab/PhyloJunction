from __future__ import annotations

import math
import typing as ty
import numpy as np

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.data.tree as pjtr
import phylojunction.calculation.clock_model as pjclock

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


def populate_gamma_spike_mol_rates_on_rec_tr_branches(
        ann_tr: pjtr.AnnotatedTree,
        clock_rate: float,
        spike_mean: float,
        spike_shape: float,
        br_specific_rel_rates: ty.List[float],
        strict_clock_when_no_stubs: bool = False) -> None:
    """Populate molecular-rate attrs on rec. tree branches (Gamma-spike clock).

    Implements the punctuated ('Gamma-spike') relaxed-clock model of
    Douglas et al. (2025): each branch's molecular rate is a gradual
    (log-normal relaxed-clock) component plus an abrupt ('spike')
    component whose size depends on the number of stubs (unsampled
    bifurcations) on that branch.

    Under the model exactly as described by Douglas et al.
    ('strict_clock_when_no_stubs=False', the default), a branch leaving
    a genuine speciation with zero stubs still gets independent gradual and
    spike components (a relaxed clock varies rates branch-to-branch
    regardless of stubs; zero stubs just means one baseline spike --
    see eqn 2.7 -- rather than a larger, stub-inflated one). A branch
    leaving a sampled-ancestor observation has no baseline spike. Setting
    'strict_clock_when_no_stubs=True' deliberately departs from the
    paper: a branch with zero stubs instead gets a fixed, deterministic
    rate (equivalent to a strict clock: 'branch_raw_rate = 1.0',
    'spike_effect = 0.0', so 'mol_rate' is exactly 'clock_rate' for
    every such branch), and only a branch with at least one stub keeps
    the usual relaxed-clock-plus-spike draw.

    NOTE: 'strict_clock_when_no_stubs=True' is NOT a realistic
    alternative model of molecular evolution, and should not be used
    to represent a biological hypothesis -- pinning every stub-free
    branch to a fixed rate has no basis in Douglas et al. (2025) or
    elsewhere. It exists purely as a debugging/diagnostic device: by
    making every stub-free branch identical (e.g., the same color, if
    colored by 'mol_rate'), it isolates the spike component's effect
    from ordinary relaxed-clock rate variation, which is useful for
    probing/visualizing how stubs alone affect rates, or for sanity-
    checking a clock-model implementation, but it is not meant to be
    used as a clock model for actual analyses.

    This function is deliberately independent of PhyloJunction's DAG
    and grammar machinery (see 'GammaSpikeClockModel' below for that):
    it operates on a single 'AnnotatedTree' instance, in place, and
    can be imported and used on its own.

    Random draws are made with NumPy's legacy, global random number
    generator (module-level 'numpy.random.*' functions), the same
    convention used throughout PhyloJunction (e.g.,
    'distribution.dn_discrete_sse.DnSSE'). There is no 'seed'
    parameter here: reproducibility is controlled once, globally, by
    whichever caller already seeded that generator -- the DAG's own
    'random_seed' (set from a '.pj' script's '-r' CLI flag or the
    GUI's seed field; see 'pgm.DirectedAcyclicGraph.random_seed').
    Callers that are not going through the DAG (e.g., unit tests)
    must call 'numpy.random.seed(...)' themselves beforehand for
    reproducible draws.

    Side-effect: populates the attributes 'branch_raw_rate' (the
    log-normal relaxed-clock draw), 'spike_effect' (the Gamma-spike
    draw), and 'mol_rate' (their combination, see below) on every
    node of 'ann_tr.tree_reconstructed' with a non-zero-length
    branch. The root, and any sampled-ancestor branch (both of which
    have no meaningful rate of their own), get 'None' for all three.

    This function also (re-)populates the 'n_stub' attribute these
    branches require, by calling
    'AnnotatedTree.populate_n_stubs_on_branches_of_rec_tr' itself --
    counting stubs is pure tree-topology bookkeeping, shared by any
    future stub-aware clock model, and cheap/idempotent to redo, so
    callers do not need to remember to call it first.

    References:
        Douglas, J., Bouckaert, R., Harris, S.C., Carter Jr, C.W.,
        Wills, P.R. (2025) Evolution is coupled with branching
        across many granularities of life. Proceedings of the Royal
        Society B, 292: 20250182.
        https://doi.org/10.1098/rspb.2025.0182

    Args:
        ann_tr (AnnotatedTree): Tree instance whose reconstructed
            tree's branches will be annotated, in place.
        clock_rate (float): Mean clock (substitution) rate (mu_c in
            Douglas et al. 2025). Must be finite and > 0.0.
        spike_mean (float): Mean spike size (S_mu). Must be finite and >= 0.0.
        spike_shape (float): Shape (S_alpha) of the Gamma
            distribution governing spike sizes (eqn 2.7). Must be
            finite and > 0.0.
        br_specific_rel_rates (list of float): Pre-generated
            per-branch *relative* rates (the gradual, relaxed-clock
            component), e.g. produced by a separate
            'distribution.dn_parametric.DnLogNormal' node via the
            grammar's 'n_nr_from' argument pointed at an
            'n_rec_branches' det-fn node for this same tree (see
            'interface.grammar.dn_grammar'). Conventionally centered
            at 1.0 in real space (e.g., a log-normal with
            'meanlog = -sdlog**2 / 2'), so that 'clock_rate' alone
            remains the mean absolute rate -- the same
            identifiability convention BEAST2's UCLN relaxed clock
            uses -- though this function itself does not enforce that
            centering; it simply multiplies whatever is given by
            'clock_rate'. Entries must be finite and nonnegative. Must
            have exactly 'ann_tr.n_rec_branches()'
            entries, one per real branch of the reconstructed tree, in
            the same (preorder) order 'n_rec_branches' counts them in
            -- checked directly against the tree's own actual real-
            branch count, regardless of how 'br_specific_rel_rates'
            was itself produced. The spike component is always still
            drawn internally here.
        strict_clock_when_no_stubs (bool): See above -- NOT a
            realistic model, a debugging/diagnostic device only.
            Defaults to 'False' (the model as described by Douglas et
            al.).

    Raises:
        ValueError: Is raised if any of 'clock_rate', 'spike_mean',
            or 'spike_shape' is out of its valid range, or if
            'br_specific_rel_rates' contains invalid values or has the
            wrong number of entries.
    """

    if not math.isfinite(clock_rate) or clock_rate <= 0.0:
        raise ValueError("clock_rate must be finite and greater than 0.0")

    if not math.isfinite(spike_mean) or spike_mean < 0.0:
        raise ValueError("spike_mean must be finite and nonnegative")

    if not math.isfinite(spike_shape) or spike_shape <= 0.0:
        raise ValueError("spike_shape must be finite and greater than 0.0")

    # Validate the entire vector before annotating any branch or drawing any spikes.
    # Zero relative rates are valid: the gradual component can be absent on an edge.
    if any(not math.isfinite(rate) or rate < 0.0 for rate in br_specific_rel_rates):
        raise ValueError("br_specific_rel_rates must all be finite and nonnegative")

    if not ann_tr.tree_reconstructed:
        ann_tr.extract_reconstructed_tree(plotting_overhead=False)

    n_rec_branches = ann_tr.n_rec_branches()

    if len(br_specific_rel_rates) != n_rec_branches:
        raise ValueError(
            "'br_specific_rel_rates' had " +
            str(len(br_specific_rel_rates)) +
            " entries, but the reconstructed tree has " +
            str(n_rec_branches) + " real branches -- they must match.")

    # idempotent (recomputed from tree topology, no randomness involved),
    # so it is safe to call again even if the caller already has
    ann_tr.populate_n_stubs_on_branches_of_rec_tr()

    # advances once per real branch (see 'pjtr.is_real_branch'),
    # independently of 'strict_clock_when_no_stubs' overriding that
    # branch's rate below -- so it always lines up with
    # 'br_specific_rel_rates', which 'n_rec_branches()' (above) sizes
    # the same way
    branch_idx = 0

    for rec_nd in ann_tr.tree_reconstructed.preorder_node_iter():
        branch_time = rec_nd.edge_length

        # if origin or branch subtending node is effectively
        # 0.0 (direct ancestor node), then ignore this node and
        # move to the next -- see 'pjtr.is_real_branch' for exactly
        # which branches this excludes
        if not pjtr.is_real_branch(rec_nd, ann_tr.epsilon):
            rec_nd.branch_raw_rate = None
            rec_nd.spike_effect = None
            rec_nd.mol_rate = None
            continue

        this_branch_idx = branch_idx
        branch_idx += 1

        # deliberate departure from Douglas et al. -- see this
        # function's docstring; a stub-free branch gets a fixed,
        # strict-clock rate instead of a relaxed-clock-plus-spike
        # draw, and consumes no entry of 'br_specific_rel_rates'
        if strict_clock_when_no_stubs and rec_nd.n_stub == 0:
            rec_nd.branch_raw_rate = 1.0
            rec_nd.spike_effect = 0.0
            rec_nd.mol_rate = clock_rate
            continue

        # branch-specific relative rate, externally provided (e.g.,
        # by a 'lognormal(n_nr_from=...)' DAG node)
        rec_nd.branch_raw_rate = br_specific_rel_rates[this_branch_idx]

        # Number of Gamma-distributed spikes: one per hidden split on this segment,
        # plus one for the PARENT if it is a genuine speciation, not an SA marker.
        # Whether this node itself bifurcates is irrelevant to its incoming edge;
        # that bifurcation contributes to each of its daughters' edges instead.
        # This also covers a genuine root reached through a chain of SA observations.
        # Inspect reconstructed topology directly: suppression can leave lineage flags stale.
        # Douglas et al. (2025), eqn 2.7, attribute the extra spike to the parent;
        # likewise, getNSpikes() in GammaSpikeModel's BranchSpikePrior.java checks
        # node.getParent().isFake(), not the node's own children.
        n_spikes = rec_nd.n_stub
        if not rec_nd.parent_node.is_sa_dummy_parent:
            n_spikes += 1

        # if there are no spikes to add (no stubs, and no baseline
        # spike for this branch), the spike effect is 0.0 -- a
        # Gamma distribution with shape 0.0 is a degenerate point
        # mass at 0.0, so we special-case it rather than draw from
        # rng.gamma()
        if n_spikes == 0:
            rec_nd.spike_effect = 0.0

        else:
            gamma_shape = spike_shape * n_spikes
            gamma_scale = 1.0 / spike_shape

            spike_effect_prime = np.random.gamma(
                shape=gamma_shape,
                scale=gamma_scale
            )

            # recall from Douglas et al.
            # spike_mean = S_mu (e.g., normally distributed around 0.01)

            # final spike effect
            rec_nd.spike_effect = spike_mean * spike_effect_prime

        # (mu_c * r_i * branch_time) + spike_effect

        # final mol rate for this branch
        rec_nd.mol_rate = (
            clock_rate * rec_nd.branch_raw_rate
            + rec_nd.spike_effect / branch_time
            # we divide spike_effect by branch_time because the
            # spike effect is added to the final rate (see below)
            #
            # ((mu_c * r_i) + spike_effect) * branch_time =
            # mu_c * branch_time + r_i * branch_time + spike_effect * branch_time
        )


class GammaSpikeClockModel(pjclock.ClockModel):
    """Punctuated ('Gamma-spike') relaxed molecular clock model.

    Each branch's molecular rate is the sum of a gradual component
    (an i.i.d. log-normal relaxed clock) and an abrupt, 'spike'
    component, whose size is drawn from a Gamma distribution
    conditioned on the number of stubs (unsampled bifurcations) on
    that branch. See 'populate_gamma_spike_mol_rates_on_rec_tr_branches'
    (this module) for the underlying algorithm, which this class
    wraps as a PhyloJunction clock model.

    References:
        Douglas, J., Bouckaert, R., Harris, S.C., Carter Jr, C.W.,
        Wills, P.R. (2025) Evolution is coupled with branching
        across many granularities of life. Proceedings of the Royal
        Society B, 292: 20250182.
        https://doi.org/10.1098/rspb.2025.0182
    """

    CLOCK_MODEL_NAME = "gamma_spike_clock"

    clock_rate: float
    spike_mean: float
    spike_shape: float
    strict_clock_when_no_stubs: bool
    br_specific_rel_rate_dag_node: pgm.NodeDAG

    def __init__(
            self,
            *,
            tree_dag_node: pgm.NodeDAG,
            attr_name: str,
            clock_rate: float,
            spike_mean: float,
            spike_shape: float,
            br_specific_rel_rate_dag_node: pgm.NodeDAG,
            strict_clock_when_no_stubs: bool = False) -> None:
        """Store Gamma-spike clock-model hyperparameters.

        Hyperparameter validity (e.g., 'clock_rate' > 0.0) is checked
        when 'populate_branch_rates()' is called (via
        'populate_gamma_spike_mol_rates_on_rec_tr_branches'), not
        here -- except for 'br_specific_rel_rate_dag_node's own shape,
        checked eagerly below, since a mismatch there reflects a
        DAG-wiring mistake rather than an ordinary bad hyperparameter
        value.

        There is no 'seed' parameter -- see
        'populate_gamma_spike_mol_rates_on_rec_tr_branches' (this
        module) for why.

        Args:
            tree_dag_node (pgm.NodeDAG): See 'ClockModel.__init__'.
            attr_name (str): See 'ClockModel.__init__'.
            clock_rate (float): Mean clock (substitution) rate
                (mu_c in Douglas et al. 2025).
            spike_mean (float): Mean spike size (S_mu).
            spike_shape (float): Shape (S_alpha) of the Gamma
                distribution governing spike sizes (eqn 2.7).
            br_specific_rel_rate_dag_node (pgm.NodeDAG): A DAG node
                (e.g., a 'lognormal(n_nr_from=...)'-sampled one) whose
                value holds one pre-generated list of per-branch
                *relative* rates per tree in 'tree_dag_node', in the
                same (sample, replicate) order -- see
                'populate_gamma_spike_mol_rates_on_rec_tr_branches's
                'br_specific_rel_rates' parameter for exactly how
                these are consumed, and for the conventional (but not
                enforced) centering at 1.0 in real space. Must have
                the same 'sample_size'/'repl_size' as 'tree_dag_node'.
                This model does not draw the gradual, relaxed-clock
                component internally -- it must always be supplied
                this way.
            strict_clock_when_no_stubs (bool): See
                'populate_gamma_spike_mol_rates_on_rec_tr_branches'
                (this module) -- NOT a realistic model, a
                debugging/diagnostic device only. Defaults to 'False'
                (the model as described by Douglas et al.).

        Raises:
            ValueError: Is raised if 'br_specific_rel_rate_dag_node'
                does not have the same 'sample_size'/'repl_size' as
                'tree_dag_node'.
        """

        super().__init__(tree_dag_node=tree_dag_node, attr_name=attr_name)

        if br_specific_rel_rate_dag_node.sample_size != \
                tree_dag_node.sample_size \
                or br_specific_rel_rate_dag_node.repl_size != \
                tree_dag_node.repl_size:
            raise ValueError(
                "'br_specific_rel_rate_dag_node' has sample_size=" +
                str(br_specific_rel_rate_dag_node.sample_size) +
                ", repl_size=" +
                str(br_specific_rel_rate_dag_node.repl_size) +
                ", but 'tree_dag_node' has sample_size=" +
                str(tree_dag_node.sample_size) + ", repl_size=" +
                str(tree_dag_node.repl_size) +
                " -- they must match, one relative-rate vector per "
                "tree.")

        self.clock_rate = clock_rate
        self.spike_mean = spike_mean
        self.spike_shape = spike_shape
        self.strict_clock_when_no_stubs = strict_clock_when_no_stubs
        self.br_specific_rel_rate_dag_node = br_specific_rel_rate_dag_node

    def populate_branch_rates(
            self, ann_tr: pjtr.AnnotatedTree, tree_idx: int) -> None:
        """See 'ClockModel.populate_branch_rates'.

        In addition to 'self.attr_name', this model also populates
        'branch_raw_rate' and 'spike_effect' on every non-zero-length
        reconstructed-tree branch -- see
        'populate_gamma_spike_mol_rates_on_rec_tr_branches' (this
        module) for their exact definitions.
        """

        br_specific_rel_rates = \
            self.br_specific_rel_rate_dag_node.value[tree_idx]

        populate_gamma_spike_mol_rates_on_rec_tr_branches(
            ann_tr,
            clock_rate=self.clock_rate,
            spike_mean=self.spike_mean,
            spike_shape=self.spike_shape,
            br_specific_rel_rates=br_specific_rel_rates,
            strict_clock_when_no_stubs=self.strict_clock_when_no_stubs)

        for rec_nd in ann_tr.tree_reconstructed.preorder_node_iter():
            setattr(rec_nd, self.attr_name, rec_nd.mol_rate)

            # without this, 'self.attr_name' would be a plain Python
            # attribute, invisible to dendropy's Newick/Nexus writers
            # (e.g., 'AnnotatedTree.rec_str()', used by PhyloJunction's
            # GUI to display a tree node's string value) -- those only
            # serialize attributes a node has been explicitly told to
            # bind, the same way 'distribution.dn_discrete_sse' binds
            # "state"/"index" at simulation time. A branch with no
            # rate of its own (the root, or a sampled-ancestor branch;
            # see 'populate_gamma_spike_mol_rates_on_rec_tr_branches')
            # serializes as '[&...,mutation_rate=None]', which is
            # accurate, not an error.
            rec_nd.annotations.add_bound_attribute(self.attr_name)

    def init_check_vectorize_sample_size(
            self,
            param_list: ty.List[ty.Any] = []) -> None:
        """Validate Gamma-spike clock hyperparameters (currently scalar-only).

        FOR LATER: once per-sample vectorized clock hyperparameters
        are supported (e.g., a different 'clock_rate' per tree
        sample, drawn from an upstream DAG node), this method should
        grow each scalar hyperparameter into a list of length
        'self.n_samples', mirroring how distributions in
        'distribution.dn_parametric' do it today. For now, every
        hyperparameter is a plain scalar float (enforced by this
        class' type hints), so there is nothing to vectorize.
        """

        pass

    def get_rev_inference_spec_info(self) -> ty.List[str]:
        """See 'pgm.ValueGenerator.get_rev_inference_spec_info'.

        Not yet implemented: RevBayes has no known equivalent of the
        Gamma-spike clock model (it is a 2025 BEAST2-only addition at
        the time of writing), and PhyloJunction's own RevBayes
        inference-script export ('pjcli --inference-output') is
        itself still an unfinished feature (see 'pj_cli.py'). Revisit
        this once that feature is implemented for other distributions.
        """

        raise NotImplementedError(
            "RevBayes inference-script export is not yet implemented "
            "for GammaSpikeClockModel.")


if __name__ == "__main__":
    pass
