from __future__ import annotations

import math
import typing as ty
import numpy as np

# pj imports
import phylojunction.data.tree as pjtr

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
    and grammar machinery (the distribution wrapper handles those):
    it operates on a single 'AnnotatedTree' instance, in place, and
    can be imported and used on its own.

    Random draws are made with NumPy's legacy, global random number
    generator (module-level 'numpy.random.*' functions), the same
    convention used throughout PhyloJunction (e.g.,
    'distribution.dn_discrete_sse.DnSSE'). There is no 'seed'
    parameter here: reproducibility is controlled once, globally, by
    whichever caller already seeded that generator -- model execution's
    'random_seed' (set from a '.pj' script's '-r' CLI flag or the
    GUI's seed field).
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
        br_specific_rel_rates (list of float): One finite nonnegative relative rate per
            eligible reconstructed edge, in preorder. Conventionally these are centered
            at 1 (e.g. lognormal meanlog=-sdlog**2/2), so clock_rate is the mean absolute
            gradual rate. The model does not require that centering. Vector length must
            match ann_tr.n_rec_branches(); sampling the vector is the caller's responsibility.
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
        # draw. Its supplied relative-rate entry is ignored; no Gamma draw is made.
        if strict_clock_when_no_stubs and rec_nd.n_stub == 0:
            rec_nd.branch_raw_rate = 1.0
            rec_nd.spike_effect = 0.0
            rec_nd.mol_rate = clock_rate
            continue

        # branch-specific relative rate, externally provided (e.g.,
        # by an 'iid(n_rec_branches(tree), lognormal(...))' draw)
        rec_nd.branch_raw_rate = br_specific_rel_rates[this_branch_idx]

        # Number of Gamma-distributed spikes: one per hidden split on this segment,
        # plus one for the PARENT if it is a genuine speciation, not an SA marker.
        # Whether this node itself bifurcates is irrelevant to its incoming edge;
        # that bifurcation contributes to each of its daughters' edges instead.
        # This also covers a genuine root reached through a chain of SA observations.
        # Inspect reconstructed topology directly: suppression can leave lineage flags stale.
        # Douglas et al. (2025), eqn 2.7, attribute the extra spike to the parent;
        # likewise, getNSpikes() in GammaSpikeModel's BranchSpikePrior.java checks
        # node.getParent().isFake(), not the node's own children. An origin/unary
        # stem is not a speciation either (it can occur in singleton reconstructions).
        n_spikes = rec_nd.n_stub
        if not rec_nd.parent_node.is_sa_dummy_parent and rec_nd.parent_node.num_child_nodes() >= 2:
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
        )


