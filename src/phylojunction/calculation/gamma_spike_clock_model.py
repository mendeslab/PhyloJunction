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


def populate_gamma_spike_mol_rates_on_complete_tree_branches(
        ann_tr: pjtr.AnnotatedTree,
        clock_rate: float,
        spike_mean: float,
        spike_shape: float,
        br_specific_rel_rates_complete: ty.List[float]) -> None:
    """Populate molecular-rate attrs on complete-tree branches too.

    'populate_gamma_spike_mol_rates_on_rec_tr_branches' (this module)
    annotates the *reconstructed* tree only, treating every stub
    (unsampled bifurcation) hidden along a reconstructed branch as a
    single, lumped contribution. This function spreads that same
    information out over the *complete* tree's own, finer-grained
    branches, and also rates every branch that has no reconstructed
    counterpart at all (one belonging to an extinct lineage, or to a
    lineage still alive today but never sampled):

    (i)  A reconstructed branch with no hidden stubs maps exactly
         onto one complete-tree branch; that one keeps the
         reconstructed branch's own, already-computed values,
         unchanged.
    (ii) A reconstructed branch with one or more hidden stubs maps
         onto several complete-tree branches in a row (its
         "backbone"). Each one keeps the reconstructed branch's own
         'branch_raw_rate' (a rate is constant along the whole
         branch, by this model's own definition, so there is nothing
         to divide there), but the reconstructed branch's single,
         already-drawn 'spike_effect' total is divided across them --
         using a Dirichlet distribution, since a sum of 'n' i.i.d.
         Gamma(alpha, theta) draws is, conditional on that sum, split
         into proportions distributed exactly as Dirichlet(alpha,
         ..., alpha) (with 'n' repeats of 'alpha'), independently of
         the sum itself. This gives each backbone segment a share
         that is statistically indistinguishable from one drawn
         directly, rather than an arbitrary split that merely adds up
         to the right total.
    (iii) A branch with no reconstructed counterpart at all gets a
          fresh, independent draw for both components: the gradual,
          relaxed-clock component from 'br_specific_rel_rates_complete'
          (externally supplied, exactly like 'br_specific_rel_rates' is
          for reconstructed branches), and the spike component from
          this same Gamma mechanism, with exactly one spike (or zero,
          if the branch is an 'is_sa_lineage' one) -- a complete-tree
          branch is, by construction, never itself subdivided by
          further hidden stubs, so there is no "+1 baseline, plus one
          per stub" sum to make here; it is just the ordinary single-
          branch rule with zero stubs.

    This function must be called after
    'populate_gamma_spike_mol_rates_on_rec_tr_branches' has already
    populated 'ann_tr.tree_reconstructed' -- it reads that tree's own
    'branch_raw_rate', 'spike_effect', and 'n_stub' values rather than
    recomputing them. Every node of 'ann_tr.tree' ends up with
    'branch_raw_rate', 'spike_effect', and 'mol_rate' populated (set
    to 'None' on the root and any zero-length branch, same as the
    reconstructed-tree function does). This also leaves the stem
    leading to the reconstructed tree's own seed node without a rate
    -- that node's own label is present in the reconstructed tree
    (it is not "absent" by this function's own definition), but the
    reconstructed tree's seed node has no incoming edge of its own to
    inherit a rate from either, so there is nothing to copy or split
    there.

    See that function's own docstring for every argument except
    'br_specific_rel_rates_complete', which plays the same role as its
    'br_specific_rel_rates', just for the branches absent from the
    reconstructed tree -- it must have exactly
    'ann_tr.n_complete_only_branches()' entries, one per branch, in
    the same (preorder) order 'AnnotatedTree.complete_only_branch_nodes'
    returns them in.

    Args:
        ann_tr (AnnotatedTree): Tree instance whose complete tree's
            branches will be annotated, in place.
        clock_rate (float): See
            'populate_gamma_spike_mol_rates_on_rec_tr_branches'.
        spike_mean (float): See
            'populate_gamma_spike_mol_rates_on_rec_tr_branches'.
        spike_shape (float): See
            'populate_gamma_spike_mol_rates_on_rec_tr_branches'.
        br_specific_rel_rates_complete (list of float): Pre-generated
            per-branch relative rates for the branches absent from
            the reconstructed tree. Must have exactly
            'ann_tr.n_complete_only_branches()' entries.

    Raises:
        ValueError: Is raised if 'br_specific_rel_rates_complete' is
            given with the wrong number of entries.
    """

    complete_only_nodes = ann_tr.complete_only_branch_nodes()

    if len(br_specific_rel_rates_complete) != len(complete_only_nodes):
        raise ValueError(
            "'br_specific_rel_rates_complete' had " +
            str(len(br_specific_rel_rates_complete)) +
            " entries, but the complete tree has " +
            str(len(complete_only_nodes)) +
            " branches absent from the reconstructed tree -- they "
            "must match.")

    for nd in ann_tr.tree.preorder_node_iter():
        nd.branch_raw_rate = None
        nd.spike_effect = None
        nd.mol_rate = None

    # (iii) branches with no reconstructed counterpart at all
    for idx, nd in enumerate(complete_only_nodes):
        nd.branch_raw_rate = br_specific_rel_rates_complete[idx]

        if nd.is_sa_lineage:
            nd.spike_effect = 0.0

        else:
            nd.spike_effect = spike_mean * np.random.gamma(
                shape=spike_shape, scale=1.0 / spike_shape)

        nd.mol_rate = (
            clock_rate * nd.branch_raw_rate
            + nd.spike_effect / nd.edge_length)

    # (i)/(ii) branches that do have a reconstructed counterpart
    for rec_nd in ann_tr.tree_reconstructed.preorder_node_iter():
        if not pjtr.is_real_branch(rec_nd, ann_tr.epsilon):
            continue

        complete_nd = ann_tr.tree.find_node_with_label(rec_nd.label)
        complete_parent = ann_tr.tree.find_node_with_label(
            rec_nd.parent_node.label)

        # the hidden, stub-bearing nodes between 'complete_parent' and
        # 'complete_nd', oldest first, mirroring the walk in
        # 'AnnotatedTree.populate_n_stubs_on_branches_of_rec_tr' -- a
        # dummy node (the parent of a sampled ancestor) counts here
        # exactly like any other collapsed node; it gets no exemption
        # (see that method's own docstring for why)
        stub_segment_nodes = []
        current_nd = complete_nd.parent_node

        while current_nd != complete_parent:
            stub_segment_nodes.append(current_nd)
            current_nd = current_nd.parent_node

        stub_segment_nodes.reverse()

        # 'complete_nd' itself is the segment carrying the branch's
        # own "+1 baseline" -- included last, unconditionally; see
        # below for why it may end up with no Dirichlet share of its
        # own despite being part of the backbone
        backbone_segment_nodes = stub_segment_nodes + [complete_nd]

        n_spikes = (
            rec_nd.n_stub if rec_nd.is_sa_lineage else rec_nd.n_stub + 1)

        # when 'rec_nd.is_sa_lineage', the "+1 baseline" is dropped,
        # so only the 'n_stub' hidden segments (never 'complete_nd'
        # itself, which sits last) are conceptually due a share of
        # the spike total. A zero-length one among them (e.g. a
        # sampled ancestor whose own dummy parent was collapsed into
        # it, leaving it with no surviving sibling lineage) cannot
        # itself carry a rate -- there is no time for one to apply
        # over -- so it is dropped here, and the *real*-length
        # segments alone share the whole total instead. This is not
        # an approximation: the Dirichlet distribution's own
        # aggregation property means the real segments' own
        # proportions, conditional on summing to the whole total (as
        # they must, since a zero-length segment can only ever
        # contribute exactly 0), are still Dirichlet-distributed with
        # one entry removed per zero-length segment dropped.
        candidate_sharing_segments = backbone_segment_nodes[:n_spikes]
        real_sharing_segments = [
            seg_nd for seg_nd in candidate_sharing_segments
            if pjtr.is_real_branch(seg_nd, ann_tr.epsilon)]

        n_real = len(real_sharing_segments)

        if n_real == 0:
            shares: ty.List[float] = []

        elif n_real == 1:
            # a single segment needs no splitting -- it gets the
            # reconstructed branch's own total, unchanged (this is
            # also exactly the no-stub case, (i) above)
            shares = [rec_nd.spike_effect]

        else:
            shares = (
                np.random.dirichlet([spike_shape] * n_real)
                * rec_nd.spike_effect
            ).tolist()

        for seg_nd, share in zip(real_sharing_segments, shares):
            seg_nd.branch_raw_rate = rec_nd.branch_raw_rate
            seg_nd.spike_effect = share
            seg_nd.mol_rate = (
                clock_rate * seg_nd.branch_raw_rate
                + share / seg_nd.edge_length)

        # any segment not among 'real_sharing_segments' -- whether
        # beyond the first 'n_spikes' (only ever 'complete_nd'
        # itself, and only when 'rec_nd.is_sa_lineage'), or a zero-
        # length one dropped above -- still needs its own rate if it
        # is itself a real branch (it carries no spike of its own);
        # a zero-length one is left with no rate at all, same as
        # 'is_real_branch' excludes it everywhere else in this module
        for seg_nd in backbone_segment_nodes:
            if seg_nd in real_sharing_segments:
                continue

            if not pjtr.is_real_branch(seg_nd, ann_tr.epsilon):
                continue

            seg_nd.branch_raw_rate = rec_nd.branch_raw_rate
            seg_nd.spike_effect = 0.0
            seg_nd.mol_rate = clock_rate * seg_nd.branch_raw_rate


