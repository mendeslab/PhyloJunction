QuaSSE
------

QuaSSE couples a real-valued Brownian trait to trait-dependent speciation, extinction and optional fossil sampling.
Two simulation methods are available: the default ``diversitree`` approximation and
``local_thinning``, which evaluates rates at Brownian-bridge candidate traits and has an
explicit probability-error budget under the assumptions below.

With PhyloJunction installed, run these commands from the repository root (where
``examples/`` is located):

.. code-block:: bash

   pjcli examples/quasse.pj -d -f 'trees;1-2' -r 42 -o /tmp/quasse-output
   pjcli examples/quasse_size.pj -d -r 42 -o /tmp/quasse-size
   pjcli examples/quasse_conditioned.pj -d -r 42 -o /tmp/quasse-conditioned
   pjcli examples/quasse_local_thinning.pj -d -r 42 -o /tmp/quasse-local


Each ``.pj`` statement must occupy one line. ``-d`` writes data; ``-f`` selects figures.
The GUI's normal tree display also colors continuous terminal traits.

Model
^^^^^

.. literalinclude:: ../../../examples/quasse.pj
   :language: text

Download `quasse.pj <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/quasse.pj>`_.

``quasse_logistic`` matches Java SSE's link function:

.. code-block:: text

   f(x) = y0 + (y1 − y0) / (1 + exp(slope × (midpoint − x)))


The plateau values are finite and nonnegative. The midpoint and slope are finite. Either
slope sign and either plateau ordering are supported. Equal plateaus produce a rate independent
of the trait. Rate functions have no explicit time dependence; rates change as traits evolve.

Other rate families are available for speciation, extinction or fossil sampling:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Constructor and required arguments
     - Rate at trait x
   * - ``quasse_constant(rate)``
     - rate
   * - ``quasse_gaussian(baseline, center_rate, center, width)``
     - baseline + (center_rate − baseline) exp(−z²/2), z = (x − center)/width
   * - ``quasse_step(left, right, threshold)``
     - left when x < threshold; right otherwise, including equality
   * - ``quasse_linear(intercept, slope)``
     - max(0, intercept + slope × x)
   * - ``quasse_skew_gaussian(baseline, amplitude, location, width, skew)``
     - baseline + 2 amplitude exp(−z²/2) Φ(skew × z), z = (x − location)/width
   * - ``quasse_quadratic(baseline, strength, center)``
     - baseline + strength × (x − center)²
   * - ``quasse_absolute(baseline, strength, center)``
     - baseline + strength × abs(x − center)

Supply arguments by name in scripts. All parameters must be finite; rate levels, amplitudes
and strengths must be nonnegative, and widths positive. Locations, slopes, skew and the linear
intercept can have either sign. Gaussian curves include troughs when ``center_rate < baseline``.
Φ is the standard normal cumulative distribution function. The asymmetric peak's ``location``
and ``amplitude`` are not generally its peak position and height; ``skew=0`` gives a symmetric
Gaussian peak. All constructor parameters are constant through time within a simulation.

Every constructor accepts optional ``cap``: the evaluated rate is min(cap, f(x)). The cap must
be finite and nonnegative, may be below the baseline, and may be zero. It broadcasts with the
other parameters. Without it, linear, quadratic and absolute-value rates remain unbounded
when their slope or strength is nonzero. There is no hidden upper truncation.

An explicit cap changes the biological model; computational resource limits do not cap rates.
For example, pure branching Brownian motion with quadratic speciation can have an infinite
expected population after a finite time even though each finite-time population is finite
almost surely. A large realization or exceeded node budget is not proof of infinitely many
nodes. General birth–death models should not be assigned an explosion classification from
these runtime checks.

The trait increment over duration δ has mean ``drift × δ`` and variance ``diffusion × δ``.
``diffusion`` is variance per unit time, not a standard deviation. Drift and diffusion are
independent of the trait and constant through time within each simulation.
Zero diffusion is supported.
Each tree starts with one lineage at the origin. Daughters inherit the parent's trait and
then evolve independently. Optional ``fossil_rate`` supplies a rate object for non-destructive
sampling along each lineage, in observations per lineage per unit time. It defaults to zero
(``None`` in Python), accepts every rate family, and broadcasts at length one or ``n``.
Fossils record the exact trait at sampling; they do not remove the lineage or reset its trait.
Sampling at death, removal probabilities and measurement error are not supported.

Rate-constructor parameters broadcast to a common vector length. Distribution model
parameters, rate objects, ``stop_value``, and ``sampling_prob`` have length one or ``n``.
Each of the ``n`` samples has ``nr`` independent replicates using the same parameter values;
output ordering is sample first, then replicate. Rejected trees do not resample DAG parents.

Diversitree method and stopping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At the current living traits, compute total birth–death–sampling rate R. A full step has duration
1/(kR), and one event occurs with probability 1/k. Choose the event type and lineage using
the current rates. Apply the event first, then extend and evolve the resulting living
lineages. A newborn daughter therefore evolves during that step, while a dead lineage does not.
A fossil records the start-of-step trait and its continuing lineage evolves for the whole step.

For ``stop="age"``, shorten the final interval to the remaining duration δ and use event
probability Rδ. Surviving branches end exactly at the observation age. Completely extinct
trees retain the actual extinction times, with node ages measured relative to observation.

For ``stop="size"``, N counts **sampled living tips**. The initial lineage receives a
Bernoulli(``sampling_prob``) mark. At birth, one randomly chosen daughter inherits its parent's
mark and the other receives an independent fresh mark. Thus a birth adds zero or one to
the sampled count. A sampled death subtracts one; fossils preserve the continuing mark.

Stop immediately before a proposed birth would increase that count from N to N+1, omitting
that birth and its following trait update. Merely reaching N does not stop the simulation:
zero-increment births can occur at N and sampled deaths can lower the count again. Do not
resample marks at the end. ``0 < sampling_prob <= 1`` is required for size stopping; zero
is permitted for age stopping. All living lineages, including unsampled ones, evolve and
count toward ``max_alive``. Zero sampled tips with living lineages is not extinction.
If survival conditioning is disabled, extinct outcomes can also be returned.

On a fixed biological tree independent of the marks, these leaf marks are independent
Bernoulli draws. The flag-dependent stopping time conditions their distribution; the
returned size-stopped tree is defined by this marked process, not independent post-sampling
of a tree stopped by its total living count. At complete sampling the old behavior is preserved.

The interior algorithm follows ``diversitree/R/simulate-quasse.R``. Intentional differences are:

- The final age-limited interval and event probability are shortened, preventing overshoot.
- Size stopping omits diffusion after the terminating birth, matching PhyloJunction's convention.
- No artificial initial branch length of 1e-8 is added; genuine zero-length branches are valid.
- Tree storage, observation sampling, and conditioning use PhyloJunction's conventions.

The method is approximate. Increasing ``k`` generally improves resolution at greater cost;
``k`` does not specify an absolute time-step tolerance. It defaults to 500 and must be an integer
at least one. The default method is ``"diversitree"``.

Optional ``dt_max`` is a positive finite duration. When supplied, each interval is at most
min(1/(kR), dt_max, remaining age), with event probability R times that interval. This can
resolve trait changes where present rates are small, but is not an error tolerance or an exact
simulation method. Step discontinuities and rapidly varying or unbounded rates require
resolution checks using smaller ``dt_max`` and larger ``k``. Omitting it preserves the original
step selection for positive total rates.

Brownian bridge local thinning
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Select ``method="local_thinning"`` and supply a positive finite ``block_duration`` in the
model's time units. Each living lineage samples a Brownian endpoint up to that duration away.
The endpoints and a Brownian-bridge tail bound define a likely containing trait interval.
The sum of bounds on birth, death and fossil rates over that interval bounds a Poisson proposal process.
At each candidate time the simulator samples the actual intermediate trait from the stored
bridge and selects birth, death, fossil sampling or rejection using actual rates.
Rejected candidates and fossil events retain the same endpoint and proposal bound. Sampling
only replaces the living node by its continuing child; it does not allocate another bridge
or error allowance. At birth the parent ends at that time and its daughters evolve independently.
There is no one-event-per-step restriction and no requirement that rates change slowly.

Longer blocks save endpoint work but can enlarge proposal bounds and increase rejected
candidates. Thus ``block_duration`` primarily controls cost. It need not be extremely small
for accuracy. Local thinning does not use ``dt_max`` or ``k``: a supplied ``dt_max`` or
nondefault ``k`` is rejected; the inherited default ``k=500`` is ignored.

``bridge_error`` defaults to ``1e-8`` and must be finite and strictly between zero and one.
It budgets excursions outside the likely trait regions across a complete ``generate()``
call, including every sample, replicate, biological rejection and numerical restart.
Block j receives probability allowance δ/[j(j+1)], whose infinite sum is δ.
A direct ``simulate()`` call gets its own budget. Separate calls start new budgets.
All existing rate families are supported, including uncapped linear, quadratic and
absolute-value functions: a finite local bound suffices, without a global rate cap.

A finite evaluated total rate exceeding its proposal bound automatically restarts the
current tree with the same model parameters. RNG progression, error budget, deadline and
remaining work allowance persist. Completed earlier trees are retained. This does not
consume a biological rejection attempt. A loose valid bound merely wastes candidate work;
an invalid bound can miss events, so restarting does not make this an exact sampler.

With exact arithmetic and sampling, valid interval bounds, and a terminating, nonexplosive
reference process, the total variation distance from the exact complete output law is at
most ``bridge_error``. The coupling agrees until the first region excursion, which also
covers detected violations and their restarts. This is an absolute probability bound,
not a relative bound on population means or other unbounded statistics. It is not a claim
of second-order convergence in block duration. Floating-point error and computational
interruptions are outside this guarantee. Biological rejection-limit failure is included
as an output; discarding failures and conditioning on success needs a separate bound
(at most 2δ/p, capped at one, where p is exact success probability and both are positive).

Age stopping ends living branches exactly at the requested age. Size stopping retains the
same N-to-N+1 convention as diversitree, but uses the continuous candidate birth time and
advances every unfinished bridge to that time. The stored future endpoints are retained
when sampling final traits. Both methods use the same inherited marks for size stopping and
independent end-of-simulation sampling for age stopping. The probability bound includes the
marked stopping law: before a region excursion, shared mark draws give identical stopping
decisions. Adding fossil intensity to the same region's bound needs no extra error allowance.

The derivation and assumptions are in ``TODO/simulator.tex`` in the source repository.

Observation and conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Argument
     - Default
     - Meaning
   * - ``sampling_prob``
     - 1
     - Independent terminal sampling for age stopping; inherited Bernoulli marks for size stopping
   * - ``cond_surv``
     - ``"true"``
     - Require at least one living lineage
   * - ``cond_spn``
     - ``"false"``
     - Require a retained speciation event
   * - ``cond_obs_both_sides``
     - ``"false"``
     - Require living sampled descendants on both sides of the original root
   * - ``min_rec_taxa``
     - 0
     - Minimum living sampled-tip count, for age stopping
   * - ``max_rec_taxa``
     - unbounded
     - Maximum living sampled-tip count, for age stopping

Equal count bounds request an exact observed count. Sampling precedes these checks.
Survival does not imply observation: with ``sampling_prob=0``, a surviving tree has no sampled living
tips, but may retain fossil observations. Root-side checks use the original root, before pruning can replace it.
For unconditional age-based simulation set ``cond_surv="false"`` and leave other conditions unset.
Explicit reconstructed-taxon count bounds with size stopping are errors. Fossils do not
satisfy living-tip count, root-side or survival conditions; sampling alone is not speciation.

Rejected attempts use fresh simulation randomness with unchanged model parameters. Rare or
impossible conditions can exhaust the rejection budget. Limits are:

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Argument
     - Default
     - Scope
   * - ``max_n_attempts``
     - 200
     - Total rejected attempts across one ``generate()`` call
   * - ``runtime_limit``
     - 300 seconds
     - Entire ``generate()`` call, including retries
   * - ``max_steps``
     - 1,000,000
     - Diversitree steps or local queue actions per attempt, shared across numerical restarts
   * - ``max_alive``
     - 100,000
     - Living lineages per attempt
   * - ``max_nodes``
     - unset
     - All nodes per attempt, including extinct nodes, fossils, sampling dummy nodes, origin and initial lineage

All limits are positive; count limits are integers. Reaching the rejection budget raises an
error. Except for the finite local bound-exceedance recovery described above, numerical
failures and exceeded resource limits raise errors without retry. In particular, transiently
exceeding ``max_rec_taxa`` does not reject an attempt.

``max_nodes`` is checked before creating the initial two nodes or the two new nodes at a birth or fossil event.
The omitted size-stopping birth does not consume nodes. This is a memory safeguard, not an
infinite-population detector; hitting it raises an error without returning a truncated tree.

When all three rate functions are identically zero, age stopping uses a direct Brownian
transition, even with ``dt_max`` set. Identically zero birth and death rates make size stopping
unreachable even if fossil sampling is positive.
For diversitree, if the total rate is zero only at the current traits (including underflow), supply
``dt_max`` so traits can advance and rates be reevaluated; without it, simulation raises an error.
Local thinning instead renews Brownian blocks even when a local proposal bound is zero.

Trees, plots, and output
^^^^^^^^^^^^^^^^^^^^^^^^

Complete trees retain extinct and unsampled tips. Reconstructed trees retain paths to sampled
living tips and fossils, including fossils with extinct or unobserved descendants and trees
with fossils alone. A lone fossil retains its original age and origin-to-observation stem. All complete-tree nodes carry the trait at that node's time; no full trait trajectory is
stored. Black branches and colored terminal markers show traits using a continuous colorbar.
Ancestral fossils appear as trait-colored markers along their continuing branches; terminal
fossils have their own rows. ``sa_along_branches=False`` gives every fossil its own row.
The complete tree's terminal range, including fossils, determines colors in both views.
In complete-tree plots, paths from retained observations back to the full-tree root are drawn
normally; other edges are lighter (25% opacity) and thinner (65% linewidth). This includes any
displayed ancestral stem above the reconstructed root. Fossil markers, tip colors and labels
remain unchanged, and reconstructed plots retain their usual styling.
Plotting uses no randomness. An empty reconstructed tree displays ``No sampled tips``.

Data output includes complete/reconstructed Newick tables, annotated variants, generic tree
summaries, and ``<node>_traits.tsv``. The trait table has columns ``sample``, ``replicate``, ``node``,
``trait``, ``alive``, and ``sampled``; sample and replicate indices start at one. It includes internal
nodes and the origin as well as living and extinct tips. Discrete state tables and discrete
NEXUS character matrices are not generated for these trees.

``<node>_fossils.tsv`` contains ``sample``, ``replicate``, ``node``, ``time``, ``age``, ``trait``.
Indices are one-based, ``time`` is forward from the origin, and ``age`` is relative to the
simulation horizon. There is one row per fossil in the complete tree, or a header-only table
when none occur. Annotated Newick identifies fossils with ``is_sa`` and sampling nodes with
``is_sa_dummy_parent``. The generic ``Direct ancestor count`` counts all fossil observations,
including those made terminal by reconstruction. ``Total taxon count`` retains its existing
biological-terminal meaning.

QuaSSE RevBayes inference export is not supported and fails before script generation.

Python interface
^^^^^^^^^^^^^^^^

.. code-block:: python

   from phylojunction.calculation.continuous_sse import LogisticRate
   from phylojunction.distribution.dn_quasse import DnQuaSSE
   
   simulator = DnQuaSSE(
       birth_rate=LogisticRate(0.4, 0.8, 0.0, 2.0),
       death_rate=LogisticRate(0.02, 0.1, 1.0, 2.0),
       stop="size", stop_value=10, diffusion=0.1, rng_seed=42,
   )
   trees = simulator.generate()


``generate()`` applies conditions and returns a list. ``simulate(sample_idx=0)`` returns one
attempt after observation sampling, without rejection conditioning. The Python constructor's
``rng_seed`` and CLI ``-r`` use NumPy's existing global random stream; instances do not own
independent random generators.

The Python rate classes are ``ConstantRate``, ``LogisticRate``, ``GaussianRate``, ``StepRate``,
``LinearRate``, ``SkewGaussianRate``, ``QuadraticRate`` and ``AbsoluteRate``, all in
``phylojunction.calculation.continuous_sse``. They accept the same named arguments as the
script constructors and evaluate NumPy arrays or scalars. Their shared ``RateFunction``
interface exposes ``is_identically_zero`` and ``upper_bound`` (a conservative finite bound,
or ``None`` when unavailable). A zero rate at one trait is not an identically zero function.
``bound_on_interval(left, right)`` returns a finite conservative bound on a closed finite
interval, including zero-width intervals, or raises ``ValueError`` if unavailable.
Positive bounds have an outward floating-point guard; this is not certified interval arithmetic.

Additional scripts
^^^^^^^^^^^^^^^^^^

Size stopping (`quasse_size.pj <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/quasse_size.pj>`_):

.. literalinclude:: ../../../examples/quasse_size.pj
   :language: text

Conditioning on observed tip count
(`quasse_conditioned.pj <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/quasse_conditioned.pj>`_):

.. literalinclude:: ../../../examples/quasse_conditioned.pj
   :language: text

Bounded and explicitly unbounded rates
(`quasse_rates.pj <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/quasse_rates.pj>`_):

.. literalinclude:: ../../../examples/quasse_rates.pj
   :language: text

Local thinning with uncapped quadratic speciation
(`quasse_local_thinning.pj <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/quasse_local_thinning.pj>`_):

.. literalinclude:: ../../../examples/quasse_local_thinning.pj
   :language: text

Fossil sampling and incomplete size sampling
(`quasse_fossils.pj <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/quasse_fossils.pj>`_):

.. literalinclude:: ../../../examples/quasse_fossils.pj
   :language: text

The three calls produce mixed observations, fossil-only observations, and size-stopped trees
with exactly five sampled living tips. The last call can retain additional unsampled living
tips and any number of fossils; neither contributes to its target count.
