QuaSSE
------

QuaSSE couples a real-valued Brownian trait to trait-dependent speciation and extinction.
Simulation uses an approximate method based on ``diversitree``.

With PhyloJunction installed, run these commands from the repository root (where
``examples/`` is located):

.. code-block:: bash

   pjcli examples/quasse.pj -d -f 'trees;1-2' -r 42 -o /tmp/quasse-output
   pjcli examples/quasse_size.pj -d -r 42 -o /tmp/quasse-size
   pjcli examples/quasse_conditioned.pj -d -r 42 -o /tmp/quasse-conditioned


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

The trait increment over duration δ has mean ``drift × δ`` and variance ``diffusion × δ``.
``diffusion`` is variance per unit time, not a standard deviation. Drift and diffusion are
independent of the trait and constant through time within each simulation.
Zero diffusion is supported.
Each tree starts with one lineage at the origin. Daughters inherit the parent's trait and
then evolve independently. There are no fossil samples or sampled ancestors.

Rate-constructor parameters broadcast to a common vector length. Distribution model
parameters, rate objects, ``stop_value``, and ``sampling_prob`` have length one or ``n``.
Each of the ``n`` samples has ``nr`` independent replicates using the same parameter values;
output ordering is sample first, then replicate. Rejected trees do not resample DAG parents.

Diversitree method and stopping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At the current living traits, compute total birth–death rate R. A full step has duration
1/(kR), and one event occurs with probability 1/k. Choose the event type and lineage using
the current rates. Apply the event first, then extend and evolve the resulting living
lineages. A newborn daughter therefore evolves during that step, while a dead lineage does not.

For ``stop="age"``, shorten the final interval to the remaining duration δ and use event
probability Rδ. Surviving branches end exactly at the observation age. Completely extinct
trees retain the actual extinction times, with node ages measured relative to observation.

For ``stop="size"``, stop immediately before the birth that would increase the living count
from N to N+1. Omit that birth and its following trait update. The count can fall below N
before reaching this terminating birth. The stopping age is random. Complete sampling is
required. If survival conditioning is disabled, extinct outcomes can also be returned.

The interior algorithm follows ``diversitree/R/simulate-quasse.R``. Intentional differences are:

- The final age-limited interval and event probability are shortened, preventing overshoot.
- Size stopping omits diffusion after the terminating birth, matching PhyloJunction's convention.
- No artificial initial branch length of 1e-8 is added; genuine zero-length branches are valid.
- Tree storage, observation sampling, and conditioning use PhyloJunction's conventions.

The method is approximate. Increasing ``k`` generally improves resolution at greater cost;
``k`` does not specify an absolute time-step tolerance. It defaults to 500 and must be an integer
at least one. The default and currently only method is ``"diversitree"``.

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
     - Independent sampling probability for each living terminal
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
Survival does not imply observation: with ``sampling_prob=0``, a surviving tree has no observed
terminals. Root-side checks use the original root, before pruning can replace it.
For unconditional age-based simulation set ``cond_surv="false"`` and leave other conditions unset.
Explicit count bounds with size stopping, or incomplete sampling with size stopping, are errors.

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
     - Iterations per attempt
   * - ``max_alive``
     - 100,000
     - Living lineages per attempt

All limits are positive; count limits are integers. Reaching the rejection budget raises an
error. Numerical failures and exceeded resource limits also raise errors, without rejection
and retry. In particular, transiently exceeding ``max_rec_taxa`` does not reject an attempt.

When both rate functions are identically zero, age stopping uses a direct Brownian transition.
Size stopping then fails because its terminating birth is unreachable. Numerically zero rates
from nonzero functions raise an underflow error, rather than assuming rates stay zero forever.

Trees, plots, and output
^^^^^^^^^^^^^^^^^^^^^^^^

Complete trees retain extinct and unsampled tips. Reconstructed trees retain living sampled
paths. All complete-tree nodes carry the trait at that node's time; no full trait trajectory is
stored. Black branches and colored terminal markers show traits using a continuous colorbar.
The complete tree's terminal range determines colors in both complete and reconstructed views.
Plotting uses no randomness. An empty reconstructed tree displays ``No sampled tips``.

Data output includes complete/reconstructed Newick tables, annotated variants, generic tree
summaries, and ``<node>_traits.tsv``. The trait table has columns ``sample``, ``replicate``, ``node``,
``trait``, ``alive``, and ``sampled``; sample and replicate indices start at one. It includes internal
nodes and the origin as well as living and extinct tips. Discrete state tables and discrete
NEXUS character matrices are not generated for these trees.

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

Additional scripts
^^^^^^^^^^^^^^^^^^

Size stopping (`quasse_size.pj <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/quasse_size.pj>`_):

.. literalinclude:: ../../../examples/quasse_size.pj
   :language: text

Conditioning on observed tip count
(`quasse_conditioned.pj <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/quasse_conditioned.pj>`_):

.. literalinclude:: ../../../examples/quasse_conditioned.pj
   :language: text
