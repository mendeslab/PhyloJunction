.. role:: underline
    :class: underline

Phylogenetic tree stochastic nodes can be added to a model's DAG using tree distributions.
|pj| provides discrete SSE and continuous-trait QuaSSE tree distributions.

Discrete SSE
------------

Discrete SSE (``discrete_sse``) stands for discrete state-dependent speciation and extinction process, a general model that includes many popular special cases documented below.

:underline:`Rate types`

The implementation and usage of ``discrete_sse`` in |pj| is very flexible and takes inspiration from Tim Vaughan's great `ReMASTER <https://tgvaughan.github.io/remaster/>`_ simulator program: a model with any number of rates of whatever types may be specified.
Of course, this flexibility comes at the price of user responsibility!

SSE rates are added to a DAG by passing deterministic function ``sse_rate`` a stochastic or constant node containing a rate's value (the value can also be passed directly):

.. code-block:: 
    :caption: **Script excerpt 1.** Command in *phylojunction* to add a birth rate (deterministic node) to the DAG.

    birth_rate ~ exponential(rate=1.0)

    det_birth_rate := sse_rate(name="lambda", value=birth_rate, event="speciation")
    
    # alternatively, one could pass values directly to sse_rate()
    # det_birth_rate := sse_rate(name="lambda", value=[0.95, 1.05], event="speciation")

Parameters for ``sse_rate`` include:

    | ``name`` (string, required): The name of the rate parameter. Not particularly critical other than for printing sensical messages to the screen, when using the GUI. The argument to ``name`` does not have to match anything else in the model, and can in practice be any string enclosed in double quotes.
    | ``value`` (positive real, required): The rate value. This parameter can also take the name of other nodes in the DAG (in which case the value of those nodes are extracted and used).
    | ``event`` (string, required): A string specifying the rate type (see below).
    | ``epoch`` (positive integer, optional): Integer specifying the epoch to which this rate applies. See more in the "Time-heterogeneity" section below. Defaults to 1.

    | Arguments for the ``event`` parameter should be enclosed in double quotes and can be:

        i. ``"speciation"`` or ``"w_speciation"`` for within-state birth rates;
        ii. ``"b_speciation"`` for between-state birth rates;
        iii. ``"asym_speciation"`` for asymmetric birth rates;
        iv. ``"extinction"`` for extinction rates;
        v. ``"transition"`` for anagenetic transition rates;
        vi. ``"anc_sampling"`` for ancestor sampling (e.g., "fossilization") rates.

:underline:`Sampling probabilities`

The discrete SSE model in |pj| supports not only state-dependent rates, but also state-dependent sampling probabilities.
Specifying SSE probabilities can be done like so:

.. code-block:: 
    :caption: **Script excerpt 2.** Command in the *phylojunction* language to add a sampling probability parameter (deterministic node) to the DAG.

    det_sampling_prob := sse_prob(name="rho", value=1.0, state=0)
    # by setting the value to 1.0, we assume complete sampling!

The parameters for ``sse_prob`` are:

    | ``name`` (string, required): The name of the probability parameter. Behaves the same as the name of SSE rate parameters.
    | ``value`` (real between 0.0 and 1.0, required): The probability value. As with SSE rates, this parameter can also take the name of DAG nodes that hold probability values.
    | ``state`` (positive integer starting at 0, optional): This integer represents the state to which the probability is associated. Defaults to 0.
    | ``epoch`` (positive integer, optional): Integer specifying the epoch to which this probability applies. See more in the "Time-heterogeneity" section below. Defaults to 1.

Note that the ``state`` parameter defaults to 0, and is the state value used in non-state-dependent tree models, like the Yule or birth-death models.

:underline:`The SSE stash (and time-heterogeneity)`

Before ``discrete_sse`` can be assigned as the tree distribution to a DAG stochastic node, state-dependent rate and probability parameters must be stashed together and combined with some more information. 
This is accomplished with deterministic function ``sse_stash``, whose parameters are:

    | ``flat_rate_mat`` (string vector, required): A vector (deterministic) DAG node names holding the objects created with ``sse_rate``. DAG node names do not have to appear in any specific order.
    | ``flat_prob_mat`` (string vector, required): A vector (deterministic) DAG node names holding the objects created with ``sse_prob``. DAG node names do not have to appear in any specific order. 
    | ``n_states`` (positive integer, optional): How many states the SSE process has. Defaults to 1.
    | ``n_epochs`` (positive integer, optional): How many epochs characterize the skyline of the SSE process. Defaults to 1.
    | ``seed_age`` (positive real, optional): The age of the node at which the SSE process starts (either the origin or root node). This parameter must receive an argument if parameters are time-heterogeneous, and it must match the value given to ``stop_value`` when specifying ``discrete_sse`` with ``stop="age"`` (see more below).
    | ``epoch_age_ends`` (positive real vector, optional): The age of the ends of each epoch (i.e., the "younger" limit of each epoch, measured as an age), going from the oldest to the youngest epoch. The end of the youngest epoch (which contains the present) should be ommitted as it is added by default. The number of provided age ends should thus be the number of epochs (in ``n_epochs``) minus one. Defaults to [0.0].

The discrete SSE model in |pj| is a "skyline" model, that is, it allows both state-dependent rates and parameters to vary over time (i.e., to be time-heterogeneous).
In the presence of time-heterogeneity, these parameters are assumed to vary according to a user-defined piecewise constant function (this is coded into the model via the ``value`` and ``epoch`` parameters of ``sse_rate`` and ``sse_prob``).
The limits of each "piece" of that function are given by the argument to ``epoch_age_ends`` (see above).

.. code-block:: 
    :caption: **Script excerpt 3.** Command in the *phylojunction* language to organize rate and probability parameters into a stash (a deterministic node in the DAG).

    det_birth_rate1 := sse_rate(name="lambda1", value=1.0, event="speciation", epoch=1)
    det_birth_rate2 := sse_rate(name="lambda2", value=2.0, event="speciation", epoch=2)

    stash := sse_stash(flat_rate_mat=[det_birth_rate1, det_birth_rate2], n_epochs=2, seed_age=2.0, epoch_age_ends=[1.0])
    # there are two epochs, the oldest one ends at age 1.0, the youngest one ends at the present
    # the birth-rate of the oldest epoch is 1.0, but that rate doubles in the youngest epoch

:underline:`The discrete SSE distribution`

We have seen how to add to the DAG a series of ingredients required by the ``discrete_sse`` distribution.
The parameters of this function are:

    | ``stash`` (string, required): Name of DAG node holding the object created by ``sse_stash``.
    | ``start_state`` (nonnegative integer, optional): State at which each tree starts. Defaults to 0.
    | ``stop`` (string, required): String specifying the simulation stopping condition. If set to ``"age"``, the simulation will stop when the origin or root's age (see parameter ``origin``) reaches the value provided in ``stop_value``. If set to ``"size"``, the simulation will stop when the tree has the number of sampled nodes provided in ``stop_value``.
    | ``stop_value`` (positive numeric, required): The stopping condition threshold. If ``stop="age"``, this can be any positive real number. If ``stop="size"``, this can be any positive integer.
    | ``origin`` (string boolean, required): Flag specifying if process starts at origin node (as opposed to the root node). Defaults to "true".
    | ``cond_spn`` (string boolean, optional): Condition the process on at least one speciation event happening. Note that this speciation event may or not be the event represented by the root of the reconstructed tree. Defaults to "false".
    | ``cond_surv`` (string boolean, optional): Condition the process on having at least one surviving lineage at the present moment. This conditioning only makes sense when ``stop="age"``. Defaults to "true".
    | ``cond_obs_both_sides`` (string boolean, optional): Condition the process on having sampled nodes on both sides of the complete tree's root node. Defaults to "false".
    | ``min_rec_taxa`` (positive integer, optional): Minimum number of sampled (extant and sampled ancestors) nodes the tree should have. Defaults to 0.
    | ``max_rec_taxa`` (positive integer, optional): Maximum number of sampled (extant and sampled ancestors) nodes the tree should have. Defaults to 1e12.
    | ``abort_at_alive_count`` (positive integer, optional): The count of a tree's living nodes (irrespective of node sampling) at which point the tree is considered too large for simulation to continue. Trees can grow out of control when ``stop="age"`` and one of the rates allows it to grow too big before the stopping condition is met. This parameter can be used to throw away such trees and preventing PhyloJunction from crashing. Defaults to 1e12.
    | ``runtime_limit`` (positive integer, optional): Maximum number of seconds to obtain one accepted tree. Sampling is aborted at this point. Defaults to 300.
    | ``max_n_attempts`` (positive integer, optional): Maximum number of (failed) attempts before sampling is aborted. Defaults to 200.

Below, an example of assigning a discrete SSE to a stochastic DAG node representing phylogenetic tree samples:

.. code-block:: 
    :caption: **Script excerpt 4.** Commands in *phylojunction* to build a discrete SSE distribution and assign it to a stochastic DAG node.

    det_birth_rate1 := sse_rate(name="lambda1", value=1.0, event="speciation", epoch=1)
    det_birth_rate2 := sse_rate(name="lambda2", value=2.0, event="speciation", epoch=2)

    stash := sse_stash(flat_rate_mat=[det_birth_rate1, det_birth_rate2], n_epochs=2, seed_age=2.0, epoch_age_ends=[1.0])
    
    trs ~ discrete_sse(stash=stash, start_state=0, stop="age", stop_value=2.0, origin="true")

Each distribution draws one tree. Use ``iid(2, discrete_sse(...))`` for two trees with the
same parameters, and **Replicates** or ``--runs`` to repeat the complete model.
Vectorization and comprehensions allow different parameters for each tree.
Each tree has its own runtime and failed-attempt allowance, including retries.

:underline:`Special cases of discrete SSE tree models`

    1. `Yule <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/yule.pj>`_ (pure-birth) model;
    2. `Birth-death <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/birthdeath.pj>`_ model;
    3. `Birth-death with incomplete sampling <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/birthdeath_incomplete_sampling.pj>`_ model;
    4. `Fossilized birth-death <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/fbd.pj>`_ (FBD) model;
    5. `Binary state-dependent speciation and extinction <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/bisse.pj>`_ (BiSSE) model;
    6. `Time-heterogeneous binary state-dependent speciation and extinction <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/bisse_timehet.pj>`_ (skyline BiSSE) model;
    7. `Geographic state-dependent speciation and extinction <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/geosse.pj>`_ (GeoSSE) model;
    8. `Time-heterogeneous geographic state-dependent speciation and extinction <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/geosse_timehet.pj>`_ (skyline GeoSSE) model;

Reading trees and selecting elements
------------------------------------

``read_tree`` always returns a list, even when the input contains one tree. Select a scalar
tree explicitly with ``[1]``:

.. code-block:: text

    trees <- read_tree(string="((A:1,B:1):1,C:2);")
    first <- trees[1]

Reading and applying stochastic histories
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A history records the state changes along a tree's branches for one stochastic-map iteration.
Read the trees and histories, then attach each history to a copy of the corresponding tree:

.. code-block:: text

    trees := read_tree(file_path="examples/trees_maps_files/geosse_dummy_tree1.tre", node_name_attr="index")
    histories := read_maps(file_path="examples/trees_maps_files/geosse_dummy_tree1_maps.tsv", tip_attr_file_path="examples/trees_maps_files/geosse_dummy_tree1_tip_states.tsv", n_regions=2, geo=true)
    mapped := map_attr(tree=trees[1], history=histories, attr_name="state")

``read_maps`` returns one history per iteration identifier in first-appearance order.
``map_attr`` copies its tree before attaching the selected history. A scalar tree is shared
across a vector of histories; equally sized tree and history lists pair elementwise. Singleton
lists are not expanded: ``trees[1]`` selects the one tree explicitly. Both operations also work
inside comprehensions. History node identifiers must match the tree's ``index`` attributes;
``node_name_attr="index"`` reads those identifiers from the Newick annotations.

.. include:: quasse.rst

Reading complete trees
----------------------

``read_tree(..., is_complete=true)`` restores biological status from Newick annotations.
Each ordinary tip requires ``alive=true`` or ``alive=false``; dead tips additionally require
``sampled=true`` or ``sampled=false``. Sampled fossil tips survive reconstruction. Explicit
status resolves near-present fossils even when branch lengths are rounded; lengths are unchanged.
``keep_unsampled_alive=true`` treats explicitly unsampled living tips as sampled for reconstruction.

The default, ``is_complete=false``, reads an already reconstructed tree: every tip is sampled,
and living status is inferred from depth. Explicit living-status annotations or unsampled tips
are rejected in this mode. Both modes return a list of trees; select one with ``[1]``.

``n_rec_branches(tree)`` counts reconstructed edges with a parent and duration greater than
that tree's numerical tolerance. Sampled-ancestor marker edges are excluded. Edges on either
side of an observation remain separate. Branch-specific vectors follow preorder traversal.
Hidden-speciation counts (``n_stub``) belong to each reconstructed edge separately; observations
do not transfer counts from their incoming edge to their continuing edge.

Gamma-spike branch rates
------------------------

``gamma_spike_clock(tree, attr_name, clock_rate, spike_mean, spike_shape,
br_specific_rel_rate, strict_clock_when_no_stubs=false)`` describes a distribution on a copy
of one tree with rates on its reconstructed edges. Draw with ``~``. The input's times,
traits and stochastic maps are preserved. ``br_specific_rel_rate`` is a vector in the
preorder used by ``n_rec_branches``; generate it separately with ``iid``.

The gradual component is ``clock_rate * relative_rate``. Each hidden speciation on an edge
and each genuine parent speciation contributes a Gamma spike with mean ``spike_mean`` and
shape ``spike_shape``. Sampling observations contribute no spike. Dividing the total spike
size by edge duration gives its additional rate. ``spike_mean=0`` removes the spike effect.
All parameters must be finite; relative rates and spike_mean may be zero, whereas clock_rate
and spike_shape must be positive. The optional ``strict_clock_when_no_stubs`` forces stub-free
edges to the constant clock_rate; it is a diagnostic option, not the usual statistical model.

Displaying rates
~~~~~~~~~~~~~~~~

A tree with rate annotations initially displays the reconstructed tree in evolutionary-distance
units (time × rate). Uncheck **Evolutionary distance** to use time units and enable the
**Reconstructed** checkbox. Switching back to time restores its previous reconstruction choice.
Rates are currently assigned only to reconstructed edges, so complete-tree distance plots are
unavailable. The present-time reference line appears only in time views.

GUI summaries describe the selected complete or reconstructed topology. Ages remain in time
units in both views. Exported statistics and comparison/coverage calculations are independent
of the display selection. Python plotting also accepts ``branch_rate_attr`` to color edges by
a continuous rate attribute without replacing biological trait annotations.
