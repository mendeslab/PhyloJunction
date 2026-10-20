|pj| also ships with a number of deterministic functions (assigned with the ``:=`` operator), which introduce no stochasticity of their own -- they transform, combine, or annotate values already stored elsewhere in the DAG.
On this page you will find the deterministic functions involved in annotating a tree and applying a relaxed-clock model to it.
Every function below takes an already-built tree as its ``tree`` argument -- normally one read in with the :ref:`read-tree` constant function, on the constant functions page.

Mapping external annotations onto a tree
------------------------------------------

Deterministic function ``map_attr`` takes a tree (read in with :ref:`read-tree` and its ``node_name_attr`` argument) and "paints" it with a discrete-state history computed by an outside tool -- a stochastic character map, for example -- whose mid-branch (anagenetic) and cladogenetic state changes are given in a separate maps file.
The tree itself is left untouched topologically; only its ``at_dict``/``clado_at_dict``/``node_attr_dict`` members (the per-node attribute bookkeeping used when plotting, or when combined with a clock model, see ``gamma_spike_clock`` below) are populated.

.. code-block::
    :caption: **Example command.** Reading a tree by node index, then mapping a GeoSSE stochastic map onto it.

    trs <- read_tree(n=1, file_path="examples/trees_maps_files/geosse_dummy_tree2.tre", node_name_attr="index")
    mapped_trs := map_attr(tree=trs, fun="smap", maps_file_path="examples/trees_maps_files/geosse_dummy_tree2_maps.tsv", tip_attr_file_path="examples/trees_maps_files/geosse_dummy_tree2_tip_states.tsv", n_regions=2, geo="true")

The parameters for ``map_attr`` are:

    | ``tree`` (string, required): Name of the DAG node holding the tree(s) to be mapped (normally a :ref:`read-tree` node, with ``node_name_attr`` set so nodes can be cross-referenced by the index the external tool used).
    | ``fun`` (string, required): The kind of mapping being applied. Currently, the only supported value is ``"smap"`` (a stochastic character map).
    | ``maps_file_path`` (string, required): Path to a tab-separated file listing every anagenetic and cladogenetic state change along the tree (one row per change), by node index.
    | ``tip_attr_file_path`` (string, required): Path to a tab-separated file listing each tip's own (node-index, state) pair.
    | ``n_regions`` (positive integer, optional): Number of discrete areas/regions, if the state space is a GeoSSE-style combination of presence/absence across regions. Defaults to 1.
    | ``geo`` (string boolean, optional): Whether the state space should be interpreted as GeoSSE-style (presence/absence across ``n_regions`` regions) rather than a plain, flat set of discrete states. Defaults to ``"false"``.

.. warning::
    The resulting mapped attribute is currently always named ``"state"`` (hard-coded), regardless of any other argument -- there is no way yet to map onto, or name, a different per-node attribute.

Counting observable branches
----------------------------

Deterministic function ``n_rec_branches`` counts a tree's observable branches (every branch of its *reconstructed* topology, excluding the root and any zero-length branch -- see ``is_real_branch`` in the code).
The output of this function is a simple tree summary statistic, but as a DAG node, it does carry information about the number of samples and replicates of the trees that got summarized. This information can then be used in later model specification steps (e.g., ``n_nr_from`` argument of parametric distributions such as ``lognormal``; see Gamma-spike example below).

.. code-block::
    :caption: **Example command.** Counting a tree's real branches, then drawing one log-normal value per branch.

    n_rec_br := n_rec_branches(tree=trs)
    br_specific_rel_rate ~ lognormal(n_nr_from=n_rec_br, meanlog=-0.125, sdlog=0.5, log_space="true")

The only parameter for ``n_rec_branches`` is:

    | ``tree`` (string, required): Name of the DAG node holding the tree(s) whose branches are to be counted.

Gamma-spike relaxed clock
--------------------------

Deterministic function ``gamma_spike_clock`` implements the Gamma-spike punctuated relaxed-clock model of Douglas et al. (2025), annotating a tree's *reconstructed* branches with a new, continuous rate attribute (e.g., ``"mutation_rate"``) that combines:

    (i) branch-specific "raw" rates from the relaxed-clock component, supplied via ``br_specific_rel_rate``, and
    (ii) a punctuated component, the "spike effect", from hidden speciation events ("stubs") on
         each reconstructed edge and the speciation event at its parent, if any.

A sampled ancestor divides a lineage into two edge segments. Each segment keeps its own stub
count; sampling contributes no spike and does not transfer stubs between segments. A branch
leaving a genuine speciation receives one additional spike contribution; a branch leaving a
sampled-ancestor observation receives none.

The clock model reads a tree's topology and branch lengths (the tree should be a chronogram with branch lengths measured in absolute time). Any pre-existing discrete-state annotations (e.g., from ``map_attr``, above) are left untouched on the returned, independently-copied tree.

.. code-block::
    :caption: **Example command.** Drawing a relaxed-clock rate first (via 'n_rec_branches'/'lognormal'), then applying the Gamma-spike clock model.

    n_rec_br := n_rec_branches(tree=trs)
    br_specific_rel_rate ~ lognormal(n_nr_from=n_rec_br, meanlog=-0.125, sdlog=0.5, log_space="true")
    phylo := gamma_spike_clock(tree=trs, attr_name="mutation_rate", clock_rate=0.001, spike_mean=0.01, spike_shape=2.0, br_specific_rel_rate=br_specific_rel_rate)

The parameters for ``gamma_spike_clock`` are:

    | ``tree`` (string, required): Name of the DAG node holding the tree(s) to annotate with clock-model rates.
    | ``attr_name`` (string, required): Name of the new per-branch attribute to populate (e.g., ``"mutation_rate"``).
    | ``clock_rate`` (positive real, required): The mean absolute substitution rate (to scale the branch-specific raw rates).
    | ``br_specific_rel_rate`` (string, required): Name of a DAG node holding branch-specific raw rates (one raw rate per branch; conventionally centered at 1.0 in real space, e.g. from ``lognormal(n_nr_from=...)``, above) -- this is how the relaxed-clock component's own distribution and hyperparameters are chosen, independently, in the script.
    | ``spike_mean`` (nonnegative real, required): Mean of the Gamma-distributed punctuated "spike" component.
    | ``spike_shape`` (positive real, required): Shape of the Gamma-distributed punctuated "spike" component.
    | ``strict_clock_when_no_stubs`` (string boolean, optional): When ``"true"``, a branch with no accumulated stubs is pinned to a fixed, strict-clock rate. This is **not** a realistic model of molecular evolution -- it exists only as a debugging/diagnostic device (e.g., to make the rarer, stub-affected branches stand out when coloring a tree by ``attr_name``). Defaults to ``"false"``.

All numeric clock parameters and branch-specific rates must be finite. Branch-specific rates
must be nonnegative (zero is allowed); ``clock_rate`` and ``spike_shape`` must be strictly positive.

Below, a complete worked example: two FBD trees in Newick format stored in a file to be read as one sample with two replicates (``n=1, nr=2`` -- see :ref:`read-tree`), with the Gamma-spike clock model applied to each, as in `examples/fbd_read_gamma_spike.pj <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/fbd_read_gamma_spike.pj>`_.
A version of this same pipeline applied to a tree additionally "painted" with a stochastic map (via ``map_attr``) can be found in `examples/sse_smap_gamma_spike_clock.pj <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/sse_smap_gamma_spike_clock.pj>`_.
