In addition to plain constant values (e.g., ``birth_rate <- 1.0``), the ``<-`` operator can also call a small number of "named" constant functions, which build something more involved than a single literal value from the arguments given to them.

.. _read-tree:

Reading a tree from Newick
---------------------------

Constant function ``read_tree`` builds a tree from a Newick string (perhaps obtained from a different simulator, or from an empirical study).
After a tree is read in with this function, it can be annotated (e.g., with stochastic mappings using ``map_attr``) or transformed (e.g., scaled by a clock model).

.. code-block::
    :caption: **Example command.** Reading a tree from a Newick string, provided directly inline.

    trs <- read_tree(n=1, string="(nd1:2.0,nd2:2.0):0.0;")

The parameters for ``read_tree`` are:

    | ``n`` (integer, optional): Number of samples in the provided file (see ``file_path``below). Behaves in one of exactly two ways: if neither ``n`` nor ``nr`` is given at all, every tree string found (see ``file_path`` below) goes straight into ``n`` (never into ``nr``) with no further input needed; if ``n`` and/or ``nr`` *is* given (even as ``n=1``), ``n`` :math:`\times` ``nr`` must equal how many tree strings are provided in ``file_path``. Defaults to 1.
    | ``nr`` (integer, optional): Number of replicates per sample. See ``n`` above for the same two possible behaviors that will also apply to ``nr``. Defaults to 1.
    | ``file_path`` (string, required unless ``string`` is given): Path to a file containing one or more Newick tree strings, one per line. Mutually exclusive with ``string`` (which only ever holds a single tree).
    | ``string`` (string, required unless ``file_path`` is given): A single Newick tree string, provided directly. Mutually exclusive with ``file_path``.
    | ``node_name_attr`` (string, optional): Name of an NHX annotation (e.g., ``[&index=9]``) to relabel every node in the tree by. This is needed to cross-reference the tree's nodes against a user-provided annotation, such as a stochastic map (see ``map_attr`` on the deterministic functions page), by whatever node indexing the user employed.
    | ``epoch_age_ends`` (positive real vector, optional): The ages of epoch ends, oldest to youngest, omitting the youngest (which is added by default). This parameter is needed for plotting epoch boundaries, but for nothing else. Also, see ``sse_stash``'s parameter of the same name (see the tree distributions page).
    | ``is_complete`` (string boolean, optional): Whether the Newick string represents a *complete* tree (both sampled and unsampled lineages), rather than an already-reconstructed tree (|pj|'s default assumption, under which every terminal node is assumed sampled). When ``"true"``, every leaf **must** carry explicit ``&alive=`` or ``&sampled=`` annotations (e.g., ``nd1[&alive=True,sampled=True]:2.0``), since branch lengths and depth alone cannot otherwise distinguish a fossil tip (dead, but sampled) from an unsampled extinction. Defaults to ``"false"``. When left at its default, neither annotation is ever consulted -- ``alive``/``sampled`` are deduced purely from depth/construction instead -- so an explicit ``&alive=`` of either value, or an explicit ``&sampled=False``, on any leaf is rejected outright (an error is raised) rather than silently ignored, since either is most likely a sign that ``is_complete="true"`` should have been passed instead.
    | ``keep_unsampled_alive`` (string boolean, optional): Only meaningful when ``is_complete="true"``. A terminal node that is alive may or not have been sampled, so if ``&sampled=False``, that terminal node is removed by default. Users may on occasion, however, keep a terminal node that they know exist but was not sampled -- ``keep_unsampled_alive`` makes it so unsampled but living terminal nodes are not ignored (as an extinct or unsampled node would normally be) in any manipulation involving the reconstructed tree associated with the (complete) tree provided to ``read_tree()``. Defaults to ``"false"``.

The example below shows how to read several trees written in a file specified by ``file_path``. Because neither ``n`` nor ``nr`` are specified, all trees are read into ``n``:

.. code-block::
    :caption: **Example command.** A two-line file read with neither 'n' nor 'nr' given: both trees become two samples (not replicates) of 'trs'.

    trs <- read_tree(file_path="examples/trees/two_FBD_trees.tre")
    # 'two_FBD_trees.tre' holds two lines, each a Newick tree string,
    # so 'trs' now has sample_size=2, repl_size=1 -- as if 'n=2' had
    # been written explicitly

Giving ``n`` and/or ``nr`` explicitly, on the other hand, allows one to build replication (per sample) into the model design:

.. code-block::
    :caption: **Example command.** The same file, but with explicit 'n=1, nr=2' specifying two replicates for one sample.

    trs <- read_tree(n=1, nr=2, file_path="examples/trees/two_FBD_trees.tre")

Below, an example of reading in a complete tree from an annotated Newick string. The tree here is a 5-taxon fossilized birth-death (FBF) tree.
It has one direct (sampled) ancestor ('sa1'), one fossil tip ('nd6', dead but sampled), one unsampled extinct lineage ('nd2'), and two living taxa ('nd3', 'nd5'):

.. code-block::
    :caption: **Example command.** Reading a complete FBD tree whose terminal nodes have been annotated for '&alive='/'&sampled=' as required.

    trs <- read_tree(n=1, string="((sa1[&alive=False,sampled=True]:0.0,((nd3[&alive=True,sampled=True]:0.7008060330779379,(nd5[&alive=True,sampled=True]:0.43018486984623605,nd6[&alive=False,sampled=True]:0.10676020530326771)nd4_nd4:0.27062116323170177)nd1_nd1:0.3777323878021377,nd2[&alive=False,sampled=False]:0.15074253577213667)root_root:0.41924114926716083)brosc_dummy1:0.5022204298527637)origin_origin:0.0;", is_complete="true")
