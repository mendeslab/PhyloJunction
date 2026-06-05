%%%%%%%%%%%%%
Documentation
%%%%%%%%%%%%%

.. contents::
    :local:
    :backlinks: none

--------
Overview
--------

Origins and design
==================

The |pj| project was born from our need to simulate arbitrarily complex SSE processes, but not being able to do so with many existing applications.
A series of related diversification models had been implemented in excellent computational methods, but each of those tools made unique modeling assumptions that we wanted to relax.

Instead of reverse-engineering and modifying code written by others, however, it seemed like the path of least resistance was to write our own simulator.
We also realized that our program could potentially come in handy in multiple projects, current and future, though only if written as more than a "one-off" chunk of code.

Ideally, |pj| would have to be "modular" in its architecture, both in terms of its codebase as well as how models should be represented and specified.
This modularity would allow future models to be more easily implemented -- potentially by other developers -- and seamlessly patched into |pj|.

With the above in mind, we reasoned our work would be more likely to bear fruits (and faster) if we used a programming language that:
    
    (i) was easy to prototype and debug model code in, and
    (ii) cross-platform,
    (iii) supported object-oriented programming, and for which
    (iv) biology and data-science libraries were available.

Python was our language of choice.
As for item (4), for example, we could make heavy use of the great `Dendropy <https://dendropy.org/>`_ library, maintained by Jeet Sukumaran and collaborators.

Modularity was achieved by building |pj| around a graphical modeling architecture.
This design (used in the past, see below) would not only make |pj| a more flexible simulator, but also more generally useful in the future, as a tool for method development and inference.

Graphical models
================

|pj| follows the graphical model paradigm popularized in the last decade by programs like `RevBayes <https://revbayes.github.io/>`__, `BEAST <https://beast.community/>`__ and `BEAST 2 <https://www.beast2.org/>`__.
In what follows, we shall offer just a brief refresher on what probabilistic graphical models are, and one simple example of their use in evolutionary modeling.
Those interested in a more detailed exposition will find it in Höhna et al. (2014) and references therein (also, take a look `here <https://revbayes.github.io/tutorials/intro/graph_models.html>`_).

Just as with the RevBayes and BEAST programs, specifying a model in |pj| amounts to writing down any number of conditional dependencies among random variables.
When looked at collectively, random variables comprise a type of probabilistic graphical model called a Bayesian network.
As the name suggest, such models can be described as graphs -- directed acyclic graphs (**DAG**), more specifically -- where each node (or vertex) represents a random variable, and each edge a conditional dependency between the two nodes (variables) it connects.

A Bayesian network representing a very simple model can be seen in figure 1a:

..  figure:: ../images/simple_graphical_model_both_manual.png
    :figwidth: 100%
    :scale: 40%
    :align: center

    **Figure 1.** A simple probabilistic model represented by (a) a Bayesian network, and (b) a factor graph.
    From the notation in the main text, :math:`\theta=\{sd, m\}`.
    The dashed box is called a "plate", and for each node it envelops it denotes multiple (in this case, precisely five) i.i.d. random variables.
    White and gray nodes denote random variables whose values we do (i.e., data) and do not know (i.e., parameters), respectively.

If we were to pen down what this Bayesian network represents in probability terms, we would arrive at a joint density (assuming continuous variables) given by expression:

.. math::
    :label: jointprob

    f_{\Theta,D}(\theta,d) = f_{\Theta|D}(\theta|D=d)f_D(d) = f_{D|\Theta}(d|\Theta=\theta)f_\Theta(\theta),

where :math:`\theta=\{sd, m\}`.

By comparing the Bayesian network graph and the expression above, the attentive reader will notice that functions are not represented in the graph.
We do not know exactly how random variables are related, other than :math:`d_1` through :math:`d_5` depending somehow on :math:`sd` and :math:`m`.
This is why it can be helpful to adopt a more general graphical representation: a factor graph (Fig. 1b).

The model underlying the factor graph in figure 1b is the same as that shown by the Bayesian network, the main difference being the additional factor nodes (the filled squares).
Factor nodes can make the relationships between two or more random variables more explicit, especially if their (factor) functions are annotated.

In the example above we can see a factor whose function :math:`f_{D|\Theta}(d|\Theta=\theta)` gives us the probability of :math:`\boldsymbol{d} = \{d_i: 1 \leq i \leq n\}` given :math:`\theta`.
It is annotated as "Normal", so we in fact know a close-form expression for this factor function; it is the probability density function of a normal distribution.
Random variables :math:`sd` and :math:`m` stand for the standard deviation and the mean of that normal distribution.
(For those curious to see just how complicated factor graphs can be, an example can be found in Zhang et al., 2023; see their Supplementary Fig. X.)

How to read a DAG?
------------------

Most biologists specifying a DAG on a computer do so as the first step of statistical **inference**.
They are interested in estimating unknown quantities about the natural world.
These quantities -- the random variables! -- whose values we do not know (but would like to estimate) are referred to as the **parameters** of the model.
Of course, parameter estimation requires that we know the value(s) of one or more random variables, which we naturally refer to as **data**.

In the example in figure 1, the parameters are :math:`sd` and :math:`m`, jointly referred to as :math:`\theta`; the data is node labeled :math:`d_i`.
(Note how they are colored differently, gray for data, white for parameter.)
Carrying out statistical inference thus implies reading the DAG from its data nodes, normally at the bottom, towards the parameter nodes above.

The most natural approach to parameter estimation requires that we simply take the middle and right-hand side terms in equation :eq:`jointprob`,
and solve for :math:`f_{\Theta|D}(\theta|D=d)`, the posterior density function:

.. math::
    :label: bayestheorem

    f_{\Theta|D}(\theta|D=d) = \frac{f_{D|\Theta}(d|\Theta=\theta)f_\Theta(\theta)}{f_D(d)}

This expression is known as Bayes theorem.
What programs like RevBayes, BEAST and BEAST 2 do is evaluate the posterior density function at several values of :math:`\theta`, and output the resulting (posterior) distribution.

----

|pj|, however, is among other things a collection of **simulation** (rather than estimation) engines.
Borrowing the jargon used above, we are primarily interested in generating values for random variables, some of which are parameters, some data.

Here, it makes sense to think of factors as the distributions from which values will be sampled (though as we will see below, factors can also represent deterministic functions).
As opposed to estimation, simulation starts at the top (or "outer") layers of the DAG, from parameters whose values we do know, and flows in the direction pointed to by arrows (normally downwards; Fig. 1b).

In figure 1b, for example, we start from known values for the parameters of the distributions at the top.
The exponential distribution from which we sample (i.e., simulate) :math:`sd` has a rate of 1.0; the standard normal distribution (:math:`Z`) from which we sample :math:`m`, by definition, has a mean of 0.0 and a standard deviation of 1.0.
We then define a normal distribution from the sampled values of :math:`sd` and :math:`m`, and in turn sample five times from that normal distribution, obtaining our data values in :math:`\boldsymbol{d}`.

Specifying a model (an example)
===============================

|pj| takes a model specification approach that sits between those adopted by the BEAST and RevBayes communities, and that largely intersects with the `LinguaPhylo <https://linguaphylo.github.io/>`_ project.
DAG-building instructions in |pj| are written in its own programming language, *phylojunction* (written in lowercase), whose syntax evokes RevBayes' Rev language.
But unlike Rev, *phylojunction* is not a fully fledged scripting language; it is lightweight and behaves more like a markup language such as XML, BEAST's format of choice.

Commands in *phylojunction* can be read as mathematical statements, and are naturally interpreted as instructions for building a node in a DAG.
Here is what a series of those commands would look like if written as a *phylojunction* script:

.. code-block:: 
    :caption: **Example script 1.** Script written in *phylojunction* specifying a time-homogenous birth-death model.

    # hyperprior
    m <- 0.0 # mean of log-normal below
    sd <- 0.1 # standard deviation of log-normal below
    
    # rate values
    d <- 1.0 # death
    b ~ lognormal(mean=m, sd=sd) # birth
    
    # deterministic rate containers
    dr := sse_rate(name="death_rate", value=d, event="extinction")
    br := sse_rate(name="birth_rate", value=b, event="speciation")
    
    O <- 2.0 # origin age
    
    # deterministic parameter stash
    s := sse_stash(flat_rate_mat=[dr, br], n_states=1, n_epochs=1) # parameter stash

    # phylogenetic tree
    T ~ discrete_sse(stash=s, stop="age", stop_value=O, origin="true")

Each line in the above script is a command instructing the application engine to take some form of user input, produce some value from it, and then store that value permanently in a new variable created on the spot.

Every command string consists of an assignment operator (e.g., ``<-``, ``:=``, ``~``) placed between the variable being created (on its left side) and some user input (on its right side).
We can look at some of these commands individually:

* ``d <- 1.0`` creates a variable ``d`` (the death rate) that is then passed and henceforth stores constant value ``1.0``. We use ``<-`` when we know the value of a random variable and want to create a constant node in the DAG;
* ``b ~ lognormal(mean=m, sd=sd)`` creates a variable ``b`` (the birth rate) that will store a random value drawn from a user-specified distribution. Here, that distribution is a log-normal with mean `m` and standard deviation `sd`.
* ``dr := sse_rate(name="death_rate", value=d, event="extinction")`` calls deterministic function ``sse_rate``, which creates variable ``dr`` from the value stored in ``d`` and some other user input. There is no stochasticity to deterministic assignments; they help make explicit steps where variables are transformed, annotated or combined with others.

And this is the DAG such script instantiates:

..  figure:: ../images/bd_graphical_model_manual.png
    :figwidth: 100%
    :scale: 40%
    :align: center

    **Figure 2.** Factor graph representing the time-homogenous birth-death model specified in example script 1.

Note how this DAG has a factor node characterized by a different type of function: deterministic function ``sse_rate``.
Such functions are denoted by filled diamonds instead of filled squares (those represent distributions), and will have their output enclosed in a hollow diamond.

Multiple samples and replicates
-------------------------------

Unlike other probabilistic programming languages, in *phylojunction* the output of every function is immutable and depends exclusively on that function's input.
Initialized variables cannot be altered by mutator methods or side-effects.

An important consequence of variable immutability is that it precludes loop control structures (e.g., *for* and *while* loops).
A reasonable question is then: *How does one have the model be sampled (i.e., simulated) multiple times?*

Let us look at the simple model shown in figure 1.
Ten independent samples of that model can be obtained with following *phylojunction* script:

.. code-block::
    :caption: **Example script 2.** Script written in *phylojunction* specifying the simple model in figure 1 and sampling (i.e., simulating) it ten times.

    sd ~ exponential(rate=1) # this value will be vectorized
    m ~ normal(mean=0.0, sd=1.0) # ... and so will this!
    d ~ normal(n=10, nr=5, mean=m, sd=sd)

As can be seen from the last line in script 2, all it took was providing argument ``n=10`` to the function call of a distribution.
But note how the first and second lines in script 2 do not set ``n=10``.
There is just a single value being drawn from exponential and standard normal distributions; implicitly, those commands are setting ``n=1``.

|pj| deals with this discrepancy in the requested number of samples by **vectorizing** the single values stored in :math:`sd` and :math:`m`.
For example, if the sampled value of :math:`m` is 0.1, under the hood |pj| converts :math:`m=1.0` into :math:`\boldsymbol{m}=\{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1\}`.
In such case, all ten samples requested by command ``d ~ normal(n=10, nr=5, mean=m, sd=sd)`` would thus come from normal distributions with the same mean of 0.1.

.. warning::
    |pj| will raise an error if the number of samples ``n`` specified in two commands are both different and greater than 1.
    In other words, vectorization can only be applied on variables holding a single value.
    (This type of behavior should be familiar to R users.)

----

In addition to simulating an entire model multiple times as explained above, it is also possible to specify a model with multiple i.i.d. random variables, alternatively referred to as **replicates**.
Multiple replicates can be specified by providing an additional argument ``nr`` to distribution calls, e.g., ``nr=5`` as in the last line of script 2.

In a DAG, i.i.d. random variables can be collectively represented by a single plated node (Fig. 1), or each be represented by an individual node.
There are in principle no constraints on the number of replicates a plated node may represent.
Errors will be thrown only if such nodes are used as input for a function, and the replicate count somehow violates that function's signature.

------------------------------
Graphical user interface (GUI)
------------------------------

.. _GUI:

|pj| comes with a (simple!) graphical user interface (GUI; see Fig. 3) that allows users to incrementally build a model while simultaneously inspecting any simulation output (among other things).

..  figure:: ../images/pjgui_model.png
    :figwidth: 100%
    :align: center

    **Figure 3.** |pj|'s GUI main window.

Like any modern computer application, |pj|'s GUI exposes its features to users via a menu (Fig. 3, number 2).
On the main tab ("Model"), one can build a model (DAG) with commands entered in the prompt (Fig. 3, number 1), or by loading a *phylojunction* script or a model from a previous |pj| section.
The latter can be done by clicking the application's menu (top-left corner of the Desktop on Apple machines), and clicking "File > Read script" or "File > Load model".
Model building commands can be revisited by going into the menu's "Command log" tab.

.. warning::
    If mistakes are made when entering a command, |pj| will flash a "Warning produced" message in red, at the bottom of the main window.
    The warning message can be seen by visiting the "Warnings" tab.

Given a model, DAG node values can be navigated as a plot, text string, or both (Fig. 3, numbers 3 and 4).
Users can cycle through replicated simulations (Fig. 3, number 5), and examine node-value summaries computed for individual simulations or across replicates (Fig. 3, number 6).
Node-value summaries include the mean and standard deviation for scalar variables, and statistics like the root age and number of tips for phylogenetic trees.

After a model is specified, saving the main plot (from the "Model" tab) to disk can be done by clicking "Save plot as".
Simulated data tables, or the model itself, in turn, can be saved to disk by clicking the application menu (top menu), and clicking "File > Save data as" or "File > Save model as".

Implementation comparison
=========================

The "Compare" tab (under "Model") streamlines the comparison of |pj|'s simulations with those from an external simulator, with respect to a summary statistic.
There are two main moments when carrying out this type of comparison may be useful or required.
The first is when a model is being first implemented in |pj|, and there happens to exist an independent implementation elsewhere.
This latter implementation can be used to provide a baseline expectation for |pj|.

Beyond model development within the platform, the "Compare" tab was further designed to expedite activities related to teaching.
Students in an evolutionary modeling course may be asked to code their own Yule model, for example, in whatever programming language.
After Yule trees are obtained, all they need to do is collect summary statistics that are also monitored by |pj| (e.g., the root age), and upload them into the "Compare" tab.
|pj| then automatically plots the two distributions (one from the student, one from |pj| itself) side by side, revealing if implementations behave similarly.

In order to use the functionalities exposed by the "Compare" tab, first one must to build a model in |pj|.
We can start with the simple model in `examples/multiple_scalar_tree_plated.pj <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/multiple_scalar_tree_plated.pj>`_.
As the name suggests, this model has a few different scalar and tree random variables, some of them plated.
We can load it by clicking "File > Read script" and selecting the aforementioned file.
By navigating to the "Compare" tab, you will see that |pj| has loaded the DAG nodes that could potentially be compared on the left-hand "Nodes to compare" menu.

In the "Compare" tab, now click "Compare to .csv (...)".
This button allows user to input a .csv file containing summary statistics computed from external simulations.
We will choose file `examples/compare_files/comparison_multiple_scalar_tree_plated_rv5.csv <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/compare_files/comparison_multiple_scalar_tree_plated_rv5.csv>`_.
Once the .csv file is loaded, |pj| should display the contents of it in the top window.

..  figure:: ../images/pjgui_compare.png
    :figwidth: 100%
    :align: center

    **Figure 4.** |pj|'s GUI "Compare" tab.

Note the "rv5" column name at header of the table -- these are the variable names that |pj| will try to look up (and compare with) in its own model's list of nodes.
Accordingly, select "rv5" from the top-left "Node to compare" menu.
Now we are ready to draw, click "Draw".
The produced graph (Fig. 4) should indicate that |pj| and the external simulator behave similarly.

Let us now try a different random variable, a tree this time: click on node "trs2" on the "Node to compare" menu.
Because tree space is complicated (there are both a discrete and a continuous dimensions to it, the topology and branch lengths, respectively), we will need to look at tree summary statistics.
Different tree statistics should have been listed in the bottom-left "Summary statistics" menu.
Choose one of them, "root age", for example.

Now we must then load a different .csv file, `examples/compare_files/comparison_multiple_scalar_tree_plated_trs2.csv <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/compare_files/comparison_multiple_scalar_tree_plated_trs2.csv>`_.
This file contains tree summary statistics presumably produced by a different external simulator.
Click "Draw".
Again, the produced graph suggests the model implemented in both simulators behave similarly.

Coverage validation
===================

The "Coverage" tab was developed for automating coverage validation, a procedure for verifying the correctness of a model implemented within a Bayesian framework.
Coverage validation in fact validates multiple things simultaneously: the generating (simulator) and inferential implementations of the model, and any machinery used in inference (e.g., MCMC proposals).
For the sake of clarity and conciseness, we will assume in this section that we are validating a model's inferential engine, and that all other involved code works as intended.

The gist of coverage validation is straightforward: the :math:`\alpha` %-HPD over a model parameter produced during Bayesian inference should contain the "true" (i.e., simulated) parameter value approximately :math:`\alpha` % of the time.
For example, if a stochastic node in the DAG is sampled 100 times (i.e., 100 i.i.d. random variables), then for approximately 95 of those simulations the 95%-HPD interval will contain the true sampled value.

In order to carry out coverage validation with the GUI, first we must build a model in |pj|.
Click "File > Read script" and select `examples/coverage_files/r_b_exp.pj <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/coverage_files/r_b_exp.pj>`_.
In the loaded model, one hundred independent samples for a parameter :math:`r_b` were assigned to a constant node ``r_b`` named after it.
These 100 values were drawn from an exponential distribution implemented elsewhere, but could also have been sampled directly in |pj|.

In the "Coverage" tab, one can now load a .csv file containing a table with :math:`r_b` 's posterior distribution's mean value, and the lower and upper bounds of an HPD interval.
(This posterior distribution came from a Bayesian analysis done with a different program.)
The table headers must be ``posterior_mean``, ``lower_95hpd``, and ``higher_95hpd``.
Click "Read HPDs .csv (...)", and select `examples/coverage_files/r_b_exp.csv <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/examples/coverage_files/r_b_exp.csv>`_.
The table should have been shown in the top middle window.

..  figure:: ../images/pjgui_coverage.png
    :figwidth: 100%
    :align: center

    **Figure 5.** |pj|'s GUI "Coverage" tab.

To inspect coverage, select "r_b" from the top-left menu "Non-det. nodes", and then click on "Draw".
The x- and y-axes in the plot (the main panel in Fig. 5) show to the true simulated values and the posterior means, respectively.
Vertical bars denote the 95%-HPDs over :math:`r_b`, with blue bars (95 of them) indicating intervals containing the true value, and red bars (5 of them) not containing it.
In other words, one can deduce from the graph that coverage for :math:`r_b` was 0.95 (though this is also shown in the top-right panel; Fig. 5); this coverage is adequate and indicative of a correctly implemented model.

In addition to loading a .csv file with pre-calculated posterior distribution summaries, users can also tell |pj| the path of a directory containing MCMC output .log files.
After loading the same model as before and clicking on "r_b" in the top-left menu, click on "Directory to .log's (...)", and choose `examples/coverage_files/logfiles/`.
If in the .log files the node being examined were to be named something different from "r_b", users can tell |pj| what to look for by entering the alternative name in "Parameter name in .log".
For example, one could type "r_b" in "Parameter name in .log" (in this case it is not necessary).

Now click on "Draw".
An almost identical plot should have been produced, and coverage is again appropriate (0.94).
(The difference of 1 in coverage has to do with the behavior of different libraries when computing HPD boundaries.)

----------------------------
Command-line interface (CLI)
----------------------------

.. _CLI:

|pj| has a command-line interface (CLI) as an entry point for some of its functionalities -- it is particularly handy when running |pj| is a step in a longer pipeline (perhaps being run on a cluster).
The easiest way to see what is possible with the CLI is to call it on the Terminal with the ``-h`` flag, like so:

.. code-block:: bash

    pjcli -h

The help message should list a series of flags, which we expand on below:

    | ``-d``: Simulated data tables are written to disk.
    | ``-f``: Plots for the selected DAG nodes (whose names and sample ranges are passed as arguments) are written to disk. The argument for ``-f`` should be passed within single-quotes, e.g., 'birth_rate,trs;1-2,1-10'.
    | ``-o``: Directory where to place simulated data tables and figures.
    | ``-p``: Prefix to prepend output files with (e.g., "simulation1")
    | ``-nex``: Certain simulated data is written to disk in NEXUS format.
    | ``-r``: Random seed so that simulations always produce the same output.

Note that if no argument is passed to ``-o``, data tables written to disk (prompted by ``-d``) will be placed in the directory from which PhyloJunction is called.
Figures (prompted by ``-f``), in turn, will be placed in a "figures/" sub-directory, next to the data tables. 

---------------------------------
Sandbox: Bypassing the interfaces
---------------------------------

.. _bypass:

Users and developers who wish to bypass |pj|'s standalone interfaces to explore the code base can do so by using ``pj_sandbox.py``.
Of course, some familiarity with Python is required.


The ``pj_sandbox.py`` `script <https://raw.githubusercontent.com/fkmendes/PhyloJunction/main/src/phylojunction/interface/pj_sandbox.py>`_ contains a few examples of how |pj| can be imported as a module and used for different purposes.
Each example is wrapped within its own function, like so:

.. code-block:: python
    :caption: **Example 1 in pj_sandbox.py.** The first example in the sandbox script demonstrates how one can build a simple DAG with a (Yule) tree stochastic node.

    # pj imports
    import phylojunction.interface.cmdbox.cmd_parse as cmdp
    import phylojunction.pgm.pgm as pgm
    ...

    def run_example_yule_string() -> pgm.DirectedAcyclicGraph:
        """Run example 1.
    
        This example sets up a simple Yule model with two samples,
        and two Yule tree replicates per sample.
        """

        yule_model_str = \
            ('n_sim <- 2\nn_rep <- 2\nbirth_rate <- 1.0\n'
             'det_birth_rate := sse_rate(name="lambda", value=birth_rate,'
             ' event="speciation")\n'
             'stash := sse_stash(flat_rate_mat=[det_birth_rate])\n'
             'trs ~ discrete_sse(n=n_sim, nr=n_rep, stash=stash, start_state=[0],'
             ' stop="size", stop_value=10.0, origin="false")\n')
    
        dag_obj = cmdp.script2dag(yule_model_str, in_pj_file=False)

        return dag_obj

The example above is then invoked later in ``pj_sandbox.py``:

    .. code-block:: python
        :caption: **Calling example 1 in pj_sandbox.py.** Which example to execute can be specified in the script main body.

        dag_obj: pgm.DirectedAcyclicGraph = pgm.DirectedAcyclicGraph()

        ...

        # examples
        # 1: Yule model string, 2 tree samples, 2 tree replicates per sample
        # 2: Same as (1), but reading a pre-made .pj script in examples/yule.pj
        # 3: BiSSE model with incomplete sample, 2 tree samples, 2 tree replicates per sample
        # 4: Builds discrete SSE tree manually, then prints on screen
        # 5: Builds discrete SSE tree from newick string, then prints on screen
        # 6: Read .pj script examples/see_stoch_maps.pj

        example_to_run = 1
        # example_to_run = 2
        # example_to_run = 3
        # example_to_run = 4
        # example_to_run = 5
        # example_to_run = 6
            
        if example_to_run == 1:
            dag_obj = run_example_yule_string()

        ...


-------
Lexicon
-------

Parametric distributions
========================

.. include:: parametric.rst

Tree distributions
==================

.. include::  treedists.rst