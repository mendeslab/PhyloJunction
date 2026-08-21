import dendropy as dp
import matplotlib

# pj imports
import phylojunction.interface.cmdbox.cmd_parse as cmdp
import phylojunction.pgm.pgm as pgm
import phylojunction.data.tree as pjtr
import phylojunction.data.sampled_ancestor as pjsa
import phylojunction.readwrite.pj_read as pjr
import phylojunction.data.attribute_transition as pjat
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.distribution.dn_discrete_sse as distsse

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"

def run_example_yule_string() -> pgm.DirectedAcyclicGraph:
    """Run example 1.
    
    This example sets up a simple Yule model with two samples,
    and two Yule tree replicates per sample.
    """

    yule_model_str = \
        ('n_sim <- 2\nn_rep <- 2\nbirth_rate <- 1.0\n'
         'start_states <- [0, 0, 0, 0]\n'''
         'det_birth_rate := sse_rate(name="lambda", value=birth_rate,'
         ' event="speciation")\n'
         'stash := sse_stash(flat_rate_mat=[det_birth_rate])\n'
         'trs ~ discrete_sse(n=n_sim, nr=n_rep, stash=stash, start_state=start_states,'
         ' stop="size", stop_value=10.0, origin="false")\n')
    
    dag_obj = cmdp.script2dag(yule_model_str, in_pj_file=False)

    return dag_obj


def run_example_geosse_file() -> pgm.DirectedAcyclicGraph:
    """Run example 2.
    
    This example reads script examples/bisse.pj. The tree
    is a sample from a geographic state-dependent speciation
    extinction (GeoSSE) process.
    """

    dag_obj = cmdp.script2dag("examples/geosse.pj", \
                              in_pj_file=True)
    
    return dag_obj


def run_example_manual_incomplete_sampling_bisse() -> pgm.DirectedAcyclicGraph:
    """Run example 3.

    This example manually builds a BiSSE distribution with incomplete
    sampling, assigns it to a stochastic DAG node, and adds the node
    to a DAG.
    """

    total_n_states = 2

    #############################
    # SSE rates with two states #
    #############################
    rates_s0 = [
        sseobj.DiscreteStateDependentRate(
            name="lambda0",
            val=1.0,
            event=sseobj.MacroevolEvent.W_SPECIATION,
            states=[0,0,0]),
        sseobj.DiscreteStateDependentRate(
            name="mu0",
            val=0.5,
            event=sseobj.MacroevolEvent.EXTINCTION,
            states=[0]),
        sseobj.DiscreteStateDependentRate(
            name="q01",
            val=0.25,
            event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION,
            states=[0,1])
    ]

    rates_s1 = [
        sseobj.DiscreteStateDependentRate(
            name="lambda1",
            val=1.0,
            event=sseobj.MacroevolEvent.W_SPECIATION,
            states=[1,1,1]),
        sseobj.DiscreteStateDependentRate(
            name="mu1",
            val=0.5,
            event=sseobj.MacroevolEvent.EXTINCTION,
            states=[1]),
        sseobj.DiscreteStateDependentRate(
            name="q10",
            val=0.25,
            event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION,
            states=[1,0])
    ]

    rates = rates_s0 + rates_s1
    matrix_atomic_rates = [rates]
    state_dep_rate_manager = \
        sseobj.DiscreteStateDependentParameterManager(
            matrix_atomic_rates,
            total_n_states)
    event_handler = \
        sseobj.MacroevolEventHandler(state_dep_rate_manager)

    ##################################
    # Sampling probability parameter #
    ##################################    
    sampling_probs = [
        sseobj.DiscreteStateDependentProbability(
            name="rho",
            val=0.25,
            state=0),
        sseobj.DiscreteStateDependentProbability(
            name="rho", 
            val=0.75,
            state=1)
    ]
    
    matrix_state_dep_probs = [sampling_probs]
    state_dep_prob_manager = \
        sseobj.DiscreteStateDependentParameterManager(
            matrix_state_dep_probs,
            total_n_states)
    prob_handler = \
        sseobj.DiscreteStateDependentProbabilityHandler(
            state_dep_prob_manager
        )

    ######################
    # building SSE stash #
    ######################
    sse_stash = sseobj.SSEStash(event_handler, prob_handler)

    ###################################
    # Discrete SSE distribution setup #
    ###################################

    stop_condition = "age"
    stop_condition_value = [2.75]  # 2.75 time units
    start_at_origin = True
    sample_size = 2
    start_states_list = [0, 1]

    dist_obj = distsse.DnSSE(
        sse_stash,
        n=sample_size,
        origin=start_at_origin,
        start_states_list=start_states_list,
        stop=stop_condition,
        stop_value=stop_condition_value,
        epsilon=1e-8,
        runtime_limit=3000,
        debug=False)
    
    #####################
    # DAG-related setup #
    #####################
    
    dag_obj = pgm.DirectedAcyclicGraph()

    a_stoch_node_name = "trs"

    stoch_node_dag = pgm.StochasticNodeDAG(
        a_stoch_node_name,
        sample_size=sample_size,
        replicate_size=1,
        sampled_from=dist_obj,
        parent_nodes=None,
        clamped=False)
    
    dag_obj.add_node(stoch_node_dag)

    return dag_obj


def run_example_manual_tree_building(ax: matplotlib.pyplot.Axes) -> None:
    """Run example 4.

    This example manually builds an AnnotatedTree instance manually
    and plots it on provided Axes object.
    """

    def build_tree() -> pjtr.AnnotatedTree:
        origin_node = dp.Node(taxon=dp.Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node.state = 0
        origin_node.annotations.add_bound_attribute("state")
        origin_node.index = 0
        origin_node.annotations.add_bound_attribute("index")
        origin_node.alive = False
        origin_node.sampled = False
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        dummy_node = dp.Node(taxon=dp.Taxon(label="dummy1"), label="dummy1", edge_length=1.0)
        dummy_node.state = 0
        dummy_node.annotations.add_bound_attribute("state")
        dummy_node.index = 1
        dummy_node.annotations.add_bound_attribute("index")
        dummy_node.alive = False
        dummy_node.sampled = False
        dummy_node.is_sa = False
        dummy_node.is_sa_dummy_parent = True
        dummy_node.is_sa_lineage = False

        origin_node.add_child(dummy_node)

        # right child of dummy_node
        sa_node = dp.Node(taxon=dp.Taxon(label="sa1"), label="sa1", edge_length=0.0)
        sa_node.state = 0
        sa_node.annotations.add_bound_attribute("state")
        sa_node.index = 2
        sa_node.annotations.add_bound_attribute("index")
        sa_node.alive = False
        sa_node.sampled = True
        sa_node.is_sa = True
        sa_node.is_sa_dummy_parent = False
        sa_node.is_sa_lineage = False

        # left child of dummy node
        root_node = dp.Node(taxon=dp.Taxon(label="root"), label="root", edge_length=0.5)
        root_node.state = 1
        root_node.annotations.add_bound_attribute("state")
        root_node.index = 3
        root_node.annotations.add_bound_attribute("index")
        root_node.alive = False
        root_node.sampled = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = True

        dummy_node.add_child(sa_node)
        dummy_node.add_child(root_node)

        # left child of root node
        extant_sp1 = dp.Node(taxon=dp.Taxon(label="sp1"), label="sp1", edge_length=0.25)
        extant_sp1.state = 2
        extant_sp1.annotations.add_bound_attribute("state")
        extant_sp1.index = 4
        extant_sp1.annotations.add_bound_attribute("index")
        extant_sp1.alive = False
        extant_sp1.sampled = False
        extant_sp1.is_sa = False
        extant_sp1.is_sa_dummy_parent = False
        extant_sp1.is_sa_lineage = False

        # right child of root node
        extant_sp2 = dp.Node(taxon=dp.Taxon(label="sp2"), label="sp2", edge_length=0.5)
        extant_sp2.state = 3
        extant_sp2.annotations.add_bound_attribute("state")
        extant_sp2.index = 5
        extant_sp2.annotations.add_bound_attribute("index")
        extant_sp2.alive = True
        extant_sp2.sampled = True
        extant_sp2.is_sa = False
        extant_sp2.is_sa_dummy_parent = False
        extant_sp2.is_sa_lineage = False

        root_node.add_child(extant_sp1)
        root_node.add_child(extant_sp2)

        tr_sa_with_root_survives = dp.Tree(seed_node=origin_node)
        tr_sa_with_root_survives.taxon_namespace.add_taxon(origin_node.taxon)
        tr_sa_with_root_survives.taxon_namespace.add_taxon(sa_node.taxon)
        tr_sa_with_root_survives.taxon_namespace.add_taxon(dummy_node.taxon)
        tr_sa_with_root_survives.taxon_namespace.add_taxon(root_node.taxon)
        tr_sa_with_root_survives.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_sa_with_root_survives.taxon_namespace.add_taxon(extant_sp2.taxon)

        total_state_count = 4

        sa_global_time = 1.0
        time_to_sa_lineage_node = 0.5
        sa = pjsa.SampledAncestor(
            "sa1",
            "root",
            sa_global_time,
            state=0,
            time_to_lineage_node=time_to_sa_lineage_node)
        sa_lineage_dict = {"root": [sa]}

        at1 = pjat.AttributeTransition("state", "root", 1.25, 0, 1)
        at2 = pjat.AttributeTransition("state", "sp1", 1.6, 1, 2)
        at3 = pjat.AttributeTransition("state", "sp2", 1.8, 1, 3)
        at_dict = {
            "root": [at1],
            "sp1": [at2],
            "sp2": [at3]
        }

        max_age = 2.0

        ann_tr_sa_with_root_survives_max_age = \
            pjtr.AnnotatedTree(
                tr_sa_with_root_survives,
                total_state_count,
                start_at_origin=True,
                max_age=max_age,
                sa_lineage_dict=sa_lineage_dict,
                at_dict=at_dict,
                epsilon=1e-12
            )
        
        return ann_tr_sa_with_root_survives_max_age

    ann_tr_sa_with_root_survives_max_age = build_tree()

    print(ann_tr_sa_with_root_survives_max_age.tree.as_string(schema="newick"))

    # need this to paint branches
    attrs_of_interest = ["state"]
    ann_tr_sa_with_root_survives_max_age.\
        populate_nd_attr_dict(attrs_of_interest)

    pjtr.plot_ann_tree(ann_tr_sa_with_root_survives_max_age,
                       ax,
                       use_age=False,
                       start_at_origin=True,
                       sa_along_branches=True,
                       attr_of_interest="state")
                       # sa_along_branches=False


def run_example_read_tree_function(ax: matplotlib.pyplot.Axes) -> None:
    """Run example 5.

    This example reads a newick string into an AnnotatedTree instance, and
    plots it in the provided Axes object.
    """

    # NOTE: at the moment, the terminal-node states must be passed in the newick string
    # for colors to show correctly in plot; later implement deterministic function
    # that maps attributes to terminal nodes (at_dict would still be built by hand
    # in example)
    
    # tr_str = "(((nd10:0.05585634266420793,(nd12:0.01444225652017525,nd13:0.01444225652017525)nd11:0.04141408614403268)nd2:1.225238932333245,((nd6:0.15152443528855317,nd7:0.15152443528855317)nd4:0.16841576897032884,(nd8:0.11388809491478649,nd9:0.0982030908340934)nd5:0.2060521093440955)nd3:0.9611550707385709)root:0.7189047250025469)origin:0.0;"
    tr_str = "(((nd10[&state=1]:0.05585634266420793,(nd12[&state=0]:0.01444225652017525,nd13[&state=2]:0.01444225652017525)nd11[&state=2]:0.04141408614403268)nd2[&state=2]:1.225238932333245,((nd6[&state=1]:0.15152443528855317,nd7[&state=1]:0.15152443528855317)nd4[&state=1]:0.16841576897032884,(nd8[&state=1]:0.11388809491478649,nd9[&state=1]:0.0982030908340934)nd5[&state=1]:0.2060521093440955)nd3[&state=1]:0.9611550707385709)root[&state=0]:0.7189047250025469)origin[&state=0]:0.0;"

    at_dict = \
        {"nd3":
        [
            pjat.AttributeTransition(
                "state",
                "nd3",
                1.6028010851528143,
                0,
                2),
            pjat.AttributeTransition(
                "state",
                "nd3",
                1.6714902163855656,
                2,
                1)
        ],
        "nd2": [pjat.AttributeTransition("state", "nd2", 1.884598041236693, 0, 2)],
        "nd10": [pjat.AttributeTransition("state", "nd10", 1.944143657335792, 2, 1)],
        "nd12": [pjat.AttributeTransition("state", "nd12", 1.9855577434798246, 2, 0)]
        }
    
    try:
        ann_tr = pjr.read_nwk_tree_str(tr_str,
                                       in_file=False)

    except (TypeError,
            AttributeError,
            dp.dataio.tokenizer.Tokenizer.UnexpectedEndOfStreamError,
            dp.dataio.newickreader.NewickReader.NewickReaderMalformedStatementError,
            dp.dataio.newickreader.NewickReader.NewickReaderDuplicateTaxonError) as e:
        print(e)

    ann_tr.at_dict = at_dict
    ann_tr.populate_nd_attr_dict(["state"], attr_dict_added_separately_from_tree=True)
    ann_tr.state_count = 3 # need it for plotting colors

    print(ann_tr.tree.as_string(schema="newick"))

    pjtr.plot_ann_tree(
        ann_tr,
        ax,
        use_age=False,
        start_at_origin=True,
        sa_along_branches=False,
        attr_of_interest="state")


def run_example_map_attr(ax: matplotlib.pyplot.Axes) -> None:
    """Run example 6.

    Returns:
    """

    dag_obj = cmdp.script2dag("examples/see_stoch_maps.pj", \
                              in_pj_file=True)

    for node_name, node_dag in dag_obj.name_node_dict.items():
        if node_name == "mapped_trs":
            i = 0 # 1
            ann_tr = node_dag.value[i]
            pjtr.plot_ann_tree(
                ann_tr,
                ax,
                use_age=False,
                start_at_origin=False,
                sa_along_branches=False,
                attr_of_interest="state")

# TODO
# def run_example_inference_string():
    # as if we had clicked "See" in the inference tab
    # all_sims_model_spec_list, all_sims_mcmc_logging_spec_list, dir_list = \
    #     pjinf.dag_obj_to_rev_inference_spec(
    #         dag_obj,
    #         "./inference_test/",
    #         mcmc_chain_length=1000,
    #         prefix="test")
    # for i, ith_sim_model_spec in enumerate(all_sims_model_spec_list):
    #     print(ith_sim_model_spec)
    #     print(all_sims_mcmc_logging_spec_list[i])

    # pjio.output_inference_rev_scripts(
    #     all_sims_model_spec_list,
    #     all_sims_mcmc_logging_spec_list,
    #     dir_list,
    #     prefix="test")

if __name__ == "__main__":

    dag_obj: pgm.DirectedAcyclicGraph = pgm.DirectedAcyclicGraph()

    #####################################
    # Preparing plotting when necessary #
    #####################################
    # fig = Figure(figsize=(11,4.5))

    # note that pjgui uses matplotlib.figure.Figure
    # (which is part of Matplotlib's OOP class library)
    #
    # here, we instead use pyplot's figure, which is the
    # Matlab-like state-machine API
    fig = matplotlib.pyplot.figure()

    ax = fig.add_axes([0.25, 0.2, 0.5, 0.6])
    ax.patch.set_alpha(0.0)
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # examples
    # 1: Yule model string, 2 tree samples, 2 tree replicates per sample
    # 2: Same as (1), but reading a pre-made .pj script in examples/geosse.pj
    # 3: BiSSE model with incomplete sample, 2 tree samples, 2 tree replicates per sample
    # 4: Builds discrete SSE tree manually, then prints on screen
    # 5: Builds discrete SSE tree from Newick string, then prints on screen
    # 6: Read .pj script examples/see_stoch_maps.pj

    example_to_run = 1
    # example_to_run = 2
    # example_to_run = 3
    # example_to_run = 4
    # example_to_run = 5
    # example_to_run = 6
        
    if example_to_run == 1:
        dag_obj = run_example_yule_string()

    elif example_to_run == 2:
        dag_obj = run_example_geosse_file()

    elif example_to_run == 3:
        dag_obj = run_example_manual_incomplete_sampling_bisse()

    elif example_to_run == 4:
        run_example_manual_tree_building(ax)
        matplotlib.pyplot.show()

    elif example_to_run == 5:
        run_example_read_tree_function(ax)
        matplotlib.pyplot.show()

    elif example_to_run == 6:
        run_example_map_attr(ax)
        matplotlib.pyplot.show()

    # printing to screen
    if example_to_run in ([1, 2, 3]):
        for node_name, node_dag in dag_obj.name_node_dict.items():
            if isinstance(node_dag, pgm.StochasticNodeDAG):
                if isinstance(node_dag.value[0], pjtr.AnnotatedTree):
                    print(node_dag.value[0].tree.as_string(schema="newick"))
                    print(node_dag.get_node_stats_str(0, len(node_dag.value), 0))

                    # if one wants to see cladogenetic attribute transitions
                    for nd_name, clado_at in node_dag.value[0].clado_at_dict.items():
                        print(nd_name, "\n", clado_at)