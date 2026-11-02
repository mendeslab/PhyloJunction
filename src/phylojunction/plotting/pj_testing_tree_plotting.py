from phylojunction.data.trait import DiscreteTrait
# Time-heterogeneous GeoSSE tree
#
# (((nd10:0.05585634266420793[&state=1],(nd12:0.01444225652017525[&state=0],nd13:0.01444225652017525[&state=2])nd11:0.04141408614403268[&state=2])nd2:1.225238932333245[&state=2],((nd6:0.15152443528855317[&state=1],nd7:0.15152443528855317[&state=1])nd4:0.16841576897032884[&state=1],(nd8:0.11388809491478649[&state=1],nd9:0.0982030908340934[&state=1])nd5:0.2060521093440955[&state=1])nd3:0.9611550707385709[&state=1])root:0.7189047250025469[&state=0])origin:0.0[&state=0];
# at_dict = 
# nd3    Attribute ('state') transition:
#     Time: 1.6028010851528143
#     Subtending node: nd3
#     Departing state: 0
#     Arriving state: 2,
#        Attribute ('state') transition:
#     Time: 1.6714902163855656
#     Subtending node: nd3
#     Departing state: 2
#     Arriving state: 1
# nd2    Attribute ('state') transition:
#     Time: 1.884598041236693
#     Subtending node: nd2
#     Departing state: 0
#     Arriving state: 2
# nd10    Attribute ('state') transition:
#     Time: 1.944143657335792
#     Subtending node: nd10
#     Departing state: 2
#     Arriving state: 1
# nd12    Attribute ('state') transition:
#     Time: 1.9855577434798246
#     Subtending node: nd12
#     Departing state: 2
#     Arriving state: 0

import phylojunction.readwrite.pj_read as pjr
import phylojunction.data.attribute_transition as pjat
import phylojunction.data.tree as pjt

import dendropy as dp
import matplotlib
import matplotlib.pyplot as plt  # type: ignore


if __name__ == "__main__":
    # tr_str = "(((nd10:0.05585634266420793[&state=1],(nd12:0.01444225652017525[&state=0],nd13:0.01444225652017525[&state=2])nd11:0.04141408614403268[&state=2])nd2:1.225238932333245[&state=2],((nd6:0.15152443528855317[&state=1],nd7:0.15152443528855317[&state=1])nd4:0.16841576897032884[&state=1],(nd8:0.11388809491478649[&state=1],nd9:0.0982030908340934[&state=1])nd5:0.2060521093440955[&state=1])nd3:0.9611550707385709[&state=1])root:0.7189047250025469[&state=0])origin:0.0[&state=0];"
    tr_str1 = "(((nd10[&state=1]:0.05585634266420793,(nd12[&state=0]:0.01444225652017525,nd13[&state=2]:0.01444225652017525)nd11[&state=2]:0.04141408614403268)nd2[&state=2]:1.225238932333245,((nd6[&state=1]:0.15152443528855317,nd7[&state=1]:0.15152443528855317)nd4[&state=1]:0.16841576897032884,(nd8[&state=1]:0.11388809491478649,nd9[&state=1]:0.0982030908340934)nd5[&state=1]:0.2060521093440955)nd3[&state=1]:0.9611550707385709)root[&state=0]:0.7189047250025469)origin[&state=0]:0.0;"

    try:
        ann_tr1 = pjr.read_nwk_tree_str(tr_str1,
                                        in_file=False)

    except (TypeError,
            AttributeError,
            dp.dataio.tokenizer.Tokenizer.UnexpectedEndOfStreamError,
            dp.dataio.newickreader.NewickReader.NewickReaderMalformedStatementError,
            dp.dataio.newickreader.NewickReader.NewickReaderDuplicateTaxonError) as e:
        print(e)
    
    # print("ann_tr1.tree_read_as_newick", ann_tr1.tree_read_as_newick)

    # we will overwrite ann_tr1 current at_dict, which
    # is an empty dictionary at the moment
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

    ann_tr1.at_dict = at_dict
    ann_tr1.populate_nd_attr_dict(["state"], attr_dict_added_separately_from_tree=True)
    ann_tr1.trait = DiscreteTrait(3, name=ann_tr1.trait.name)
    
    fig = matplotlib.pyplot.figure()

    ax = fig.add_axes([0.25, 0.2, 0.5, 0.6])
    ax.patch.set_alpha(0.0)
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    pjt.plot_ann_tree(
        ann_tr1,
        ax,
        use_age=False,
        start_at_origin=True,
        sa_along_branches=False,
        attr_of_interest="state")
    
    plt.show()




    # tr_str2 = "((A:1.0,B:1.0):1.0,C:2.0);"
    tr_str2 = "examples/trees_maps_files/geosse_dummy_tree1.tre"

    try:
        # ann_tr2.at_dict is an empty dictionary
        ann_tr2 = pjr.read_nwk_tree_str(tr_str2,
                                        in_file=True)
    except (TypeError,
            AttributeError,
            dp.dataio.tokenizer.Tokenizer.UnexpectedEndOfStreamError,
            dp.dataio.newickreader.NewickReader.NewickReaderMalformedStatementError,
            dp.dataio.newickreader.NewickReader.NewickReaderDuplicateTaxonError) as e:
        print(e)

    # TODO: this function should take the type of the attribute value
    # e.g., "int" for discrete states (so type conversion can be carried out)
    pjr.read_node_attr_update_tree(
        "examples/trees_maps_files/geosse_dummy_tree1_tip_states.tsv",
        "state",
        ann_tr2)
    
    # for nd in ann_tr2.tree.nodes():
    #     print(nd.state)

    # TODO: should update state_count to 3 in read_node_attr_update_tree
    # print(ann_tr2.trait.states)

    # pjt.plot_ann_tree(
    #     ann_tr2,
    #     ax,
    #     use_age=False,
    #     start_at_origin=False,
    #     sa_along_branches=False,
    #     attr_of_interest="state")

    # plt.show()
