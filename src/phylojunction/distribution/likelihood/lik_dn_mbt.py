from phylojunction.data.trait import DiscreteTrait
import typing as ty
import dendropy as dp
import numpy as np
import math

# pj imports
from phylojunction.data.tree import AnnotatedTree
import phylojunction.distribution.likelihood.ODE.dn_mbt_ode as pjmbt
import phylojunction.distribution.likelihood.ODE.dn_bisse_ode as pjbisse

def recursively_do_node(nd: dp.Node,
                        nd_age_dict: ty.Dict[str, float],
                        pars: ty.Tuple[np.matrix, np.ndarray],
                        n_states: int,
                        log_norm_factors: ty.List[float]) -> np.ndarray:

    # get parameters
    qs, mus, lambdas = pars

    # get node inf
    nd_name = nd.label
    nd_age = nd_age_dict(nd_name)

    # get integration interval info

    # MBT
    t_start = nd_age
    t_end = nd.edge_length + nd_age
    # simple BiSSE
    # t_start = 0
    # t_end = nd.edge_length

    children_ds_es: ty.List[np.ndarray] = list()

    # if there are children, recur
    for idx, ch_nd in enumerate(nd.child_nodes()):
        # we ignore direct (sampled) ancestors
        if not ch_nd.is_sa:
            ch_ds_es, log_norm_factors = \
                recursively_do_node(ch_nd,
                                    nd_age_dict,
                                    pars,
                                    n_states,
                                    log_norm_factors)
            children_ds_es.append(ch_ds_es)

    ######################
    # Doing current node #
    ######################

    # debugging
    # print("doing nd", nd.label)

    # new Ds and Es for every node, terminal or internal
    ds_es = np.zeros(2 * n_states)

    # if internal node
    # (combine children if neither is a direct ancestor)
    if nd.num_child_nodes() == 2:
        left_ds_es, right_ds_es = children_ds_es

        # debugging
        # print("ds_es (before combining)\n", ds_es)
        # print("left_ds_es", left_ds_es)
        # print("right_ds_es", right_ds_es)

        for i in range(n_states):
            # combine Ds
            ds_es[i] = lambdas[i] * left_ds_es[i] * right_ds_es[i]

            # do Es
            ds_es[n_states + i] = left_ds_es[n_states + i]


    # if tip, initialize ds_es to observed state
    elif nd.is_leaf():
        st = nd.state

        # state is ambiguous!
        if st == "?":
            # all Ds is 1, Es remain zero
            for i in range(n_states):
                ds_es[i] = 1

        # state is observed
        else:
            st = int(st)
            # obs D is 1, Es remain zero
            ds_es[st] = 1

    # debugging
    # print("ds_es (after combining if internal)", ds_es)
    # print("t_start", t_start, "t_end", t_end)

    # if not root, we have a branch to integrate over
    if nd.parent_node is not None:
        # ds_es = \
        #     pjbisse.solve_bisse_ds_es(ds_es,
        #                               t_start,
        #                               t_end,
        #                               qs,
        #                               mus,
        #                               lambdas,
        #                               verbose=False)

        ds_es = \
            pjmbt.solve_mbt_ds_es(ds_es,
                                  t_start,
                                  t_end,
                                  qs,
                                  mus,
                                  lambdas,
                                  verbose=False)

        norm_factor = np.sum(ds_es)

        # normalize ds_es
        ds_es /= norm_factor

        log_norm_factors.append(math.log(norm_factor))

    # debugging
    # print("returned ds_es", ds_es, "\n")

    return ds_es, log_norm_factors


def prune_mbt(ann_tr: AnnotatedTree,
              pars: ty.Tuple[np.ndarray],
              n_states: int,
              cond_on_survival: bool = True,
              correct_tip_shuffling: bool = False):

    # organize parameters
    pi = pars[-1]
    pars = pars[:-1]

    # initialize ds_es and log_norm_factors
    ds_es = np.zeros(2 * n_states)
    log_norm_factors: ty.List[float] = list()

    # get origin or root node
    seed_nd = ann_tr.origin_node if ann_tr.with_origin \
        else ann_tr.root_node
    nd_age_dict = ann_tr.node_ages_dict

    norm_factors: ty.List[float] = list()

    # ds_es and norm_factors are initialized at tips
    ds_es, log_norm_factors = \
        recursively_do_node(seed_nd,
                            nd_age_dict,
                            pars,
                            n_states,
                            log_norm_factors)

    # TODO: now do stuff at root depending on conditioning
    lk = 0.0

    # looking at Ds
    for i in range(n_states):
        lk += pi[i] * ds_es[i]

    log_lk = math.log(lk)

    # put back normalization factors
    # print("log-norm factor sum", sum(log_norm_factors))
    log_lk += sum(log_norm_factors)

    if cond_on_survival:
        pass

    return log_lk


if __name__ == "__main__":

    n_states = 2

    #####################
    # Initializing tree #
    #####################
    # '?' signalizes ambiguous state (represented as '-1' inside AnnotatedTree's dict members)
    tr_str = "((A:0.0001[&state=1],B:0.0001[&state=0])nd1:0.0001[&state=?],C:0.0002[&state=0])root[&state=?];"
    dp_tree = dp.Tree.get(data=tr_str, schema="newick")

    # preparing DendroPy tree to be fed into AnnotatedTree
    for nd in dp_tree:
        setattr(nd, "is_sa", False)
        setattr(nd, "is_sa_dummy_parent", False)
        setattr(nd, "is_sa_lineage", False)

        if nd.is_leaf():
            setattr(nd, "alive", True)
            setattr(nd, "sampled", True)

        else:
            setattr(nd, "alive", False)
            setattr(nd, "sampled", False)

        # making sure all nodes have labels
        if nd.label is None:
            nd.label = str(nd.taxon)

        # setting states
        state_val = nd.annotations["state"].value
        setattr(nd, "state", state_val)

    ann_tr = AnnotatedTree(dp_tree,
                           DiscreteTrait(n_states),
                           start_at_origin=False,
                           max_age=2.0,
                           tree_died=False,
                           tree_invalid=False,
                           read_as_newick_string=True,
                           epsilon=1e-12)

    ###########################
    # Initializing parameters #
    ###########################

    qs = np.matrix([[-2., .9], [.001, -.2]])
    mus = np.array([.1, .1])
    lambdas = np.array([1.0, .099])
    pi = np.array([.5, .5])
    pars = (qs, mus, lambdas, pi)

    ###################
    # Get likelihood! #
    ###################

    log_lk = \
        prune_mbt(ann_tr,
                  pars,
                  n_states,
                  cond_on_survival=True,
                  correct_tip_shuffling=False)

    print("Final log-lk =", log_lk)

    # tr_str = "((A:0.0001[&state=1],B:0.0001[&state=0])nd1:0.0001[&state=?],C:0.0002[&state=0])root[&state=?];"
    # log-lk = -10.00966707670628
