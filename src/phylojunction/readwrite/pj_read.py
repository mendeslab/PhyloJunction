import os
import typing as ty
import numpy as np
import pandas as pd  # type: ignore
import pickle  # type: ignore
import dendropy as dp  # type: ignore
import string
import csv
import collections

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.utility.exception_classes as ec
import phylojunction.utility.helper_functions as pjh
import phylojunction.data.sampled_ancestor as pjsa
from phylojunction.data.tree import AnnotatedTree

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


def read_text_file(fp_string: str) -> ty.List[str]:
    """Read and parse text file into list of strings (one per line).

    Args:
        fp_string (str): String containing file path to text file being
            read

    Returns:
        str: List of strings, each being a line of the input text file
    """

    cmd_line_list = list()

    with open(fp_string, "r") as infile:
        for line in infile:
            line = line.rstrip()

            if line:
                cmd_line_list.append(line)

    return cmd_line_list


def read_serialized_pgm(fp_string: str) -> pgm.DirectedAcyclicGraph:
    """Read binary file storing PGM from a previous PJ session.

    Args:
        fp_string (str): String containing file path to binary file
            storing DAG from previous PJ session.

    Returns:
        (DirectedAcyclicGraph): DAG object to be initialized.
    """

    with open(fp_string, "rb") as picklefile:
        dag_obj: pgm.DirectedAcyclicGraph = pickle.load(picklefile)

    return dag_obj


def read_csv_tsv_into_dataframe(fp_string: str, is_file_csv: bool=True) -> pd.DataFrame:
    """Read .csv/.tsv file into a pandas DataFrame.

    Args:
        fp_string (str): String containing file path to text file being
            read
        is_csv (bool): Flag specifying if file is .csv (otherwise it is)
            a .tsv.

    Returns:
        pd.DataFrame: pandas DataFrame object (empty DataFrame if not
            CSV nor TSV)
    """

    if is_file_csv:
        if is_csv(fp_string):
            return pd.read_csv(fp_string)

        else:
            return pd.DataFrame()
        
    else:
        if is_tsv(fp_string):
            return pd.read_csv(fp_string, delimiter="\t")

        else:
            return pd.DataFrame()    


def is_csv(fp_string: str) -> bool:
    """Check if file in provided path is in CSV format.

    Args:
        fp_string (str): String containing file path to text file being
            read

    Returns:
        bool: True if file is in CSV format, False otherwise
    """

    try:
        with open(fp_string, newline="") as csvfile:
            start = csvfile.read(4096)  # reading file in 4 KiB chunks

            # isprintable does not allow newlines, printable does not allow
            # umlauts...
            if not all([c in string.printable or c.isprintable()
                        for c in start]):
                print(("Inside readwrite.pj_read.is_csv():\n    ERROR: "
                       "Attempted to read CSV file, but it does not seem "
                       "to be in CSV format."))

                return False

            dialect = csv.Sniffer().sniff(start)

            return True

    except csv.Error:
        # could not get a csv dialect -> probably not a csv.
        print(("Inside readwrite.pj_read.is_csv():\n    ERROR: Attempted to "
               "read CSV file, but it does not seem to be in CSV format."))

        return False


def is_tsv(fp_string: str) -> bool:
    try:
        with open(fp_string, "r") as file_in:
            reader = csv.reader(file_in, delimiter="\t")
            headers = next(reader)
            column_count = len(headers)
            
            if column_count == 1:
                raise ValueError("Column count is 1")
            
            for line in reader:
                if len(line) != column_count:
                    raise ValueError("Column count was variable")
                
            return True

    except ValueError:
        print(f"file \"{fp_string}\" is probably not a .tsv.")

        return False


def parse_cli_str_write_fig(str_write_fig: str) \
        -> ty.Dict[str, ty.Tuple[int]]:
    """Parse command-line string argument for generating node plots.

    Args:
        str_write_fig (str): User-provided string as argument to command-line
            interface -f parameter (e.g., 'tr;0-10')
    
    Returns:
        Dictionary with node names as keys (str), tuple with range start
            (int, int)
    """

    node_range_dict: ty.Dict[str, ty.Tuple[int]] = dict()

    node_range_list = str_write_fig.split(";")

    # only one node name passed, we assume all samples
    # are needed, range tuple is empty
    if len(node_range_list) == 1:
        if node_range_list[0].isnumeric():
            raise ec.PJCLIInvalidInput(
                "-f",
                ("If no range is provided, you must input a node name "
                 "(string) whose figure to plot, instead of an integer.")
            )
            # raise RuntimeError

        node_range_dict[node_range_list[0]] = tuple([])

    else:
        node_names_str, node_ranges_str = node_range_list
        node_names_list = node_names_str.split(",")
        node_ranges_list = node_ranges_str.split(",")

        if len(node_names_list) != len(node_ranges_list):
            raise ec.PJCLIInvalidInputError(
                "-f",
                ("If ranges are provided, the number of ranges must "
                 "match the number of node names to plot figures for.")
            )
            # raise RuntimeError

        for idx, node_name in enumerate(node_names_list):
            range_str = node_ranges_list[idx]

            try:
                range_tup = tuple(int(i) for i in range_str.split("-"))

            except ValueError:
                raise ec.PJCLIInvalidInputError(
                    "-f",
                    ("Ranges must be defined by non-negative integers."))
                # raise RuntimeError

            if node_name in node_range_dict:
                raise ec.PJCLIInvalidInputError(
                    "-f",
                    ("Node names appeared more than once. Nodes to plot "
                     "figures for must be unique.")
                )
                # raise RuntimeError

            node_range_dict[node_name] = range_tup

    return node_range_dict


def read_nwk_tree_str(nwk_tree_path_or_str: str,
                      fn_name: str = "read_tree",
                      in_file: bool = True,
                      node_names_attribute: str = "",
                      n_states: int = 1,
                      is_complete: bool = False,
                      keep_unsampled_alive: bool = False,
                      epsilon: float = 1e-12) \
        -> AnnotatedTree:
    """Read Newick tree string directly or in provided file.

    Args:
        nwk_tree_path_or_str (str): Tree Newick string, or path to file
            containing single tree Newick string.
        fn_name (str): Name of the function (in the .pj script) being
            called.
        in_file (bool): If tree string is in a file being passed as
            argument (True) or if Newick string is being passed directly
            (False). Defaults to 'True'.
        node_names_attribute (str): Defaults to empty string "".
        is_complete (bool): If 'False' (the default), the Newick string
            is treated as already being the *reconstructed* tree --
            every leaf is sampled by construction, no annotations
            needed *or* consulted at all: 'alive' is inferred purely
            from depth (a leaf at the tree's maximum depth is extant;
            any other leaf is a fossil sample -- no '&alive='
            annotation required, or looked at, either way), and
            'is_sa' from branch length (a zero-length leaf is a direct
            ancestor). Because neither is ever actually read from an
            annotation in this mode, an explicit '&alive=' of either
            value on any leaf is rejected outright (not only when it
            happens to disagree with the depth-based inference), and
            so is an explicit '[&sampled=False]' (every leaf of an
            already-reconstructed tree is sampled, by definition) --
            both raise 'AnnotatedTreeIncorrectAnnotationError' rather
            than being silently ignored, since either is most likely
            a sign that 'is_complete="true"' should have been passed
            instead. If 'True', the
            Newick string is treated as a *complete* tree (as a
            simulator would produce, including lineages that were
            never sampled) -- PhyloJunction will not try to infer
            which terminal branches are alive or sampled from branch
            lengths alone (not reliable for arbitrary, possibly non-
            ultrametric, externally-provided data); instead, every
            non-direct-ancestor leaf must carry an explicit
            '[&alive=True]'/'[&alive=False]' NHX annotation, and, for
            any leaf so annotated 'alive=False', an explicit
            '[&sampled=True]'/'[&sampled=False]' one too (to
            distinguish a sampled fossil tip from an unsampled
            extinction event) -- see 'keep_unsampled_alive' for the
            symmetric case on the 'alive=True' side. A missing
            required annotation raises
            'AnnotatedTreeMissingNodeAnnotationError'; an explicit
            annotation that contradicts a structural invariant (an
            internal node marked alive/sampled, or a direct ancestor
            marked alive) raises 'AnnotatedTreeIncorrectAnnotationError'.
        keep_unsampled_alive (bool): Only consulted when
            'is_complete' is 'True'. An alive leaf explicitly
            annotated '[&sampled=False]' is, by default ('False'),
            treated the same way an unsampled extinction would be
            (i.e., it will be pruned away when the reconstructed tree
            is extracted). If 'True', such a leaf is instead kept
            ('sampled' forced back to 'True'), letting a known-but-
            unsampled living taxon be represented without it being
            discarded from the reconstructed tree.

    Returns:
        AnnotatedTree with populated attributes dictionary
            (i.e., a dendropy.Tree that has been annotated)
    """

    def _get_bool_annotation(
            nd: dp.Node, attr_name: str, nd_name: str) -> ty.Optional[bool]:
        """Return a node's boolean NHX annotation, or 'None' if absent.

        Dendropy auto-casts bare 'True'/'False' Newick comment values
        (e.g. '[&alive=True]') into real Python bools already; this
        also tolerates the equivalent quoted/differently-cased string
        forms.

        Uses 'get_value' rather than bracket indexing
        ('nd.annotations[attr_name]') deliberately: the latter
        auto-*creates* a new, blank annotation on 'nd' as a side
        effect of merely looking it up when absent (dendropy's
        'AnnotationSet.__getitem__' does this), which would otherwise
        leave every internal node (this function is also called, for
        validation, on internal nodes that never had '&alive='/
        '&sampled=' in the Newick string at all) with a stray,
        empty '[&alive=,sampled=]' annotation the next time this
        tree's string representation is generated.
        """

        ann_val = nd.annotations.get_value(attr_name, "")

        if ann_val == "":
            return None

        if isinstance(ann_val, bool):
            return ann_val

        ann_str = str(ann_val).strip().lower()

        if ann_str == "true":
            return True

        if ann_str == "false":
            return False

        raise ec.AnnotatedTreeIncorrectAnnotationError(
            "Node '" + nd_name + "'s '&" + attr_name + "=' annotation "
            "must be 'True' or 'False', but was '" + str(ann_val) + "'.")

    def _reject_contradictory_annotation(
            nd: dp.Node, attr_name: str, nd_name: str,
            forbidden_value: bool, reason: str) -> None:
        """Raise if 'nd' has an explicit, contradictory annotation."""

        if _get_bool_annotation(nd, attr_name, nd_name) is forbidden_value:
            raise ec.AnnotatedTreeIncorrectAnnotationError(
                "Node '" + nd_name + "' was annotated '&" + attr_name +
                "=" + str(forbidden_value) + "', but " + reason + ".")

    tr_str: str = str()
    node_attr_dict = collections.defaultdict(
        pjh.create_str_defaultdict)
    node_id_to_name: dict[str, str] = dict()

    # debugging
    # print("node_names_attribute =", node_names_attribute)

    # if tree provided in separate file
    if in_file:
        if not os.path.isfile(nwk_tree_path_or_str):
            raise ec.PJIOFileDoesNotExistError(fn_name,
                                               nwk_tree_path_or_str)
        
        with open(nwk_tree_path_or_str, "r") as infile:
            tr_str = infile.readline()

        dp_tr = dp.Tree.get(data=tr_str, schema="newick")

    # if tree is provided as newick string
    else:
        dp_tr = dp.Tree.get(data=nwk_tree_path_or_str, schema="newick")

    # first we deal with having a root vs an origin
    is_origin = (True if (len(dp_tr.seed_node.child_nodes()) == 1)# or \
                          # dp_tr.seed_node.edge_length > 0.0)
                 else False)
    
    root_height = dp_tr.max_distance_from_root()
    seen_root = False
    root_alternative_label = ""

    # debugging
    # print("root height = " + str(root_height))
    # print("is_origin = " + str(is_origin))

    nd_count: int = 1

    # now we annotate nodes
    for nd in dp_tr.preorder_node_iter():
        # if node ids (e.g. index) are provided
        # as node_names_attribute, we will use
        # them to name every single node in the tree
        with_attribute: bool = False
        nd_name: str = "nd" + str(nd_count)
        node_names_attribute_val: str = ""
        if node_names_attribute != "":
            node_names_attribute_val = nd.annotations[node_names_attribute].value

            if node_names_attribute_val == "":
                raise ec.ObjInitInvalidArgError(
                    fn_name,
                    "'node_name_attribute'",
                    message=("The attribute for node names "
                             "was not in the newick strings."))

            # we make the node have that attribute!
            setattr(nd, node_names_attribute, node_names_attribute_val)

            # now we also use that attribute value as a node name
            nd_name = "nd" + node_names_attribute_val
            with_attribute = True

        # used in error messages below -- prefer the node's own,
        # already-parsed taxon label (e.g. "fossiltip") over the
        # generic, auto-numbered 'nd_name' placeholder, which is only
        # ever turned into this node's actual label/taxon later on,
        # well after the validation calls that need to identify it
        nd_identity = nd.taxon.label if nd.taxon is not None else nd_name

        nd.is_sa_lineage = False
        nd.is_sa_dummy_parent = False

        # (1) seed node (origin and/or root)
        if nd is dp_tr.seed_node:
            nd.state = 0
            nd.is_sa = False
            nd.alive = False
            nd.sampled = False

            # same dummy-parent-of-a-direct-ancestor detection as any
            # other internal node gets below (see "(2) not seed node")
            # -- a direct ancestor can perfectly well branch off right
            # at the seed itself (e.g. a fossil sampled very early on,
            # with the rest of the tree as its only sibling), and
            # leaving this hardcoded to 'False' regardless left
            # 'get_y_coord_from_n_obs_nodes' with no node left to
            # (mis-)compute a y-coordinate for: a seed node with
            # exactly one non-SA child is neither recognized as a
            # dummy parent (so the "average with origin" branch does
            # not apply) nor has two children itself (so the
            # bifurcating-node branch does not apply either), leaving
            # it silently missing from the returned dict -- surfacing
            # downstream as a 'KeyError' in '_draw_clade' when
            # plotting the complete tree
            nd.is_sa_dummy_parent = False
            for ch_nd in nd.child_node_iter():
                if abs(ch_nd.edge_length) <= epsilon:
                    nd.is_sa_dummy_parent = True
                    break

            # the seed node is itself an internal node (the origin or
            # the root) -- it can never be alive nor sampled, so an
            # explicit, contradictory annotation on it is rejected
            # the same way as for any other internal node
            if is_complete:
                _reject_contradictory_annotation(
                    nd, "alive", nd_identity, True,
                    "internal nodes (including the origin/root) are "
                    "never alive")
                _reject_contradictory_annotation(
                    nd, "sampled", nd_identity, True,
                    "internal nodes (including the origin/root) are "
                    "never sampled")

            if not with_attribute:
                nd_name = ("origin" if is_origin else "root")

            # the seed node is the root itself whenever there is no
            # separate origin above it -- true regardless of
            # 'with_attribute', unlike the 'nd_name == "root"' check
            # this replaces, which only ever matched when
            # 'with_attribute' was False (since 'nd_name' is instead
            # the node_names_attribute-derived name, e.g. "nd0", in
            # that case) -- leaving 'seen_root' permanently False and
            # letting the "root!" block below (for sub-seed internal
            # nodes) misfire on every later non-dummy-parent internal
            # node it saw, repeatedly overwriting
            # 'root_alternative_label' with the wrong node
            if not is_origin:
                seen_root = True

                if with_attribute:
                    root_alternative_label = nd_name

            if nd.taxon is None:
                nd.taxon = dp.Taxon(label=nd_name)
                nd.label = nd_name
                dp_tr.taxon_namespace.add_taxon(nd.taxon)

        # (2) not seed node
        else:
            nd_count += 1

            # internal node
            if not nd.is_leaf():
                # initial values
                nd.state = 0
                nd.is_sa_dummy_parent = False
                nd.is_sa = False
                nd.alive = False
                nd.sampled = False

                # seeing if origin child is root
                # or the dummy parent of a sampled ancestor
                for ch_nd in nd.child_node_iter():
                    if abs(ch_nd.edge_length) <= epsilon:
                        nd.is_sa_dummy_parent = True
                        nd.is_sa = False
                        break

                # root!
                if not nd.is_sa_dummy_parent and \
                        not seen_root:
                    if not with_attribute:
                        nd.taxon = dp.Taxon(label="root")
                        nd.label = "root"
                        seen_root = True
                    # dp_tr.taxon_namespace.add_taxon(nd.taxon)
                    root_alternative_label = nd_name

                # internal nodes are never alive nor sampled -- reject
                # an explicit, contradictory annotation on this node
                # (never inferred from one: 'nd.alive'/'nd.sampled'
                # above are always 'False' for internal nodes,
                # regardless)
                if is_complete:
                    _reject_contradictory_annotation(
                        nd, "alive", nd_identity, True,
                        "internal nodes are never alive")
                    _reject_contradictory_annotation(
                        nd, "sampled", nd_identity, True,
                        "internal nodes are never sampled")

            # leaf node
            else:
                # direct (sampled) ancestor: zero-length branch --
                # by definition always sampled, never alive,
                # regardless of 'is_complete'/annotations
                if abs(nd.edge_length) <= epsilon:
                    nd.is_sa_dummy_parent = False
                    nd.is_sa = True
                    nd.sampled = True
                    nd.alive = False

                    if is_complete:
                        _reject_contradictory_annotation(
                            nd, "alive", nd_identity, True,
                            "direct (sampled) ancestors are never "
                            "alive")
                        _reject_contradictory_annotation(
                            nd, "sampled", nd_identity, False,
                            "direct (sampled) ancestors are always "
                            "sampled, by definition")

                    else:
                        _reject_contradictory_annotation(
                            nd, "sampled", nd_identity, False,
                            "a non-complete tree ('is_complete=\"false\"', "
                            "the default) represents an already-"
                            "reconstructed process, so every leaf -- "
                            "including a direct (sampled) ancestor -- "
                            "must be sampled, by definition")

                # ordinary terminal branch (non-zero length)
                else:
                    nd.is_sa = False

                    # complete tree: branch lengths/depth alone are
                    # not trustworthy (e.g. non-ultrametric real
                    # data), so 'alive'/'sampled' must come from
                    # explicit annotations instead of being inferred
                    if is_complete:
                        alive_ann = _get_bool_annotation(
                            nd, "alive", nd_identity)

                        if alive_ann is None:
                            raise ec.AnnotatedTreeMissingNodeAnnotationError(
                                nd_identity, "alive")

                        nd.alive = alive_ann

                        if nd.alive:
                            sampled_ann = _get_bool_annotation(
                                nd, "sampled", nd_identity)

                            # an explicit '&sampled=False' on an alive
                            # leaf is honored (i.e., this leaf will be
                            # pruned away like an unsampled extinction
                            # would be) unless the caller explicitly
                            # asked to keep it anyway
                            if sampled_ann is False and \
                                    not keep_unsampled_alive:
                                nd.sampled = False

                            else:
                                nd.sampled = True

                        else:
                            sampled_ann = _get_bool_annotation(
                                nd, "sampled", nd_identity)

                            if sampled_ann is None:
                                raise \
                                    ec.AnnotatedTreeMissingNodeAnnotationError(
                                        nd_identity, "sampled")

                            nd.sampled = sampled_ann

                    # reconstructed tree (default): every leaf is
                    # already sampled by construction; 'alive' is
                    # inferred from depth (True at the tree's maximum
                    # depth, False otherwise) -- neither is ever
                    # consulted from an annotation in this mode, so an
                    # explicit '&alive=' of either value is rejected
                    # outright (not just when it happens to disagree
                    # with the depth-based inference), the same way
                    # an explicit '&sampled=False' already is -- both
                    # are most likely a sign that 'is_complete="true"'
                    # should have been passed instead
                    else:
                        if _get_bool_annotation(
                                nd, "alive", nd_identity) is not None:
                            raise ec.AnnotatedTreeIncorrectAnnotationError(
                                "Node '" + nd_identity + "' was annotated "
                                "with an explicit '&alive=' value, but a "
                                "non-complete tree ('is_complete=\"false\"', "
                                "the default) never consults it -- 'alive' "
                                "is always deduced from depth instead. "
                                "This usually means 'is_complete=\"true\"' "
                                "should have been passed instead.")

                        _reject_contradictory_annotation(
                            nd, "sampled", nd_identity, False,
                            "a non-complete tree ('is_complete=\"false\"', "
                            "the default) represents an already-"
                            "reconstructed process, so every leaf must "
                            "be sampled, by definition")

                        nd.sampled = True

                        if abs(nd.distance_from_root() - root_height) \
                                <= epsilon:
                            nd.alive = True

                        else:
                            nd.alive = False

            # Only the continuing child is an SA lineage. It can itself be a dummy
            # parent when two sampling events occur on the same lineage.
            if nd.parent_node.is_sa_dummy_parent and not nd.is_sa:
                nd.is_sa_lineage = True

            # internal and no name provided
            if nd.taxon is None:
                if nd.label is None:
                    nd.label = nd_name
                    nd.taxon = dp.Taxon(label=nd_name)

                else:
                    nd.taxon = dp.Taxon(label=nd.label)

                dp_tr.taxon_namespace.add_taxon(nd.taxon)

            if nd.is_leaf():
                nd.state = 0

                # if node_names_attribute is passed,
                # we rename all taxa, including leaves
                if with_attribute:
                    nd.label = nd_name
                    nd.taxon.label = nd_name
                
                else:
                    nd.label = nd.taxon.label

            # # internal and no name provided
            # elif nd.taxon is None:                
            #     if nd.label == "None":
            #         nd.label = nd_name
            #         nd.taxon = dp.Taxon(label=nd_name)

            #     else:
            #         nd.taxon = dp.Taxon(label=nd.label)

            #     print(dp_tr.taxon_namespace)
            #     dp_tr.taxon_namespace.add_taxon(nd.taxon)
            #     print(dp_tr.taxon_namespace)

            # node_id_to_name[nd_id] = nd.label
            node_attr_dict[nd.label][node_names_attribute] = \
                node_names_attribute_val

    # debugging
    # print(dp_tr.taxon_namespace)

    # for nd in dp_tr.postorder_node_iter():
    #     print(nd is dp_tr.seed_node)
    #     print(nd.parent_node)

    # build 'sa_lineage_dict' *before* constructing 'ann_tr' -- every
    # direct (sampled) ancestor was already identified above (a
    # zero-length leaf). 'AnnotatedTree.__init__' derives
    # 'n_sa_nodes'/'sa_obs_nodes_labels' from this dict automatically
    # (see '_count_sampled_ancestors'); without it, both would
    # silently stay at 0/empty, which has a real, visible consequence
    # beyond just those two attributes' own display:
    # 'tree.get_y_coord_from_n_obs_nodes's phylogram-specific
    # ('sa_along_branches=False') branches used to size their plot
    # directly off of 'n_sa_nodes'.
    #
    # 'tree.plot_ann_tree's own SA-dot-drawing code ('_draw_clade')
    # also looks up 'sa_lineage_dict'/'rec_tr_sa_lineage_dict' (the
    # latter derived from the former during
    # 'extract_reconstructed_tree') for both the complete and
    # reconstructed tree views -- direct ancestors are deliberately
    # excluded from getting their own row there (on the assumption
    # that they are drawn as dots along their lineage's branch
    # instead -- see 'sa_along_branches'), so without a populated
    # 'sa_lineage_dict', a Newick-read tree's direct ancestors would
    # silently disappear from both views entirely (no row, and no
    # dot either).
    sa_lineage_dict: ty.DefaultDict[str, ty.List[pjsa.SampledAncestor]] = \
        collections.defaultdict(list)

    for nd in dp_tr.leaf_node_iter():
        if not nd.is_sa:
            continue

        # the lineage node is this SA's sibling under their shared
        # dummy parent -- the *other* child, whichever it is
        lineage_sibling = None
        for sib in nd.parent_node.child_nodes():
            if sib is not nd:
                lineage_sibling = sib
                break

        if lineage_sibling is None:
            continue

        global_time = nd.distance_from_root()

        sa_lineage_dict[lineage_sibling.label].append(
            pjsa.SampledAncestor(
                nd.label,
                lineage_sibling.label,
                global_time,
                age=root_height - global_time,
                time_to_lineage_node=lineage_sibling.edge_length))

    # note the empty attribute transition dict members
    ann_tr = AnnotatedTree(dp_tr,
                           n_states,
                           alternative_root_label=root_alternative_label,
                           start_at_origin=is_origin,
                           tree_died=False,
                           read_as_newick_string=True,
                           at_dict=dict(),
                           clado_at_dict=dict(),
                           sa_lineage_dict=dict(sa_lineage_dict))

    ann_tr.node_attr_dict = node_attr_dict

    return ann_tr

    # print(dp_tr.taxon_namespace)
    # print(node_attr_dict)
    # print(node_id_to_name)


def read_node_attr_update_tree(attr_tsv_path: str,
                               attr_name: str,
                               attr_cast: ty.Callable,
                               ann_tr: AnnotatedTree) -> None:
    """Update AnnotatedTree members with attribute information.

    Args:
        attr_tsv_path (string): Path to .tsv file containing attribute
            values for nodes in the tree.
        attr_name: Name of the attribute.
        attr_cast (ty.Callable): Type to be used to cast value
            as (e.g., int, float).
        ann_tr (AnnotatedTree): AnnotatedTree object to be updated
            with attribute name and values
    """

    # nodes may be just tips, or also internal nodes, but they
    # must match the names in the namespace below
    #
    # NOTE: if the AnnotatedTree was read in by passing a tree
    # in Newick string, the node label for the root will be 'root'.
    # This may cause a problem if the .tsv file below provides
    # an attribute for the root, but the root node in that file
    # is not labeled 'root'
    node_name_attr_str_list = read_text_file(attr_tsv_path)

    ann_tr_nd_names = \
        set([nd.label for nd in ann_tr.tree.taxon_namespace])
    
    for line in node_name_attr_str_list:
        nd_name, nd_attr_val = line.split("\t")

        # if node name is not in tree, something must be wrong
        if nd_name not in ann_tr_nd_names:
            exit("Need to add exception here.")

        nd = ann_tr.tree.find_node_with_label(nd_name)
        # set the attribute value
        nd.__setattr__(attr_name, attr_cast(nd_attr_val))

    ann_tr.populate_nd_attr_dict([attr_name],
                                 attr_dict_added_separately_from_tree=False)


if __name__ == "__main__":

    # should be ok
    # print(parse_cli_str_write_fig("tr"))
    # print(parse_cli_str_write_fig("tr;0"))
    # print(parse_cli_str_write_fig("tr;0-10"))
    # print(parse_cli_str_write_fig("rv1,tr;0,0"))
    # print(parse_cli_str_write_fig("rv1,tr;0-1,0"))
    # print(parse_cli_str_write_fig("rv1,tr;0-1,0-1"))

    # should raise hell (try one at a time)
    # parse_cli_str_write_fig("0")
    # parse_cli_str_write_fig("rv;tr")
    # parse_cli_str_write_fig("rv,tr;0")
    # parse_cli_str_write_fig("rv,tr;0,")

    # should all work!
    # tr = read_nwk_tree_str("examples/trees_maps_files/turtle.tre",
    # tr = read_nwk_tree_str("examples/trees_maps_files/geosse_dummy_tree1.tre",
    #                        node_names_attribute="index")
    tr = read_nwk_tree_str("examples/trees_maps_files/geosse_dummy_tree3.tre",
                           node_names_attribute="index")

    # print(tr.node_attr_dict)

    read_node_attr_update_tree("examples/trees_maps_files/geosse_dummy_tree3_tip_states.tsv",
    # read_node_attr_update_tree("examples/trees_maps_files/geosse_dummy_tree1_tip_states.tsv",
                               "state",
                               int,
                               tr)

    
    # dummy trees should be:
    # (((sp1:1.0,sp2:1.5)nd4:1.0,sp3:1.7)root:0.5)origin:0.5
    # ((sp1:1.0,sp2:1.5)nd4:1.0,sp3:1.7)root
    # (((sp1:1.0,sp2:1.5)root:1.0,sp3:0.0)nd2:0.5)origin:0.5