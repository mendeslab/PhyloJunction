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

    Returns:
        AnnotatedTree with populated attributes dictionary
            (i.e., a dendropy.Tree that has been annotated)
    """

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

        nd.is_sa_lineage = False
        nd.is_sa_dummy_parent = False

        # (1) seed node (origin and/or root)
        if nd is dp_tr.seed_node:
            nd.state = 0
            nd.is_sa_dummy_parent = False
            nd.is_sa = False
            nd.alive = False
            nd.sampled = False
            
            if not with_attribute:
                nd_name = ("origin" if is_origin else "root")

            if nd_name == "root":
                seen_root = True

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

            # annotate node as sampled ancestor or not
            if abs(nd.edge_length) <= epsilon:
                nd.is_sa_dummy_parent = False
                nd.is_sa = True
                nd.sampled = True

            else:
                nd.is_sa = False
                nd.sampled = False
        
            # annotate node as alive or not
            if abs(nd.distance_from_root() - root_height) <= epsilon:
                nd.alive = True
                nd.sampled = True

            else:
                nd.alive = False
                nd.sampled = False

            # annotate as sa_lineage
            if nd.parent_node.is_sa_dummy_parent:
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

    # note the empty attribute transition dict members
    ann_tr = AnnotatedTree(dp_tr,
                           n_states,
                           alternative_root_label=root_alternative_label,
                           start_at_origin=is_origin,
                           tree_died=False,
                           read_as_newick_string=True,
                           at_dict=dict(),
                           clado_at_dict=dict())
    
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