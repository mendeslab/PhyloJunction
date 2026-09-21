from phylojunction.data.trait import DiscreteTrait
import copy
import os.path
import typing as ty

# pj imports
import phylojunction.utility.exception_classes as ec
import phylojunction.functionality.stoch_map as pjsmap
import phylojunction.pgm.pgm as pgm
import phylojunction.data.tree as pjtr
import phylojunction.functionality.biogeo as pjbio
import phylojunction.calculation.discrete_sse as sseobj  # type: ignore

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


def make_mapped_ann_tree(det_fn_name: str,
                         det_fn_param_dict:
                         ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) -> \
                            pjtr.AnnotatedTree:

    tip_attr_file_path = str()
    maps_file_path = str()
    mapping_function = str()
    mapped_ann_tr_list: ty.List[pjtr.AnnotatedTree] = list()
    n_regions: int = 1
    geosse: bool = False

    # val is a list of strings or nodes
    for arg, val in det_fn_param_dict.items():
        if not val:
            raise ec.ParseMissingArgumentError(arg)

        if arg == "tree":
            if not isinstance(val[0], pgm.NodeDAG):
                raise(ec.ParseInvalidArgumentError(
                    "tree",
                    "'" + val[0].node_name + "'",
                    "Argument must be a DAG node holding a phylogenetic tree."
                ))

            # making deep copy of trees so their members can be
            # updated and the tree can be "scaled" independently
            # of the original unscaled trees
            for ann_tr in val[0].value:
                mapped_ann_tr_list.append(copy.deepcopy(ann_tr))

        elif arg == "tip_attr_file_path":
            if isinstance(val[0], pgm.NodeDAG):
                raise ec.ParseInvalidArgumentError(
                    "tip_attr_file_path",
                    "'" + val[0].node_name + "'",
                    "Argument must be a string containing a file path.")

            tip_attr_file_path = val[0].replace('"', "")

            if not os.path.isfile(tip_attr_file_path):
                raise ec.PJIOFileDoesNotExistError(
                    "map_attr",
                    tip_attr_file_path
                )

        elif arg == "maps_file_path":
            if isinstance(val[0], pgm.NodeDAG):
                raise ec.ParseInvalidArgumentError(
                    "maps_file_path",
                    "'" + val[0].node_name + "'",
                    "Argument must be a string containing a file path.")

            maps_file_path = val[0].replace('"', "")

            if not os.path.isfile(maps_file_path):
                raise ec.PJIOFileDoesNotExistError(
                    "map_attr",
                    maps_file_path
                )

        elif arg == "fun":
            if isinstance(val[0], pgm.NodeDAG):
                raise(ec.ParseInvalidArgumentError(
                    "fun",
                    "'" + val[0].node_name + "'",
                    "Argument must be one of the following strings: \"smap\"."
                ))

            elif isinstance(val[0], str) and \
                    val[0] not in ('"smap"'):
                raise ec.ParseInvalidArgumentError(
                    "fun",
                    "'" + val[0] + "'",
                    "Argument must be one of the following strings: \"smap\"."
                )

            mapping_function = val[0]

        elif arg == "n_regions":
            if isinstance(val[0], pgm.NodeDAG):
                # need to declare cast_val separately so mypy won't complain
                cast_val1: ty.List[pgm.NodeDAG] = \
                    ty.cast(ty.List[pgm.NodeDAG], val)

                if len(cast_val1) > 1:
                    raise ec.ParseRequireIntegerError("map_attr", "n_regions")

                n_regions = int(pgm.extract_value_from_dagnodes(cast_val1)[0])

            # val is a list of strings
            elif isinstance(val[0], str):
                cast_val2: ty.List[str] = ty.cast(ty.List[str], val)

                if len(cast_val2) > 1:
                    raise ec.ParseRequireIntegerError("map_attr", "n_regions")

                n_regions = int(int(cast_val2[0]))

        elif arg == "geo":
            if isinstance(val[0], pgm.NodeDAG):
                raise(ec.ParseInvalidArgumentError(
                    "fun",
                    "'" + val[0].node_name + "'",
                    "Argument must be \"true\" or \"false\"."
                ))

            elif isinstance(val[0], str) and \
                    val[0] not in ('"true"', '"false"'):
                raise (ec.ParseInvalidArgumentError(
                    "fun",
                    "'" + val[0] + "'",
                    "Argument must be \"true\" or \"false\"."
                ))

            geosse = True if val[0] == '"true"' else False

    ################################
    # Setting up stoch map parsing #
    ################################

    if mapping_function == '"smap"':
        n_states_per_char = 2 if geosse else -1
        state2bit_lookup = pjbio.State2BitLookup(n_regions, n_states_per_char, geosse=geosse)

        # updating the number of states in trees
        for ann_tr in mapped_ann_tr_list:
            ann_tr.trait = DiscreteTrait(state2bit_lookup.n_states, name=ann_tr.trait.name)

        # side-effect:
        # (i) annotated trees in mapped_ann_tr_list can now be plotted
        #     because they have had their members updated (at_dict and
        #     node_attr_dict)
        stoch_mapcoll = \
            pjsmap.StochMapsOnTreeCollection(maps_file_path,
                                             mapped_ann_tr_list,
                                             state2bit_lookup,
                                             node_states_file_path=tip_attr_file_path,
                                             stoch_map_attr_name=None)

        return mapped_ann_tr_list
