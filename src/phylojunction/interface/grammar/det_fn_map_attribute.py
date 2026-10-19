from phylojunction.data.trait import DiscreteTrait
import copy
import os.path
import typing as ty

# pj imports
import phylojunction.utility.exception_classes as ec
import phylojunction.functionality.stoch_map as pjsmap
import phylojunction.data.tree as pjtr
import phylojunction.functionality.biogeo as pjbio
import phylojunction.calculation.discrete_sse as sseobj  # type: ignore

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


# Map copied tree values while preserving collection shape. A scalar tree stays scalar;
# a read_tree list stays a list, including when it contains only one tree.
def map_attr(tree, fun, tip_attr_file_path='', maps_file_path='', attr_name=None,
             n_regions=1, geo=False):
    from phylojunction.interface.grammar.constructors import flag, integer
    geosse, n_regions = flag(geo), integer(n_regions, 1)
    if fun != 'smap':
        raise ValueError("map_attr currently supports fun='smap'.")
    for path in (tip_attr_file_path, maps_file_path):
        if path and not os.path.isfile(path):
            raise ec.PJIOFileDoesNotExistError('map_attr', path)
    mapped = copy.deepcopy(tree)
    mapped_ann_tr_list = []

    # Gather tree leaves for the existing collection mapper without flattening the result.
    def collect(value):
        if isinstance(value, list):
            for item in value:
                collect(item)
        elif isinstance(value, pjtr.AnnotatedTree):
            mapped_ann_tr_list.append(value)
        else:
            raise TypeError("map_attr requires a tree or a nested list of trees.")

    collect(mapped)
    ################################
    # Setting up stoch map parsing #
    ################################

    if fun == 'smap':
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
                                             stoch_map_attr_name=attr_name)

        return mapped
