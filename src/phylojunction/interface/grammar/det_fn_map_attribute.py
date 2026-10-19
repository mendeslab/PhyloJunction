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


# File reading returns independent history values; map_attr applies one of them to one tree.
def read_maps(file_path, tip_attr_file_path='', n_regions=1, geo=False):
    from phylojunction.interface.grammar.constructors import flag, integer
    geosse, n_regions = flag(geo), integer(n_regions, 1)
    state2bit_lookup = pjbio.State2BitLookup(n_regions, 2 if geosse else -1, geosse=geosse)
    return pjsmap.read_histories(file_path, tip_attr_file_path, state2bit_lookup)


# One mapping operation always returns a copy; shared scalar trees remain reusable across
# histories. Collection pairing is performed by the ordinary language call dispatcher.
def map_attr(tree, history, attr_name=None):
    if not isinstance(tree, pjtr.AnnotatedTree) or not isinstance(history, pjsmap.StochasticHistory):
        raise TypeError("map_attr requires one tree and one history from read_maps.")
    mapped = copy.deepcopy(tree)
    mapped.trait = DiscreteTrait(history.state2bit_lookup.n_states, name=mapped.trait.name)
    history.apply(mapped, attr_name)
    return mapped
