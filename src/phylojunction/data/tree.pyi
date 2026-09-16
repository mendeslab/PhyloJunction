import dendropy as dp # type: ignore
import matplotlib.pyplot as plt # type: ignore
import typing as ty
from matplotlib import colors
from matplotlib.figure import Figure as Figure # type: ignore

import phylojunction.data.sampled_ancestor as pjsa
import phylojunction.data.attribute_transition as pjat

class AnnotatedTree(dp.Tree):
    tree: dp.Tree
    tree_reconstructed: dp.Tree
    origin_node: ty.Optional[dp.Node]
    root_node: ty.Optional[dp.Node]
    rec_tr_root_node: ty.Optional[dp.Node]
    rec_tr_root_node_label: ty.Optional[str]
    brosc_node: ty.Optional[dp.Node]
    with_origin: bool
    tree_read_as_newick: bool
    condition_on_obs_both_sides_root: bool
    tree_died: ty.Optional[bool]
    tree_invalid: ty.Optional[bool]
    seed_age: float
    max_age: ty.Optional[float]
    origin_age: ty.Optional[float]
    origin_edge_length: float
    root_age: float
    rec_tr_root_age: ty.Optional[float]
    node_heights_dict: ty.Dict[str, float]
    rec_tr_node_heights_dict: ty.Dict[str, float]
    node_ages_dict: ty.Dict[str, float]
    rec_tr_node_ages_dict: ty.Dict[str, float]
    slice_t_ends: ty.Optional[ty.List[float]]
    slice_age_ends: ty.Optional[ty.List[float]]
    state_count: int
    state_count_dict: ty.Dict[int, int]
    extant_terminal_state_count_dict: ty.Dict[int, int]
    extant_terminal_sampled_state_count_dict: ty.Dict[int, int]
    extinct_terminal_state_count_dict: ty.Dict[int, int]
    node_attr_dict: ty.Dict[str, ty.Dict[str, ty.Any]]
    n_extant_terminal_nodes: int
    n_extinct_terminal_nodes: int
    n_extant_sampled_terminal_nodes: int
    n_extinct_sampled_terminal_nodes: int
    n_sa_nodes: int
    extant_terminal_nodes_labels: ty.Tuple[str, ...] # ... = all items of same kind
    extinct_terminal_nodes_labels: ty.Tuple[str, ...]
    extant_sampled_terminal_nodes_labels: ty.Tuple[str, ...]
    sa_obs_nodes_labels: ty.Tuple[str]
    sa_lineage_dict: ty.Optional[ty.Dict[str, ty.List[pjsa.SampledAncestor]]]
    at_dict: ty.Optional[ty.Dict[str, ty.List[pjat.AttributeTransition]]]
    rec_tr_at_dict: ty.Optional[ty.Dict[str, ty.List[pjat.AttributeTransition]]]
    have_updated_rec_tr_at_dict: bool
    clado_at_dict: ty.Optional[ty.Dict[str, pjat.AttributeTransition]]
    have_updated_rec_tr_at_dict: bool
    epsilon: float

    def __init__(self,
                 a_tree: dp.Tree,
                 total_state_count: int,
                 start_at_origin: bool = False,
                 condition_on_obs_both_sides_root: bool = False,
                 max_age: ty.Optional[float] = None,
                 slice_t_ends: ty.Optional[ty.List[float]] = [],
                 slice_age_ends: ty.Optional[ty.List[float]] = None,
                 sa_lineage_dict: ty.Optional[ty.Dict[str, ty.List[pjsa.SampledAncestor]]] = None,
                 at_dict: ty.Optional[ty.Dict[str, ty.List[pjat.AttributeTransition]]] = None,
                 clado_at_dict: ty.Optional[ty.Dict[str, pjat.AttributeTransition]] = None,
                 tree_died: ty.Optional[bool] = None,
                 tree_invalid: ty.Optional[bool] = None,
                 read_as_newick_string: bool = False,
                 epsilon: float = 1e-12) -> None: ...
    def _check_input_health(self) -> None: ...
    def _has_tree_died(self, origin_or_root_node: dp.Node) -> bool: ...
    def _recursively_find_node_time(self, a_node: dp.Node, reference_time: float = 0.0) -> float: ...
    def _init_and_update_origin_root_members(self) -> None: ...
    def _count_terminal_nodes(self) -> None: ...
    def _count_sampled_ancestors(self) -> None: ...
    def _count_node_states(self) -> None: ...
    def _populate_node_age_height_dicts(self, unit_branch_lengths: bool = False, do_reconstructed_tree: bool = False) -> None: ...
    def _prepare_taxon_namespace_for_nexus_printing(self) -> None: ...
    def _add_at_ages(self) -> None: ...
    def is_extant_or_sa_on_both_sides_complete_tr_root(self, a_node: dp.Node) -> bool: ...
    def update_rec_tr_at_dict(self, rec_tree_root_node: dp.Node) -> None: ...

    def update_rec_tr_sa_lineage_dict(self) -> None: ...
    def extract_reconstructed_tree(self, plotting_overhead: bool = False, require_obs_both_sides: ty.Optional[bool] = None) -> dp.Tree: ...
    def populate_nd_attr_dict(self, attrs_of_interest_list: ty.List[str], attr_dict_added_separately_from_tree: bool = False) -> None: ...
    def __str__(self) -> str: ...
    def rec_str(self) -> str: ...
    def plot_node(self, axes: plt.Axes, node_attr: str = "state", **kwargs) -> None: ...
    def get_stats_dict(self) -> ty.Dict[str, ty.Union[int, float]]: ...
    def _get_taxon_states_dict(self) -> ty.Dict[str, int]: ...
    def get_taxon_states_str(self, nexus: bool = False, bit_format: bool = False) -> str: ...

# plotting tree functions
def clear_trait_colorbar(axes: plt.Axes) -> None: ...
def get_node_name(nd: dp.Node) -> str: ...
def get_x_coord_from_nd_heights(ann_tr: AnnotatedTree, use_age: bool = False, unit_branch_lengths: bool = False) -> ty.Dict[str, float]: ...
def get_y_coord_from_n_obs_nodes(ann_tr: AnnotatedTree, start_at_origin: bool = False, sa_along_branches: bool = True, draw_reconstructed: bool = False) -> ty.Dict[str, float]: ...
def plot_ann_tree(ann_tr: AnnotatedTree, axes: plt.Axes, use_age: bool = False, start_at_origin: bool = False, attr_of_interest: str = "state", sa_along_branches: bool = True, draw_reconstructed: bool = False) -> ty.Tuple[ty.Dict[str, float]]: ...
def get_color_map(n_states: int) -> ty.Dict[int, str]: ...
def pj_get_mrca_obs_terminals(a_node: dp.Node, nd_label_list: ty.List[str]) -> dp.Node: ...
