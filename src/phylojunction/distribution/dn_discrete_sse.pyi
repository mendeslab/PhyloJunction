import dendropy as dp # type: ignore
import typing as ty

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.data.attribute_transition as pjat
from phylojunction.data.tree import AnnotatedTree
from phylojunction.data.sampled_ancestor \
    import SampledAncestor  # type: ignore
from phylojunction.data.attribute_transition \
    import AttributeTransition  # type: ignore

class DnSSE(pgm.DistrForSampling):
    DN_NAME = "DnSSE"
    n_sim: int
    n_repl: int
    with_origin: bool
    root_is_born: bool
    start_states: ty.List[int]
    seed_age: ty.Optional[float]
    condition_on_speciation: bool
    condition_on_survival: bool
    condition_on_obs_both_sides_root: bool
    stop: str
    stop_val: ty.List[float]
    min_rec_taxa: int
    max_rec_taxa: int
    abort_at_living_count: int
    sse_stash: sseobj.SSEStash
    events: sseobj.MacroevolEventHandler
    state_count: int
    n_time_slices: int
    slice_t_ends: ty.List[ty.Optional[float]]
    prob_handler: sseobj.DiscreteStateDependentProbabilityHandler
    epsilon: float
    runtime_limit: int
    max_n_failed_attempts: int
    rng_seed: int
    debug: bool
    info: bool
    def __init__(self,
                 sse_stash: sseobj.SSEStash,
                 n: int = 1,
                 n_replicates: int = 1,
                 origin: bool = True,
                 start_states_list: ty.List[int] = [],
                 stop: str = "",
                 stop_value: ty.List[float] = [],
                 condition_on_speciation: bool = False,
                 condition_on_survival: bool = True,
                 condition_on_obs_both_sides_root: bool = False,
                 min_rec_taxa: int = 0,
                 max_rec_taxa: int = int(1e12),
                 abort_at_alive_count: int = int(1e12),
                 epsilon: float = 1e-12,
                 runtime_limit: int = 300,
                 max_n_failed_attempts: int = 100,
                 rng_seed: ty.Optional[int] = None,
                 debug: ty.Optional[bool] = False,
                 info: ty.Optional[bool] = False) -> None: ...
    def init_check_vectorize_sample_size(self) -> None: ...
    def _check_input_health(self) -> None: ...
    def _germinate_tree(self,
                        a_start_state: int,
                        with_origin: bool,
                        state_representation_dict: ty.Dict[int, ty.Set[str]],
                        untargetable_node_set: ty.Set[str]) -> ty.Tuple[dp.Tree, dp.Node, int, bool]: ...
    def _update_sa_lineage_dict(self,
                                a_time: float,
                                sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]],
                                sa_lineage_node_labels: ty.List[str],
                                debug: bool=False) -> None: ...
    def _execute_birth(self,
                       tr_namespace: dp.TaxonNamespace,
                       chosen_node: dp.Node,
                       state_representation_dict: ty.Dict[int, ty.Set[str]],
                       sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]],
                       state_transition_dict: ty.Dict[str, ty.List[AttributeTransition]],
                       clado_state_transition_dict: ty.Dict[str, AttributeTransition],
                       untargetable_node_set: ty.Set[str],
                       cumulative_node_count: int,
                       node_index: int,
                       sse_birth_rate_object: sseobj.DiscreteStateDependentRate,
                       event_t: float,
                       debug=False) -> ty.Tuple[dp.Node, int]: ...
    def _execute_death(self,
                       tr_namespace: dp.TaxonNamespace,
                       chosen_node: dp.Node,
                       state_representation_dict: ty.Dict[int, ty.Set[str]],
                       sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]],
                       untargetable_node_set: ty.Set[str],
                       event_t: float,
                       debug=False) -> dp.Node: ...
    def _execute_anatrans(self,
                          tr_namespace: dp.TaxonNamespace,
                          chosen_node: dp.Node,
                          state_representation_dict: ty.Dict[int, ty.Set[str]],
                          state_transition_dict: ty.Dict[str, ty.List[AttributeTransition]],
                          untargetable_node_set: ty.Set[str],
                          sse_anatrans_rate_object: sseobj.DiscreteStateDependentRate,
                          event_t: float,
                          debug: bool = False) -> None: ...
    def _execute_sample_ancestor(self,
                                 tr_namespace: dp.TaxonNamespace,
                                 chosen_node: dp.Node,
                                 state_representation_dict: ty.Dict[int, ty.Set[str]],
                                 sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]],
                                 untargetable_node_set: ty.Set[str],
                                 cumulative_sa_count: int,
                                 node_index: int,
                                 event_t: float,
                                 debug: bool = False) -> int: ...
    def _execute_event(self,
                       tr_namespace,
                       sse_rate_object: sseobj.DiscreteStateDependentRate,
                       chosen_node: dp.Node,
                       state_representation_dict: ty.Dict[int, ty.Set[str]],
                       state_transition_dict: ty.Dict[str, ty.List[AttributeTransition]],
                       clado_state_transition_dict: ty.Dict[str, AttributeTransition],
                       sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]],
                       untargetable_node_set: ty.Set[str],
                       cumulative_node_count: int,
                       cumulative_sa_count: int,
                       node_index: int,
                       event_t: float,
                       debug: bool = False) -> ty.Tuple[dp.Node, int, int]: ...
    def _annotate_sampled(self,
                          a_time: float,
                          living_nodes: ty.List[dp.Node],
                          sample_idx: int,
                          single_node: ty.Optional[dp.Node] = None) -> None: ...
    def _get_next_event_time(total_rate: float) -> float: ...
    def simulate(self, a_start_state: int, value_idx: int = ..., a_seed: ty.Optional[int] = ...) -> AnnotatedTree: ...
    def draw(self) -> AnnotatedTree: ...
    def generate(self) -> ty.List[AnnotatedTree]: ...
    def get_rev_inference_spec_info(self) -> ty.List[str]: ...
    def is_tr_ok(self, ann_tr: AnnotatedTree) -> bool: ...
