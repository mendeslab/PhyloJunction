import typing as ty

import phylojunction.pgm.pgm as pgm
import phylojunction.data.tree as pjtr
import phylojunction.calculation.clock_model as pjclock

def populate_gamma_spike_mol_rates_on_rec_tr_branches(ann_tr: pjtr.AnnotatedTree, clock_rate: float, spike_mean: float, spike_shape: float, br_specific_rel_rates: ty.List[float], strict_clock_when_no_stubs: bool = ...) -> None: ...

class GammaSpikeClockModel(pjclock.ClockModel):
    CLOCK_MODEL_NAME: str
    clock_rate: float
    spike_mean: float
    spike_shape: float
    strict_clock_when_no_stubs: bool
    br_specific_rel_rate_dag_node: pgm.NodeDAG
    def __init__(self, *, tree_dag_node: pgm.NodeDAG, attr_name: str, clock_rate: float, spike_mean: float, spike_shape: float, br_specific_rel_rate_dag_node: pgm.NodeDAG, strict_clock_when_no_stubs: bool = ...) -> None: ...
    def populate_branch_rates(self, ann_tr: pjtr.AnnotatedTree, tree_idx: int) -> None: ...
    def init_check_vectorize_sample_size(self, param_list: ty.List[ty.Any] = ...) -> None: ...
    def get_rev_inference_spec_info(self) -> ty.List[str]: ...
