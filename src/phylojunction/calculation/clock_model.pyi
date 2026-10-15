import typing as ty
from abc import abstractmethod

import phylojunction.pgm.pgm as pgm
import phylojunction.data.tree as pjtr

class ClockModel(pgm.DeterministicFn):
    attr_name: str
    tree_dag_node: pgm.NodeDAG
    n_samples: int
    n_repl: int
    @property
    @abstractmethod
    def CLOCK_MODEL_NAME(self) -> str: ...
    @property
    def DET_FN_NAME(self) -> str: ...
    def __init__(self, *, tree_dag_node: pgm.NodeDAG, attr_name: str) -> None: ...
    @abstractmethod
    def populate_branch_rates(self, ann_tr: pjtr.AnnotatedTree, tree_idx: int) -> None: ...
    def generate(self) -> ty.List[pjtr.AnnotatedTree]: ...
