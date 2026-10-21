from __future__ import annotations
from abc import ABC, abstractmethod
from numbers import Integral
import typing as ty
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import phylojunction.data.tree as pjtr

class DirectedAcyclicGraph:
    definitions: list[NodeDAG]
    name_node_dict: dict[str, NodeDAG]
    n_nodes: int

    def __init__(self):
        ...

    def add_definition(self, source, source_name="<interactive>", line_number=1):
        ...

    def execute(self, runs=1, random_seed=None):
        ...

    def get_node_dag_by_name(self, name):
        ...

    def get_sorted_node_dag_list(self):
        ...

class DistrForSampling(ABC):

    @staticmethod
    def validate_draw_count(count):
        ...

    def draw(self):
        ...

    def draw_many(self, count):
        ...

    @property
    @abstractmethod
    def DN_NAME(self) -> str:
        ...

class ConstantFn(ABC):

    @abstractmethod
    def generate(self):
        ...

class NodeDAG:
    node_name: str
    parent_nd_list: list[NodeDAG]
    definition_id: int
    call_order_idx: int
    bindings: dict[str, int]
    expression: ty.Any
    operator: str
    source: str
    source_name: str
    line_number: int
    is_sampled: bool
    is_deterministic: bool
    is_clamped: bool

    def __init__(self, node_name, parent_nodes=None):
        ...

    def __str__(self):
        ...

class StochasticNodeDAG(NodeDAG):
    ...

class DeterministicNodeDAG(NodeDAG):

    def __init__(self, node_name, parent_nodes=None):
        ...

def plot_node_histogram(axes: plt.Axes, values_list: ty.List[float], sample_idx: ty.Optional[int]=None, repl_size: int=1) -> None:
    ...

def plot_blank(axes: plt.Axes) -> None:
    ...
