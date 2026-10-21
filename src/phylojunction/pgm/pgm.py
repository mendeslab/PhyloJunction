from __future__ import annotations
from abc import ABC, abstractmethod
from numbers import Integral
import typing as ty
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import phylojunction.data.tree as pjtr
from phylojunction.utility.exception_classes import ScriptEvaluationError

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class DirectedAcyclicGraph:
    """Shared, ordered definitions of a model, without replicate-specific values.

    name_node_dict points to the latest definition of each source name. Every
    definition has a stable position in definitions, so reassignment preserves
    the bindings of earlier statements. ExecutionResults owns sampled values.
    """

    # A graph is compiled before it is executed; compiling consumes no random numbers.
    def __init__(self):
        self.definitions = []
        self.name_node_dict = {}
        self.n_nodes = 0

    # Bind variable references when a statement is appended, preserving sequential
    # reassignment: later definitions cannot alter dependencies of earlier statements.
    def add_definition(self, source, source_name="<interactive>", line_number=1):
        from phylojunction.interface.cmdbox.expressions import parse_statement, variable_names
        name = None
        try:
            statement = parse_statement(source)
            if statement is None:
                return
            name, operator, expression = statement
            # Clamping is observation metadata, not a parameter of the distribution object.
            import ast
            clamped = False
            if operator == '~' and isinstance(expression, ast.Call):
                clamp = next((item for item in expression.keywords if item.arg == 'clamp'), None)
                if clamp is not None:
                    from phylojunction.interface.grammar.constructors import flag
                    if isinstance(clamp.value, ast.Name) and clamp.value.id in ("true", "false"):
                        clamp.value = ast.Constant(clamp.value.id == "true")
                    if not isinstance(clamp.value, ast.Constant):
                        raise ValueError("clamp must be a literal true/false flag.")
                    clamped = flag(clamp.value.value)
                    expression.keywords.remove(clamp)
            # Clamping belongs to the entire sampling statement, never to inner constructors.
            # Reject nested flags even inside empty comprehensions that will not be evaluated.
            if any(isinstance(item, ast.Call) and any(kw.arg == 'clamp' for kw in item.keywords)
                   for item in ast.walk(expression)):
                raise ValueError("clamp is allowed only on the outer call of a sampling statement.")
            names = list(dict.fromkeys(variable_names(expression)))
            missing = [name for name in names if name not in self.name_node_dict]
            if missing:
                raise ValueError(f"Undefined variable(s): {', '.join(missing)}")
            parents = [self.name_node_dict[name] for name in names]
            node_class = DeterministicNodeDAG if operator == ':=' else StochasticNodeDAG
            node = node_class(name, parent_nodes=parents)
            node.definition_id = len(self.definitions)
            node.call_order_idx = node.definition_id
            node.expression = expression
            node.operator = operator
            node.source = source
            node.source_name = source_name
            node.line_number = line_number
            node.bindings = {name: self.name_node_dict[name].definition_id for name in names}
            node.is_sampled = operator == '~'
            node.is_clamped = clamped
            self.definitions.append(node)
            self.name_node_dict[name] = node
            self.n_nodes = len(self.definitions)
        except ScriptEvaluationError:
            raise
        except Exception as error:
            raise ScriptEvaluationError(str(error), source=source, source_name=source_name,
                                        line_number=line_number, variable=name) from error

    def execute(self, runs=1, random_seed=None):
        from phylojunction.pgm.execution import ExecutionResults
        return ExecutionResults(self, runs, random_seed)

    def get_node_dag_by_name(self, name):
        return self.name_node_dict[name]

    def get_sorted_node_dag_list(self):
        return list(self.definitions)


class DistrForSampling(ABC):
    """Class for randomly generating (i.e., sampling) values.

    A distribution describes one draw. Model execution owns repetition;
    IID distributions explicitly describe vector-valued draws.
    """

    # Validate transient draw counts independently of any model replicate count.
    @staticmethod
    def validate_draw_count(count):
        if isinstance(count, bool) or not isinstance(count, Integral) or count < 0:
            raise ValueError("The number of draws must be a nonnegative integer.")

    def draw(self):
        raise NotImplementedError("This distribution does not implement single-value draws yet.")

    # Scalar distributions override this for bulk sampling; tree and nested-vector
    # distributions repeat whole draws so each tree has independent resource limits.
    def draw_many(self, count):
        self.validate_draw_count(count)
        return [self.draw() for _ in range(count)]

    # python3's way of "declaring" an (required) abstract attribute
    @property
    @abstractmethod
    def DN_NAME(self) -> str:
        raise NotImplementedError


class ConstantFn(ABC):
    """A deterministic reader/constructor, independent of replicate counts."""

    @abstractmethod
    def generate(self):
        raise NotImplementedError


class NodeDAG:
    """A named definition and its bound dependencies, shared by all replicates.

    expression is the validated syntax tree; bindings map source names to earlier
    definition IDs. parent_nd_list retains graph edges for graph consumers. Values
    and sampled-from distribution objects are kept in ExecutionResults instead.
    """

    # Definitions have identity independent of their source name, which may be reassigned.
    def __init__(self, node_name, parent_nodes=None):
        self.node_name = node_name
        self.parent_nd_list = [] if parent_nodes is None else list(parent_nodes)
        self.is_sampled = False
        self.is_deterministic = False
        self.is_clamped = False
        self.definition_id = None
        self.call_order_idx = None
        self.expression = None
        self.bindings = {}
        self.operator = None
        self.source = None
        self.source_name = "<interactive>"
        self.line_number = 1

    def __str__(self):
        return self.source or self.node_name


class StochasticNodeDAG(NodeDAG):
    """A fixed assignment or a stochastic sampling definition."""


class DeterministicNodeDAG(NodeDAG):
    """A deterministic expression, including a standalone distribution object."""

    def __init__(self, node_name, parent_nodes=None):
        super().__init__(node_name, parent_nodes)
        self.is_deterministic = True


def plot_node_histogram(axes: plt.Axes,
                        values_list: ty.List[float],
                        sample_idx: ty.Optional[int] = None,
                        repl_size: int = 1) -> None:

    values_list_to_plot = list()
    # if not sample_idx == None:
    if isinstance(sample_idx, int):
        start = sample_idx * repl_size
        end = start + repl_size
        values_list_to_plot = [float(v) for v in values_list[start:end]]

    else:
        values_list_to_plot = [float(v) for v in values_list]

    # figure canvas was created outside main loop in GUI
    pjtr.clear_trait_colorbar(axes)
    axes.cla()
    counts, bins, _ = \
        axes.hist(values_list_to_plot,
                  histtype='stepfilled',
                  color="steelblue",
                  edgecolor="none",
                  alpha=0.3)

    axes.spines['left'].set_visible(True)
    axes.spines['bottom'].set_visible(True)
    # axes.patch.set_alpha(0.0)
    # bins = axes.get_xticks().tolist()
    axes.xaxis.set_major_locator(mticker.FixedLocator(bins))
    # label_format = '{:,.0f}'
    # axes.set_xticklabels([label_format.format(x) for x in bins])
    plt.xticks(bins)
    plt.yticks(np.linspace(0, max(counts.tolist()), 5))
    axes.axes.xaxis.set_ticklabels([np.round(i, decimals=2) for i in bins])


def plot_blank(axes: plt.Axes) -> None:
    pjtr.clear_trait_colorbar(axes)
    axes.cla()
    axes.patch.set_alpha(0.0)
    axes.xaxis.set_ticks([])
    axes.yaxis.set_ticks([])
    axes.spines['left'].set_visible(False)
    axes.spines['bottom'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)


