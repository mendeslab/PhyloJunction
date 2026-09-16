from __future__ import annotations
import typing as ty
import numpy as np
import random
import matplotlib.pyplot as plt  # type: ignore
import matplotlib.ticker as mticker  # type: ignore
import statistics as stat  # type: ignore
import pandas as pd  # type: ignore
from tabulate import tabulate  # type: ignore
from abc import ABC, abstractmethod
from phylojunction.calculation.continuous_sse import LogisticRate

# pj imports
import phylojunction.utility.exception_classes as ec
import phylojunction.data.tree as pjtr

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


# code for @abstract_attribute
R = ty.TypeVar('R')
def abstract_attribute(obj: ty.Callable[[ty.Any], R] = None) -> R:

    class DummyAttribute:
        pass

    _obj = ty.cast(ty.Any, obj)

    if obj is None:
        _obj = DummyAttribute()

    _obj.__is_abstract_attribute__ = True

    return ty.cast(R, _obj)


class DirectedAcyclicGraph():
    """Directed acyclic graph (DAG) class defining the model.
    
    Parameters:
        node_val_dict (dict): Dictionary with keys being DAG nodes, and
            values being the values stored in the nodes.
        name_node_dict (dict): Dictionary with keys being DAG node
            names, and the DAG nodes themselves as values.
        n_nodes (int): Total number of DAG nodes in the graph.
        sample_size (int): How many samples (simulations) to be either
            drawn or specified in every DAG node.
    """

    node_val_dict: ty.Dict[NodeDAG, ty.Any]
    name_node_dict: ty.Dict[str, NodeDAG]
    n_nodes: int
    sample_size: int
    _random_seed: int

    def __init__(self) -> None:
        # keys are proper DAG nodes, values are their values
        self.node_val_dict = dict()
        # keys are DAG node names, vals are NodeDAG instances
        self.name_node_dict = dict()
        self.n_nodes = 0
        self.sample_size = 0  # how many simulations will be run
        self._random_seed = None

    @property
    def random_seed(self) -> int:
        return self._random_seed
    
    @random_seed.setter
    def random_seed(self, a_seed) -> None:
        # handle seed if not None, not empty string, and is integer
        if a_seed is not None and isinstance(a_seed, int):
            self._random_seed = a_seed
            
            # now we execute the seed (for two random number generators)!
            np.random.seed(seed=a_seed)
            random.seed(a_seed)

        else:
            raise ec.DAGCannotInitialize(
                "seed was 'None'. It must be an integer.")
    

    def add_node(self, node_dag: NodeDAG) -> None:
        # check that nodes carry the right number of values
        # (the number of simulations)
        if isinstance(node_dag, StochasticNodeDAG):
            #############
            # Important #
            #############

            # note how only sampled nodes have any business in setting
            # the sample size of a DAG object; this means we let the users
            # fool around with nodes with assigned (fixed, clamped) values
            # through '<-'
            if node_dag.is_sampled:
                # if the pgm's sample size is still 0,
                # or if we started off with a scalar node but then
                # added sampled node, we update the DAG's sample size
                if not self.sample_size or \
                        (self.sample_size == 1 and node_dag._sample_size > 1):
                    self.sample_size = node_dag._sample_size

                # if the number of values in a node is 1, it gets vectorized,
                # so this is allowed; but if the node is sampled and the number
                # of values is > 1 and < than that of other nodes, we have a
                # problem
                elif self.sample_size != node_dag._sample_size and \
                        node_dag._sample_size > 1:

                    raise ec.DAGCannotAddNodeError(
                        node_dag.node_name,
                        ("Specified number of simulations was both "
                         "> 1 and different from that in a previous "
                         "command."))

        # a DAG node has been specified, but it already was in the
        # DAG, meaning we need to update the DAG
        if node_dag in self.node_val_dict:
            # removes existing node with that node_name
            call_order_idx = node_dag.call_order_idx
            del self.node_val_dict[node_dag]
            node_dag.call_order_idx = call_order_idx

        # new node (random variable), add to call order
        else:
            self.n_nodes += 1
            node_dag.call_order_idx = self.n_nodes

        # if node is not deterministic and has a value and node_name
        try:
            self.node_val_dict[node_dag] = node_dag.value

        # have to double-check this can ever be a problem
        except AttributeError:
            print(("\n\n\nSomehow we could not grab a node's value."
                  "Come back to this line of code!"))
            self.node_val_dict[node_dag] = None

        self.name_node_dict[node_dag.node_name] = node_dag

    def get_node_dag_by_name(self, node_name):
        if node_name in self.name_node_dict:
            return self.name_node_dict[node_name]

    def get_display_str_by_name(
            self,
            node_name,
            sample_idx=None,
            repl_size=1):

        if node_name in self.name_node_dict:
            # calls __str__() of NodeDAG
            return str(self.name_node_dict[node_name])

    def get_sorted_node_dag_list(self) -> ty.List[NodeDAG]:
        node_dag_list: ty.List[NodeDAG] = [node_dag for node_dag in self.node_val_dict]
        node_dag_list.sort()

        return node_dag_list


class ValueGenerator(ABC):
    @abstract_attribute
    def n_samples(self):
        pass

    @abstract_attribute
    def n_repl(self):
        pass

    @abstractmethod
    def __init__(self) -> None:
        pass

    @abstractmethod
    def generate(self) -> ty.List[ty.Any]:
        pass

    @abstractmethod
    def init_check_vectorize_sample_size(
        self,
        param_list: ty.List[ty.Any] = []) \
            -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]:
        """Check sample size against number of provided parameter values

        This is the function behind the vectorization functionality
        """
        pass

    @abstractmethod
    def get_rev_inference_spec_info(self) -> ty.List[str]:
        pass


class DistrForSampling(ValueGenerator):
    """Class for randomly generating (i.e., sampling) values.

    An object of this class is required by stochastic DAG nodes so
    values can be sampled.
    """

    # python3's way of "declaring" an (required) abstract attribute
    @property
    @abstractmethod
    def DN_NAME(self) -> str:
        raise NotImplementedError


class ConstantFn(ValueGenerator):
    """Class for deterministically generating values.
    
    An object is class whenever user input must be parsed and modified
    before it can be stored in a DAG node.
    """

    @property
    @abstractmethod
    def CT_FN_NAME(self):
        raise NotImplementedError


class NodeDAG(ABC):
    """Node class for building directed acyclic graphs (DAGs).
    
    This class has abstract methods, so it cannot be instantiated.
    It is derived by StochasticNodeDAG and DeterministicNodeDAG

    Parameters:
        node_name (str): Name of the node, e.g., 'lambda'.
        sample_size (int): How many samples (simulations) to be either
            drawn or specified for the node.
        replicate_size (int): How many replicates (for each sample) to
            be drawn or specified the node. This is the size of a
            'plate' in graphical notation. Defaults to 1.
        value (object): List of values associated to node. Defaults to
            None.
        call_order_idx (int): The order at which this node was added to
            DAG. If it is the first node of the DAG, for example, this
            value is 1. Defaults to None.
        sampled (bool): Flag specifying if what the 'value' attribute
        stores are stochastic samples from a distribution. Defaults to
            'False'.
        deterministic (bool): Flag specifying if what the 'value'
            attribute stores is the output of a deterministic function.
            Defaults to 'False'.
        clamped (bool): Flag specifying if what the 'value' attribute
            stores is observed (i.e., data). Defaults to 'False'.
        parent_nodes (NodeDAG): List of (parent) NodeDAG objects that
            are in the path between this node and the outermost layer
            in the DAG. Defaults to None.
    """

    node_name: str
    _value: ty.List[ty.Any]
    _sample_size: int
    _repl_size: int
    call_order_idx: int
    is_sampled: bool
    is_deterministic: bool
    is_clamped: bool
    parent_nd_list: ty.Optional[ty.List[NodeDAG]]

    # why decorate __init__:
    #
    # (i)   don't want to allow NodeDAG to be instantiated
    # (ii)  want to provide a "default" initializer
    # (iii) want to force all daughter classes to define their own
    # __init__(), and in there have a super().__init__() call
    @abstractmethod
    def __init__(self,
                 node_name: str,
                 sample_size: ty.Optional[int],
                 replicate_size: int = 1,
                 value: ty.Optional[ty.List[ty.Any]] = None,
                 call_order_idx: ty.Optional[int] = None,
                 sampled: bool = False,
                 deterministic: bool = False,
                 clamped: bool = False,
                 parent_nodes: ty.Optional[ty.List[NodeDAG]] = None):

        self.node_name = node_name
        self._value = value
        self._sample_size = sample_size
        self._repl_size = replicate_size
        self.call_order_idx = call_order_idx
        self.is_sampled = sampled
        self.is_deterministic = deterministic
        self.is_clamped = clamped
        self.parent_nd_list = parent_nodes
        # note that when the dag_obj adds this to its list of nodes,
        # value will be None (value is populated when we call .sample()),
        # and self.length will be = 1; we nonetheless add this call here
        # for completion (useful in debugging and testing)
        #
        # as we build the PGF through a script/gui, self.populate_length()
        # is in fact called from outside through method get_length()

        # total number of values
        # self.full_length = self.length * self.repl_size
        self.param_of = None

        if isinstance(self._value, (list, np.ndarray)) and not \
           isinstance(self._value[0], (pjtr.AnnotatedTree, LogisticRate)):
            self._flatten_and_extract_values()

    @property
    def value(self):
        return self._value
    
    @value.setter
    def value(self, a_value):
        self._value = a_value

    # no setter!
    @property
    def sample_size(self):
        return self._sample_size

    # no setter!
    @property
    def repl_size(self):
        return self._repl_size

    # TODO: maybe later make n_samples and n_repl properties
    # and make it so their setters can only raise Errors telling
    # callers that they can only be set upon initialization

    # side-effect updates self._value
    def _flatten_and_extract_values(self) -> None:
        """If value member is 2D-list, flatten it to 1D-list."""

        values: ty.List[ty.Any] = []

        # so mypy won't complain
        if isinstance(self._value, list):
            for v in self._value:
                # ._value is a 2D-list
                if not isinstance(v, (int, float, str, np.float64)):
                    values_inside_nodes = v.value

                    # so mypy won't complain
                    if isinstance(values_inside_nodes, list):
                        for val in values_inside_nodes:
                            if not isinstance(val,
                                              (int, float, str, np.float64)):
                                raise ec.VariableAssignmentError(self.node_name)

                            values.append(val)

                # each value in ._value is scalar
                else:
                    values.append(v)

        self._value = values

    # called by GUI
    def get_start2end_str(self,
                          start: int,
                          end: int,
                          repl_idx: int = 0,
                          is_tree: bool = False,
                          display_reconstructed: bool = False) -> str:
        """Get string representation of values within specific range."""     

        if isinstance(self._value, np.ndarray):
            self._value = \
                ", ".join(str(v) for v in self._value.tolist()[start:end])

        if isinstance(self._value, list):
            if len(self._value) >= 2:
                if not is_tree:
                    if isinstance(self._value[0],
                                  (int, float, str, np.float64)):
                        return ", ".join(
                            str(v) for v in self._value[start:end])

                    else:
                        return "\n".join(
                            str(v) for v in self._value[start:end])

                # tree, and repl_idx is always 0
                else:
                    if display_reconstructed:
                        tr = self._value[start + repl_idx]
                        return tr.rec_str()

                    # string representation of complete tree (__str__())
                    else:
                        return str(self._value[start + repl_idx])

            # single element in value
            else:
                return str(self._value[0])

        return str(self._value)

    def __hash__(self):
        return hash(self.node_name)

    def __eq__(self, other) -> bool:
        return other.node_name == self.node_name

    # stringify _all_ values
    def __str__(self) -> str:
        # try:
        if isinstance(self.value, list):
            return self.get_start2end_str(0, len(self.value))

        # if this blows up at some point, find out
        # what kind of error is possible and
        # reimplement the except making it non-bare
        # except:
        else:
            return str(self.value)

        # cosmetic return required by mypy
        return ""

    def __lt__(self, other) -> bool:
        return self.call_order_idx < other.call_order_idx

    def __len__(self) -> int:
        try:
            n_values = 1

            if isinstance(self._value, list):
                n_values = len(self._value)

            else:
                # could be instance of:
                #     SSE rate
                #     SSE prob
                if hasattr(self._value, "value") and \
                        isinstance(self._value.value, list):
                    n_values = len(self._value.value)

            n_repls = self._repl_size

            if n_values >= 1:
                if n_values % n_repls == 0:
                    return int(n_values / n_repls)

                else:
                    raise ec.ReplicateNumberError(node_name=self.node_name)

        # just being ultra-safe...
        # self.value has no len()
        # n_repls is None
        except TypeError:
            return 1

        # cosmetic return required by mypy
        return 1

    @abstractmethod
    def plot_node(self,
                  axes: plt.Axes,
                  sample_idx: int = 0,
                  repl_idx: int = 0,
                  repl_size: int = 1,
                  branch_attr: ty.Optional[str] = "state",
                  draw_reconstructed: ty.Optional[bool] = False) -> None:
        pass

    @abstractmethod
    def populate_operator_weight(self) -> None:
        pass

    @abstractmethod
    def get_node_stats_str(self,
                           start: int,
                           end: int,
                           repl_idx: int) -> str:
        stats_str = ""

        if isinstance(self.value, list):
            # could be (i) >= 1 samples with 1 replicate
            #          (ii) >= 1 samples with > 1 replicates
            # and we do (end - start) > 1 to make sure (i) doesn't go through
            if len(self.value) >= 2 and (end - start) > 1:
                if isinstance(self.value[0], (str, int, float, np.float64)):
                    # do stuff with the below
                    repl_values = [float(v) for v in self.value[start:end]]
                    repl_value_avg = stat.mean(repl_values)
                    repl_value_sd = stat.stdev(repl_values)

                    # so that it is very explicit what we are averaging over
                    avg_str, std_str = str(), str()

                    # all samples, all replicates
                    if (end - start) == len(self.value):
                        avg_str = "Sample average   = "
                        std_str = "Sample std. dev. = "

                    # one sample (all replicates)
                    else:
                        avg_str = "Replicate average   = "
                        std_str = "Replicate std. dev. = "

                    stats_str = \
                        avg_str + str(repl_value_avg) + "\n" + \
                        std_str + str(repl_value_sd)

                # values are objects (e.g., AnnotatedTree), not scalars
                else:
                    # do stuff with the below
                    repl_objs = self.value[start:end]
                    repl_all_stats_dict: ty.Dict[str, ty.List[float]] = dict()

                    # focal replicate given by repl_idx above
                    focal_repl_stats_dict: \
                        ty.Dict[str, ty.Union[int, float]] = dict()
                    # all replicates
                    repl_value_avg_dict: ty.Dict[str, float] = dict()
                    # all replicates
                    repl_value_sd_dict: ty.Dict[str, float] = dict()

                    # reading over all obj replicates and collecting
                    # their stats
                    for idx, repl_obj in enumerate(repl_objs):
                        obj_stat_dict = repl_obj.get_stats_dict()

                        for st, stat_v in obj_stat_dict.items():
                            float_stat_v: float = 0.0

                            if stat_v != "None":
                                try:
                                    float_stat_v = float(stat_v)

                                except ValueError:
                                    raise ec.NodeDAGNodeStatCantFloatError(
                                        self.node_name)

                            try:
                                repl_all_stats_dict[st].append(float_stat_v)

                            except KeyError:
                                repl_all_stats_dict[st] = [float_stat_v]

                            # getting focal replicate stats
                            if idx == repl_idx:
                                focal_repl_stats_dict = obj_stat_dict

                    # now getting means and sds
                    for st, stat_v_list in repl_all_stats_dict.items():
                        repl_value_avg_dict[st] = stat.mean(stat_v_list)
                        repl_value_sd_dict[st] = stat.stdev(stat_v_list)

                    # focal replicate
                    focal_obj_repl_value_df = \
                        pd.DataFrame(focal_repl_stats_dict.items())
                    # all replicates
                    obj_repl_value_avg_df = \
                        pd.DataFrame(repl_value_avg_dict.items())
                    # all replicates
                    obj_repl_value_sd_df = \
                        pd.DataFrame(repl_value_sd_dict.items())

                    # putting it all together
                    obj_repl_stats_df = \
                        pd.concat([focal_obj_repl_value_df,
                                   obj_repl_value_avg_df.iloc[:, 1],
                                   obj_repl_value_sd_df.iloc[:, 1]],
                                  ignore_index=True, axis=1)
                    repl_stats_str = tabulate(
                        obj_repl_stats_df,
                        headers=["Summary stat.",
                                 "This replicate",
                                 "Replicates avg.",
                                 "Replicates std. dev."],
                        tablefmt="plain",
                        showindex=False
                    ).lstrip()

                    stats_str = repl_stats_str

            # looking at single replicate (from any one sample)
            elif (end - start) == 1:
                if isinstance(self.value[0], (str, int, float, np.float64)):
                    return stats_str

                # if it's an object (e.g., AnnotatedTree, or SequenceAlignment)
                else:
                    obj_stats_dict = self.value[start:end][0].get_stats_dict()
                    obj_stats_df = pd.DataFrame(obj_stats_dict.items())
                    stats_str = tabulate(
                        obj_stats_df,
                        headers=["", ""],
                        tablefmt="plain",
                        showindex=False
                    ).lstrip()

        return stats_str


class StochasticNodeDAG(NodeDAG):
    """Derived DAG node that can generate own value.
    
    This class is also used when a value is being specified by a user,
    or parsed from it. In other words, constant nodes in the graph
    are still stochastic DAG nodes under the hood.

    Parameters:
        sampling_dn (DistrForSampling): Distribution object for randomly
            sampling values.
        constant_fn (ConstantFn): Function object for obtaining
            constant values modified from some user input.
        operator_weight (float): Weight given to MCMC move.
    """

    sampling_dn: ty.Optional[DistrForSampling]
    constant_fn: ty.Optional[ConstantFn]
    operator_weight: float

    def __init__(self,
                 node_name: str,
                 sample_size: int,
                 sampled_from: ty.Optional[DistrForSampling] = None,
                 returned_from: ty.Optional[ConstantFn] = None,
                 value: ty.Optional[ty.List[ty.Any]] = None,
                 replicate_size: int = 1,
                 call_order_idx: ty.Optional[int] = None,
                 deterministic: bool = False,
                 clamped: bool = False,
                 parent_nodes: ty.Optional[ty.List[ty.Any]] = None):

        self.sampling_dn = sampled_from  # dn object
        self.constant_fn = returned_from  # constant fn object
        self.is_sampled = False
        # not value checks for both [] and None
        generated_or_specified_value: ty.List[ty.Any] = \
            self.generate_value() if not value else value
        
        super().__init__(node_name,
                         sample_size=sample_size,
                         value=generated_or_specified_value,
                         replicate_size=replicate_size,
                         call_order_idx=call_order_idx,
                         deterministic=deterministic,
                         clamped=clamped,
                         parent_nodes=parent_nodes)

        # used for MCMC #move/operator setup during inference
        self.operator_weight = 0.0

        self.populate_operator_weight()

        # pgm specs
        if self.sampling_dn:
            self.is_sampled = True  # r.v. value is sampled

    def generate_value(self) -> ty.List[ty.Any]:
        """Generate value."""

        if self.sampling_dn:
            return self.sampling_dn.generate()
        
        elif self.constant_fn:
            return self.constant_fn.generate()

        else:
             raise RuntimeError(("Cannot generate value. No distribution nor "
                                 "constant function provided."))

    def __str__(self) -> str:
        return super().__str__()

    def __lt__(self, other):
        return super().__lt__(other)

    def plot_node(self,
                  axes: plt.Axes,
                  sample_idx: int = 0,
                  repl_idx: int = 0,
                  repl_size: int = 1,
                  branch_attr: str = "state",
                  draw_reconstructed: ty.Optional[bool] = False) -> None:
        """Plot node (side-effect) on provided Axes object

        Args:
            axes (matplotlib.pyplot.Axes): Axes object where we are
                drawing the tree.
            sample_idx (int): Which sample to plot. Defaults to 0.
            repl_idx (int): Which replicate to plot (one
                replicated is plotted at a time). Defaults to 0.
            repl_size (int): How many scalar random variables
                to plot at a time. Defaults to 1.
            branch_attr (str, optional): If tree, which branch attribute
                to color branch according to. Defaults to 'state'.
            draw_reconstructed (bool, optional): Flag specifying if to
                draw reconstructed tree instead of complete tree.
                Defaults to 'False'.
        """

        # if list
        if super().__len__() >= 1 and isinstance(self.value, list):
            # if tree, we only plot one
            if isinstance(self.value[0], pjtr.AnnotatedTree):
                # so mypy won't complain
                if isinstance(sample_idx, int) and \
                    isinstance(self.value[sample_idx * repl_size + repl_idx],
                               pjtr.AnnotatedTree):
                    if (self.value[0].continuous_trait is not None or
                            (self.value[0].state_count is not None and self.value[0].state_count > 1)):
                        self.value[sample_idx * repl_size + repl_idx]\
                            .plot_node(axes,
                                       node_attr=branch_attr,
                                       draw_reconstructed=draw_reconstructed)

                    else:
                        self.value[sample_idx * repl_size + repl_idx] \
                            .plot_node(axes,
                                       draw_reconstructed=draw_reconstructed)

            # not a tree
            else:
                # so mypy won't complain
                hist_vals = [ty.cast(float, v) for v in self.value]

                # one sample
                if sample_idx is not None and self.sampling_dn:
                    plot_node_histogram(axes,
                                        hist_vals,
                                        sample_idx=sample_idx,
                                        repl_size=repl_size)

                # all samples
                else:
                    plot_node_histogram(axes,
                                        hist_vals,
                                        repl_size=repl_size)

        # if scalar
        elif isinstance(self.value, (int, float, str)):
            # return
            plot_node_histogram(axes, [float(self.value)])

    def populate_operator_weight(self) -> None:
        if isinstance(self.value, (list, np.ndarray)):
            if isinstance(self.value[0], (int, float, str, np.float64)):
                self.operator_weight = 1  # this is a scalar random variable

            # has objects like DiscreteStateDependentRate inside list of values
            else:
                # TODO: later see how rev moves 2D-arrays and tree nodes
                # print("value has objects inside")
                # print(self.value)
                pass

        else:
            # TODO: later see how rev moves 2D-arrays and tree nodes
            raise RuntimeError(
                ("Could not determine dimension of StochasticNodeDAG when"
                 " figuring out operator weight. Exiting..."))

    def get_node_stats_str(self, start: int, end: int, repl_idx: int) -> str:
        return super().get_node_stats_str(start, end, repl_idx)


class DeterministicNodeDAG(NodeDAG):
    """Derived DAG node that holds value dependent on another node's.
    
    Deterministic node values depend deterministically on those from
    parent nodes. These nodes are used to modify and annotate
    stochastic node values.
    """

    def __init__(self,
                 node_name: str,
                 value: ty.Optional[ty.List[ty.Any]] = None,
                 call_order_idx: ty.Optional[int] = None,
                 deterministic: bool = True,
                 parent_nodes: ty.Optional[ty.List[NodeDAG]] = None) -> None:

        sample_size = None
        if isinstance(value, list):
            sample_size = len(value)

        super().__init__(node_name,
                         sample_size=sample_size,
                         replicate_size=1,
                         value=value,
                         call_order_idx=call_order_idx,
                         deterministic=deterministic,
                         parent_nodes=parent_nodes)

        self.is_sampled = False

    def __str__(self) -> str:
        return super().__str__()

    def __lt__(self, other) -> bool:
        return super().__lt__(other)

    def __len__(self) -> int:
        if isinstance(self.value, list):
            return super().__len__()

        else:
            return 0

    def plot_node(self,
                  axes: plt.Axes,
                  sample_idx: ty.Optional[int] = 0,
                  repl_idx: ty.Optional[int] = 0,
                  repl_size: ty.Optional[int] = 1,
                  branch_attr: ty.Optional[str] = "state",
                  draw_reconstructed: ty.Optional[bool] = False) -> None:

        # if list
        if super().__len__() >= 1 and isinstance(self.value, list):
            # if tree, we only plot one
            if isinstance(self.value[0], pjtr.AnnotatedTree):
                # so mypy won't complain
                if isinstance(sample_idx, int) and \
                        isinstance(self.value[sample_idx * repl_size + repl_idx],
                                   pjtr.AnnotatedTree):
                    if (self.value[0].continuous_trait is not None or
                            (self.value[0].state_count is not None and self.value[0].state_count > 1)):
                        self.value[sample_idx * repl_size + repl_idx] \
                            .plot_node(axes, node_attr=branch_attr)

                    else:
                        self.value[sample_idx * repl_size + repl_idx] \
                            .plot_node(axes)

        else:
            plot_blank(axes)

    def populate_operator_weight(self):
        pass

    def get_node_stats_str(self, start: int, end: int, repl_idx: int) -> str:
        return ""


############
# Plotting #
############
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
    axes.cla()
    axes.patch.set_alpha(0.0)
    axes.xaxis.set_ticks([])
    axes.yaxis.set_ticks([])
    axes.spines['left'].set_visible(False)
    axes.spines['bottom'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)


####################
# Helper functions #
####################

def extract_value_from_dagnodes(dag_node_list: ty.List[NodeDAG]) \
        -> ty.List[float]:
    """_summary_

    Args:
        dag_node_list (NodeDAG): List of NodeDAG objects.

    Raises:
        ec.NoPlatingAllowedError: _description_
        ec.StateDependentParameterMisspec: _description_

    Returns:
        (list): List of values extracted from the DAG node.
    """

    many_nodes_dag = len(dag_node_list) > 1
    v_list: ty.List[float] = []

    for node_dag in dag_node_list:

        # so mypy won't complain
        if isinstance(node_dag, NodeDAG):
            # no plating supported
            if node_dag.repl_size > 1:
                raise ec.NoPlatingAllowedError(
                    "sse_rate", node_dag.node_name)

            v = node_dag.value  # list (I think before I also
                                # allowed numpy.ndarray, but not anymore)

            # so mypy won't complain
            if isinstance(v, list):
                if len(v) > 1 and many_nodes_dag:
                    raise ec.StateDependentParameterMisspec(
                        message=(
                            ("If many variables are passed as arguments "
                             "to initialize another variable, each of these "
                             "variables can contain only a single value")))

                elif len(v) == 1 and many_nodes_dag:
                    # making list longer
                    # (v should be a list, which is why I don't use append)
                    v_list += v

                elif len(v) >= 1 and not many_nodes_dag:
                    return v

    return v_list


def extract_vals_as_str_from_node_dag(
        val_list: ty.List[ty.Union[str, NodeDAG]]) -> ty.List[str]:
    """Get values from DAG node.

    If all elements are strings, returns copy of 'val_list'.
    When elements are StochasticNodeDAGs, replaces those objects
    by their values after casting to string (their values must be
    within a list).

    Raise:
        VariableMisspec: Is raised if DAG node does not have a value
            or if it cannot be string-fied.

    Returns:
        (str): List of values as strings.
    """

    extracted_val_list: ty.List[str] = list()

    for v in val_list:
        if isinstance(v, str):
            extracted_val_list.append(v)

        elif isinstance(v, StochasticNodeDAG) and v.value:
            try:
                extracted_val_list.extend([str(i) for i in v.value])

            except (AttributeError, TypeError) as e:
                raise ec.VariableMisspec(str(v))

    # will be empty if DeterministicNodeDAG is in val_list
    return extracted_val_list


if __name__ == "__main__":
    pass
