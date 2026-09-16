import typing as ty
import numpy as np  # type: ignore
import dendropy as dp  # type: ignore
import random
import time

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.utility.helper_functions as pjh
import phylojunction.utility.exception_classes as ec
import phylojunction.distribution.dn_parametric as dnpar
from phylojunction.data.tree import AnnotatedTree
from phylojunction.data.sampled_ancestor \
    import SampledAncestor  # type: ignore
from phylojunction.data.attribute_transition \
    import AttributeTransition  # type: ignore

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class DnSSE(pgm.DistrForSampling):
    """Discrete SSE distribution.

    Class for the discrete state-dependent speciation and extinction
    (SSE) distribution. Used for sampling phylogenetic trees annotated
    with discrete states.
    
    The simulate() method in this class is a 'rising tide' sampler,
    where all living lineages grow together. It is not a recursive
    sampler in which lineages take care of growing themselves and recur
    upon birth events.

    At all times we know:
        (i)  how many lineages are alive and dead
        (ii) all the character states represented by living lineages

    Parameters:
        n_sim (int): Number of trees to sample. Defaults to 1.
        n_repl (int): Number of tree replicates per sample. This
            is equivalent to the size of a plate in the DAG
            representation. Defaults to 1.
        with_origin (bool): Flag for whether the process starts at the
            origin or not.
        root_is_born (bool): Attribute that records if root node was
            created and is in tree (for when it starts at the origin).
            'True' if process starts at the root.
        start_states (int): List of integers representing the starting
            states of each of the 'n_sim' samples.
        seed_age (float, optional): Age of seed node (either origin or
            root).
        condition_on_speciation (bool): Flag for rejecting tree samples
            that do not go a single speciation event before 'stop_val'
            is met. Note that this first speciation event may or not be
            what is canonically the reconstructed tree root node. If
            'True', rejects tree sample failing to meet condition.
            Defaults to 'False'.
        condition_on_survival (bool): Flag for rejecting tree samples
            that go extinct before 'stop_val' is met. If 'True',
            rejects tree sample failing to meet condition. Defaults to
            'True'.
        condition_on_obs_both_sides_root (bool): Flag for rejecting
            tree samples that do not have observed nodes on both sides
            of the complete tree's root node. If 'True', rejects tree
            sample failing to meet condition. Defaults to 'False'.
        stop (str): Stop condition to end sampling (simulation)
            procedure and return tree. If 'age', stops when age of
            either origin or is equal to 'stop_val' (see below).
            If "size", stops when tree has 'stop_val' observed nodes.
            Size stopping requires complete sampling (all sampling probabilities equal 1).
        stop_val (float): List of values used by 'stop' (see above) to
            end each of the 'n_sim' sampling (simulation) procedures
            and return tree. Either maximum age, or maximum count of
            observable nodes.
        min_rec_taxa (int): Required minimum number of observed taxa in
            reconstructed tree. Defaults to 0.
        max_rec_taxa (int): Required maximum number of observed taxa in
            reconstructed tree. Defaults to 1e12.
        abort_at_alive_count (int): Number of living (not observed!)
            nodes at which point sample is rejected. This parameter is
            used to abort samples whose SSE parameters cause trees to
            grow too large. Defaults to 1e12.
        sse_stash (SSEStash): Object holding all discrete
            state-dependent parameters we need to sample (i.e.,
            simulate). This object holds (i) the number of discrete
            states of the process, (ii) the number of time slices if
            process is time-heterogeneous, (iii) the end age of each
            time slice.
        events (MacroevolEventHandler): Object that computes total
            rate values, and samples SSE events. It is a member of
            sse_stash.
        state_count (int): Number of states of SSE process.
        n_time_slices (int): Number of time slices (epochs).
        slice_t_ends (float, optional): List of floats with time slice
            time ends (forward!).
        prob_handler (DiscreteStateDependentProbabilityHandler): Object
            that takes care of state-dependent taxon sampling across
            time slices.
        epsilon (float, optional): Float threshold to determine if a
            tiny decimal number is to be considered 0.0 or not. In
            other words, if the difference between a tiny value 'x'
            and 0.0 is smaller than epsilon, then 'x' is set to 0.0.
            Defaults to 1e-12.
        runtime_limit (int, optional): Runtime ceiling (in seconds)
            for obtaining the 'n' tree samples. If this limit is met,
            the sampling procedure is aborted. Defaults to 300.
        max_n_failed_attempts (int, optional): Maximum number of failed
            tree sampling attempts (replicates included) before
            PhyloJunction quits. Defaults to 200.
        rng_seed (int, optional): Integer seed for the two random
            number generators used in this class. This seed is only
            ever used by user if bypassing the scripting language
            (otherwise random number generator seeds are handled in
            cmd_parse.py). Defaults to None.
        debug (bool, optional): Flag for whether to print debugging
            messages during sampling procedure.
        info (bool, optional): Flag for whether to print information
            about running sampling procedure.
    """

    DN_NAME = "DnSSE"

    # for typing #
    # basic simulation parameters
    n_sim: int
    n_repl: int

    # start conditions
    with_origin: bool
    root_is_born: bool
    start_states: ty.List[int]
    seed_age: ty.Optional[float]

    # conditioning flags and values
    condition_on_speciation: bool
    condition_on_survival: bool
    condition_on_obs_both_sides_root: bool
    
    # stop condition and values
    stop: str
    stop_val: ty.List[float]
    min_rec_taxa: int
    max_rec_taxa: int
    abort_at_alive_count: int
    
    # model parameters
    sse_stash: sseobj.SSEStash
    events: sseobj.MacroevolEventHandler
    state_count: int
    n_time_slices: int
    slice_t_ends: ty.List[ty.Optional[float]]
    prob_handler: sseobj.DiscreteStateDependentProbabilityHandler
    
    # other specs
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
                 max_n_failed_attempts: int = 200,
                 rng_seed: ty.Optional[int] = None,
                 debug: ty.Optional[bool] = False,
                 info: ty.Optional[bool] = False) -> None:

        # basic simulation parameters
        self.n_sim = int(n)  # number of full model samples
        self.n_repl = int(n_replicates)  # plating (iid r.v.'s)!
        
        # start conditions
        self.with_origin = origin
        self.root_is_born = False
        self.start_states = start_states_list
        
        # stop condition, stop value
        self.stop = stop
        self.stop_val = stop_value
        
        # conditioning flags
        self.condition_on_speciation = condition_on_speciation
        self.condition_on_survival = condition_on_survival
        self.condition_on_obs_both_sides_root = \
            condition_on_obs_both_sides_root

        # rejection sampling
        self.min_rec_taxa = min_rec_taxa
        self.max_rec_taxa = max_rec_taxa
        self.abort_at_alive_count = abort_at_alive_count

        # model parameters #
        # macroevol rate handler and prob handler inside
        self.sse_stash = sse_stash

        # all rates, number of states and of slices
        self.events = sse_stash.get_meh()
        self.state_count = self.events.state_count
        self.n_time_slices = self.events.n_time_slices
        self.slice_t_ends = self.events.slice_t_ends  # for time-het processes
        self.seed_age = self.events.seed_age  # required for checking input
        # sampling probabilities, if provided
        self.prob_handler = sse_stash.get_prob_handler()

        # other specs #
        self.epsilon = epsilon
        self.runtime_limit = runtime_limit
        self.max_n_failed_attempts = max_n_failed_attempts
        self.debug = debug
        self.info = info
        # handling random number generator seeds
        # (this is only ever used if this distribution
        # is called without using cmd_parse.py)
        self.rng_seed = rng_seed
        if self.rng_seed is not None:
            np.random.seed(seed=self.rng_seed)  # for event times
            random.seed(self.rng_seed)  # for choosing lineages

        # checking number of provided values (vectorizing if necessary)
        self.init_check_vectorize_sample_size()

        # checking input is OK
        self._check_input_health()

    ##########################
    # Initialization methods #
    ##########################

    # from parent class
    def init_check_vectorize_sample_size(self) -> None:
        """Vectorize SSE rates and probs, start values and stop values.

        This method is in the parent class.

        Raises:
            DimensionalityError: Is raised if (i) more than one SSE
                rate value is provided, but that number is smaller
                than the requested number of tree samples (because
                then we do not know how to vectorize the values), (ii)
                same as (i), but for SSE probabilities, (iii) same as
                (i), but for starting states, (iv) same as (i), but
                for stop values
        """

        # changing atomic_rate_params_matrix here
        # should affect state_dep_rate_manager
        #
        # 1D: time slices, 2D: list of atomic rate params
        sse_rate_params_mat = \
            self.events.sse_rate_manager.matrix_state_dep_params

        sse_prob_params_mat = \
            self.prob_handler.state_dep_prob_manager.matrix_state_dep_params

        ############################
        # Checking rate parameters #
        ############################
        # iterating over time slices
        for list_sse_rate_params in sse_rate_params_mat:
            for sse_rate_param in list_sse_rate_params:
                n_val = len(sse_rate_param.value)

                # if one sse rate value is provided, but we want > 1 tree
                # samples, we multiply the single provided value
                if n_val == 1 and (self.n_sim > 1):
                    sse_rate_param.value = \
                        [sse_rate_param.value[0] for i in range(self.n_sim)]

                # do not know how to multiply, error!
                elif n_val > 1 and n_val != self.n_sim:
                    raise ec.DimensionalityError(
                        self.DN_NAME,
                        par_name="'sse_rate' -- its values")

        ###################################
        # Checking sampling probabilities #
        ###################################
        # iterating over time slices
        for list_state_dep_probs in sse_prob_params_mat:
            for state_dep_prob in list_state_dep_probs:
                n_val = len(state_dep_prob.value)

                # if one sse prob value is provided, but we want > 1 tree
                # samples, we multiply the single provided value
                if n_val == 1 and self.n_sim > 1:
                    state_dep_prob.value = \
                        [state_dep_prob.value[0] for i in range(self.n_sim)]

                # do not know how to multiply, error!
                elif n_val > 1 and n_val != self.n_sim:
                    raise ec.DimensionalityError(
                        self.DN_NAME,
                        par_name="'sse_prob' -- its values")

        ############################
        # Checking starting states #
        ############################
        # Replicates are consecutive within each sample, matching generate()'s indexing.
        # Build a new list so expansion does not modify the caller's starting states.
        n_start_states = len(self.start_states)
        if n_start_states == 1:
            self.start_states = [int(self.start_states[0])] * (self.n_sim * self.n_repl)
        elif n_start_states == self.n_sim:
            self.start_states = [state for state in self.start_states for _ in range(self.n_repl)]
        elif n_start_states == self.n_sim * self.n_repl:
            self.start_states = list(self.start_states)
        elif n_start_states > 0:
            raise ec.DimensionalityError(self.DN_NAME, par_name="'start_state'")

        ########################
        # Checking stop values #
        ########################
        n_stop_vals: int = len(self.stop_val)
        # multiplying stop values if only one was passed
        # and the number of simulations is > 1
        if n_stop_vals == 1 and self.n_sim > 1:
            self.stop_val = \
                [float(self.stop_val[0]) for i in range(self.n_sim)]

        # do not know how to multiply, error!
        elif n_stop_vals > 1 and n_stop_vals != self.n_sim:
            raise ec.DimensionalityError(
                self.DN_NAME,
                par_name="'stop_value'")

        # dummy return for mypy
        return None

    def _check_input_health(self) -> None:
        """Check the validity of many DnSSE's initialization arguments.

        Raises:
            ObjInitInvalidArgError: Is raised if (i) the requested
                number of samples ('n_sim') is < 1, (ii) no
                starting states are provided, (iii) the stop condition
                was neither 'age' nor 'size', (iv) more than one time
                slice was specified, but then 'age' was not the stop
                condition (in which case, we cannot anchor the absolute
                ages), (v) the seed age provided to the distribution
                does not match the seed age provided to the SSE stash,
                (vi) number of taxa or maximum age at which to stop
                growing the tree is not an integer (or float,
                respectively), or it is < 0 (< 0.0, respectively).
            ObjInitMissingParameterError: Is raised if more than one
                time slice is specified, and no seed age is provided
                to SSEStash via 'event' (in which case we cannot anchor
                the tree onto an absolute timescale to compare it to
                the time slice age ends).
        """
        if self.n_sim < 1:
            raise ec.ObjInitInvalidArgError(
                self.DN_NAME,
                "'n_sim'",
                "Please specify a number of simulations >= 1.")

        if not self.start_states:
            raise ec.ObjInitInvalidArgError(
                self.DN_NAME,
                "'start_states_list'",
                "You must provide a list of " + str(self.n_sim)
                + " starting states.")

        if self.stop != "age" and self.stop != "size":
            raise ec.ObjInitInvalidArgError(
                self.DN_NAME,
                "'stop'",
                "Stop condition must be \"age\" or \"size\".")

        if self.n_time_slices > 1 and not self.stop == "age":
            raise ec.ObjInitInvalidArgError(
                self.DN_NAME,
                "'age'",
                ("If you specified more than a single time slice, you "
                 "need to provide a maximum age as stopping condition "
                 "to anchor the absolute slice ages."))

        if self.n_time_slices > 1 and not self.seed_age:
            raise ec.ObjInitMissingParameterError(
                self.DN_NAME,
                "'events'",
                ("When providing time-slice age ends, you must give a "
                 "seed age (origin or root) to the SSEStash to anchor "
                 "the tree and ensure the tree spans those slices."))

        if self.stop == "age" and self.seed_age:
            # TODO: make sure seed age in FIGManager is vectorized and
            # that in script, the same variable is used for both
            # DiscreteStateDependentParameterManager and stop condition
            # value in the spec of this distribution
            for idx, a_stop_val in enumerate(self.stop_val):
                # TODO: will need to do self.seed_age[idx]
                # later when vectorization is added
                if a_stop_val != self.seed_age:
                    raise ec.ObjInitInvalidArgError(
                        self.DN_NAME,
                        "'stop'",
                        ("The max. stopping age(s) for the simulation "
                         "must match the specified seed (origin or root) "
                         "age(s)."))

        # the checks below are also carried out at the grammar level,
        # but we do it again in case distribution is used without a script
        if self.stop == "size":
            # Size stopping currently has no incomplete-sampling observation step.
            probabilities = self.prob_handler.state_dep_prob_manager.matrix_state_dep_params
            if any(value != 1.0 for epoch in probabilities for prob in epoch for value in prob.value):
                raise ec.ObjInitInvalidArgError(
                    self.DN_NAME, "'stop'",
                    'Size stopping requires complete sampling (all sampling probabilities must be 1).')

            for idx, a_stop_val in enumerate(self.stop_val):
                # must be a number
                if not isinstance(a_stop_val, (int, float)):
                    raise ec.ObjInitInvalidArgError(
                        self.DN_NAME,
                        "'stop'",
                        ("Stop condition value (number of terminal nodes) "
                         "must be a number."))

                # it is a number...
                else:
                    # can be a int-convertible float
                    if isinstance(a_stop_val, float):
                        if a_stop_val.is_integer():
                            self.stop_val[idx] = int(a_stop_val)

                        # float must be convertible to integer
                        else:
                            raise ec.ObjInitInvalidArgError(
                                self.DN_NAME,
                                "'stop'",
                                ("Stop condition value (number of terminal "
                                 "nodes) could not be converted to integer."))

                    # ideally should be int
                    elif isinstance(a_stop_val, int):
                        self.stop_val[idx] = a_stop_val

                    # once self.stop_val[idx] is initialized, it must be >= 0
                    if self.stop_val[idx] < 0:
                        raise ec.ObjInitInvalidArgError(
                            self.DN_NAME,
                            "'stop'",
                            ("Stop condition value (number of terminal nodes) "
                             "cannot be negative."))

                    # if it is 1, then the process cannot start at a root node,
                    # which by definition has two children already
                    if self.stop_val[idx] == 1 and not self.with_origin:
                        # TODO: add this check to the grammar level
                        raise ec.ObjInitInvalidArgError(
                            self.DN_NAME,
                            "'stop_value'",
                            ("Stop condition value (number of terminal "
                             "nodes) cannot be 1 for a process starting "
                             "at the root. This is only allowed if the "
                             "process starts at the origin."))

        if self.stop == "age":
            for idx, a_stop_val in enumerate(self.stop_val):
                # must be an integer
                if not isinstance(a_stop_val, float):
                    raise ec.ObjInitInvalidArgError(
                        self.DN_NAME,
                        "'stop'",
                        ("Stop condition value (tree height) must be "
                         "an integer"))

                else:
                    self.stop_val[idx] = a_stop_val
                    # must be >= 0
                    if self.stop_val[idx] < 0.0:
                        raise ec.ObjInitInvalidArgError(
                            self.DN_NAME,
                            "\'stop\'",
                            ("Stop condition value (tree height) cannot be "
                             "negative."))

    ############################
    # Large simulation methods #
    ############################

    # NOTE: Python passes arguments 'by assignment', which means that for
    # mutable arguments, we are doing something akin to 'pass by reference';
    # we can mutate -- but not reassign! -- the argument variable and
    # changes will be reflected outside the function.
    #
    # In the following executers, we make use of the above Python feature.
    # Class members are passed as arguments and modified (side-effect) inside
    # the executers. We choose to do so to make the behavior of each of the
    # executers more explicit by just looking at the method's signature

    def _germinate_tree(self,
                        a_start_state: int,
                        with_origin: bool,
                        state_representation_dict: ty.Dict[int, ty.Set[str]],
                        untargetable_node_set: ty.Set[str]) \
                            -> ty.Tuple[dp.Tree, dp.Node, int, bool]:
        """
        Initialize (germinate) tree with either origin or root.

        Args:
            a_start_state (int): State at seed node (origin or root).
            with_origin (bool): Flag specifying whether tree is to be
                germinated at origin node (we automatically attach the
                'brosc' node), or at the root node (no 'brosc' node).
            state_representation_dict (dict): Dictionary that keeps
                track of all states currently represented by living
                lineages. Keys are integer representing the states,
                values are sets of taxon names (strings).
            untargetable_node_set (set): Set of Node labels that
                cannot be targeted for events anymore (went extinct).

        Returns:
            (tuple): Tuple with four items. The germinated tree
                (dendropy.Tree), the 'brosc' node if there is one
                (dendropy.Node), the cumulative node count (integer),
                and a boolean specifying if a root exists or not.
        """

        cumulative_node_count: int = 0
        brosc_node: dp.Node = None
        root_is_born: bool = False

        # simulation starts at origin
        if with_origin:
            ##################################################################
            # brosc: before-root-origin-single-child                         #
            # Node that might be necessary for tying process off if a        #
            # speciation never happens before the process either dies, or    #
            # reaches stop condition. It receives the same index as the      #
            # origin because, unless an ancestor is sampled, the brosc node  #
            # is essentially the origin -- it is a trick to tie off the tree #
            #                                                                #
            # Note: brosc can become an SA lineage node if ancestor sampling #
            # happens, but a root is never born                              #
            ##################################################################

            # this particular brosc_node instance below will be used
            # within the simulate() method, but not inside
            # "execute_sample_ancestor()", which has its own instance of
            # brosc_node being created
            brosc_node = dp.Node(taxon=dp.Taxon(label="brosc"),
                                label="brosc",
                                edge_length=0.0)
            brosc_node.state = a_start_state
            brosc_node.annotations.add_bound_attribute("state")
            brosc_node.index = cumulative_node_count  # brosc_node is
            brosc_node.annotations.add_bound_attribute("index")
            brosc_node.sampled = True
            brosc_node.is_sa = False
            brosc_node.is_sa_dummy_parent = False
            brosc_node.is_sa_lineage = False

            # origin node will be the seed_node
            origin_node = dp.Node(taxon=dp.Taxon(label="origin"),
                                label="origin",
                                edge_length=0.0)
            origin_node.state = a_start_state
            origin_node.annotations.add_bound_attribute("state")
            origin_node.index = cumulative_node_count
            origin_node.annotations.add_bound_attribute("index")
            origin_node.alive = True
            origin_node.sampled = False
            origin_node.is_sa = False
            origin_node.is_sa_dummy_parent = False
            origin_node.is_sa_lineage = False
            state_representation_dict[origin_node.state] \
                .add(origin_node.label)

            # now make tree #
            # will remain a "is_leaf() == True" until we add children
            tr = dp.Tree(seed_node=origin_node)
            tr.taxon_namespace.add_taxon(origin_node.taxon)

        # simulation starts at root
        else:
            # root node will be the seed_node,
            # and will be untargetable from the get-go
            root_node = dp.Node(taxon=dp.Taxon(label="root"),
                                label="root",
                                edge_length=0.0)
            root_node.state = a_start_state
            root_node.annotations.add_bound_attribute("state")
            root_node.index = cumulative_node_count
            root_node.annotations.add_bound_attribute("index")
            root_node.alive = False
            root_node.sampled = False
            root_node.is_sa = False
            root_node.is_sa_dummy_parent = False
            root_node.is_sa_lineage = False
            untargetable_node_set.add("root")
            root_is_born = True

            # we will actually start off with the
            # left and right nodes being targetable
            cumulative_node_count += 1
            left_node = dp.Node(taxon=dp.Taxon(label="nd1"),
                                label="nd1",
                                edge_length=0.0)
            left_node.is_sa = False
            left_node.is_sa_dummy_parent = False
            left_node.is_sa_lineage = False
            left_node.alive = True
            left_node.sampled = True
            left_node.state = root_node.state  # get state from parent (root)
            left_node.annotations.add_bound_attribute("state")
            left_node.index = cumulative_node_count
            left_node.annotations.add_bound_attribute("index")
            root_node.add_child(left_node)
            state_representation_dict[left_node.state].add(left_node.label)

            cumulative_node_count += 1
            right_node = dp.Node(taxon=dp.Taxon(label="nd2"),
                                label="nd2",
                                edge_length=0.0)
            right_node.is_sa = False
            right_node.is_sa_dummy_parent = False
            right_node.is_sa_lineage = False
            right_node.alive = True
            right_node.sampled = True
            right_node.state = root_node.state  # get state from parent (root)
            right_node.annotations.add_bound_attribute("state")
            right_node.index = cumulative_node_count
            right_node.annotations.add_bound_attribute("index")
            root_node.add_child(right_node)
            state_representation_dict[right_node.state].add(right_node.label)

            # now make tree #
            # will remain a "is_leaf() == True" until we add children
            tr = dp.Tree(seed_node=root_node)
            tr.taxon_namespace.add_taxon(root_node.taxon)
            tr.taxon_namespace.add_taxon(left_node.taxon)
            tr.taxon_namespace.add_taxon(right_node.taxon)
    
        return tr, brosc_node, cumulative_node_count, root_is_born
    
    def _update_sa_lineage_dict(
            self,
            a_time: float,
            sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]],
            sa_lineage_node_labels: ty.List[str],
            debug: bool = False) -> None:
        """
        Update sa_lineage_dict class member as a side-effect.

        This function is called whenever a branch stops growing:
            (i) When a node undergoes a birth event (effectively
            becoming an internal node);
            (ii) When living nodes reach the present moment or a stop
            condition (i.e., the tree stops growing).

        These two possibilities only matter when the node(s) involved
        have direct (sampled) ancestors on the branch leading to them
        (i.e., if for that node .is_sa_lineage == True).

        In case (i), 'sa_lineage_node_label' will have a single
        element, the label of the node undergoing birth. In case (ii),
        'sa_lineage_node_label' will have the labels of all terminal
        living nodes for which .is_sa_lineage == True.

        By updating sa_lineage_dict as the tree is built, we don't have
        to do tree traversals later to get sampled ancestor times for
        plotting (when initializing an AnnotatedTree).

        Args:
            a_time (float): Either the time of the last birth event
                undergone by an SA lineage or the simulation end time.
            sa_lineage_dict (dict): Dictionary that keeps track of
                nodes that have direct (sampled) ancestor children.
                Keys are taxon names (strings) and values are lists
                of SampledAncestor objects.
            sa_lineage_node_labels (str): List of node labels whose
                subtending branch has direct (sampled) ancestors.
        """

        if debug:
            print("\n>> Updating SA lineage dictionary. Keys are is \'" +
                  + ", ".join(sa_label
                              for sa_label in sa_lineage_node_labels) +
                  + "\' and time being added is " + str(a_time) +
                  + "\n   SA lineage dict:")
            print("      " + "\n      ".join(
                k + " has in its lineage " +
                + ", ".join(sa.label for sa in sas) for k, sas
                in sa_lineage_dict.items()))

        for sa_lineage_node_label in sa_lineage_node_labels:
            sa_list = sa_lineage_dict[sa_lineage_node_label]

            for sa in sa_list:
                sa.time_to_lineage_node = a_time - sa.global_time

    def _execute_birth(
            self,
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
            debug=False) -> ty.Tuple[dp.Node, int]:
        """Execute lineage birth.
         
        This method has both side-effects and a return.
        The side-effects include the updating of class members. These
        members are passed as arguments in the signature so that we
        can have a better idea of what to expect as behavior:
        
            (i)   self.tr.tr_namespace,
            (ii)  self.state_representation_dict
            (iii) self.sa_lineage_dict
            (iv)  self.state_transition_dict
            (v)   self.clado_state_transition_dict
            (vi)  self.root_is_born
            (vii) self.untargetable_node_set

        Args:
            tr_namespace (dendropy.TaxonNamespace): Dendropy object
                recording taxa in the tree.
            chosen_node (dendropy.Node): Node that will undergo
                birth event.
            state_representation_dict (dict): Dictionary that keeps
                track of all states currently represented by living
                lineages. Keys are integer representing the states,
                values are sets of taxon names (strings).
            sa_lineage_dict (dict): Dictionary that keeps track of
                nodes that have direct (sampled) ancestor children.
                Keys are taxon names (strings), values are lists
                of SampledAncestor objects.
            state_transition_dict (dict): Dictionary that keeps track
                of nodes subtending which state transitions happen
                (used for plotting). Keys are taxon names (strings)
                and values are lists of AttributeTransition objects.
            clado_state_transition_dict (dict): Dictionary that keeps
                track of speciating nodes whose children have different
                states (used for stochastic maps I/O). Keys are
                taxon names (strings) and values are
                AttributeTransition objects.
            untargetable_node_set (str): Set of Node labels that
                cannot be targeted for events anymore (went extinct).
            cumulative_node_count (int): Total number of non-direct
                ancestor nodes in the tree (to be used in labeling).
            node_index (int): Total node count to be used to number
                nodes (used in their 'index' attribute).
            sse_birth_rate_object (DiscreteStateDependentRate): SSE
                rate parameter object holding the departing/arriving
                states.
            event_t (float): Time of birth event taking place.
            debug (bool): Flag for printing debugging messages.
                If 'True', prints messages. Defaults to 'False'.

        Returns:
            (tuple): Tuple with two objects, the last node to speciate
            (dendropy.Node) and the tree's cumulative node count (int).
        """

        left_arriving_state, right_arriving_state = \
            sse_birth_rate_object.arriving_states

        if debug:
            print("> SPECIATION of node " + chosen_node.label
                  + " in state " + str(chosen_node.state)
                  + " into daughters with states "
                  + str(left_arriving_state) + " and "
                  + str(right_arriving_state))

        ############################################################
        # Special case: first speciation event (root must created) #
        ############################################################
        if chosen_node.label in ("origin", "brosc"):
            node_index += 1
            root_node = dp.Node(taxon=dp.Taxon(label="root"),
                                label="root")
            # root takes origin or brosc state
            root_node.state = chosen_node.state
            root_node.annotations.add_bound_attribute("state")

            # adding index (origin and brosc will have index 0, so root should be 1)
            root_node.index = node_index
            root_node.annotations.add_bound_attribute("index")

            # we will pick the root as the first event
            root_node.alive = False
            root_node.sampled = False
            root_node.is_sa = False
            root_node.is_sa_dummy_parent = False
            root_node.is_sa_lineage = chosen_node.is_sa_lineage
            tr_namespace.add_taxon(root_node.taxon)
            self.root_is_born = True

            # origin/brosc cannot be selected anymore
            chosen_node.alive = False  # origin/brosc node is no longer alive
            # origin/brosc is an internal node, cannot be sampled
            chosen_node.sampled = False
            untargetable_node_set.add(chosen_node.label)  # cannot be targeted

            state_representation_dict[chosen_node.state] \
                .remove(chosen_node.label)  # and does not represent state

            # updating origin
            if chosen_node.label == "origin":
                # at this point, the origin was chosen to undergo a birth
                # event, but this node will never be extended (the origin
                # always has an origin_edge_length = 0.0)
                #
                # the root does NOT exist yet here; for us to account for
                # the evolution (branch length) that has happened between
                # the origin/brosc to this moment (when the root is born),
                # we must add this branch length (event_t) now when creating
                # the root (unlike other nodes, which are created with
                # edge_length = 0.0)
                root_node.edge_length = event_t
                # finally link root to origin
                chosen_node.add_child(root_node)

            # root replaces brosc (dummy_node might or not exist)
            # [ori] ---> (dummy_node) ---> [brosc]
            # or
            # [ori] ---> [brosc]
            # becomes
            # [ori] ---> (dummy_node) ---> [root]
            # and
            # [ori] ---> [root], respectively
            elif chosen_node.label == "brosc":
                # must add the evolution leading up to the brosc_node to the
                # root node edge length
                #
                # (note that the brosc_node edge length will always be 0.0 if
                # it resulted from an ancestor sampling event, but it could
                # be > 0.0 if a state transition event happened)

                root_node.edge_length = chosen_node.edge_length

                # because the extend_all_lineages() method does
                # extend brosc_node
                root_node.is_sa_lineage = chosen_node.is_sa_lineage

                # if true, means brosc's parent is a dummy node,
                # and brosc's sister taxon is a SA
                if root_node.is_sa_lineage:
                    sa_lineage_dict["root"] = sa_lineage_dict.pop("brosc")

                # could be origin, or dummy node
                # (in which case, brosc would have a SA sister taxon)
                brosc_parent_node = chosen_node.parent_node

                # brosc is obliterated
                brosc_parent_node.remove_child(chosen_node)

                # we add the root node in its place
                brosc_parent_node.add_child(root_node)

                if "brosc" in state_transition_dict:
                    state_transition_dict["root"] = \
                        state_transition_dict["brosc"]
                    del state_transition_dict["brosc"]

            # now make chosen_node become the root,
            # so the rest of the birth event can be executed
            chosen_node = root_node

        # normal birth event code (assuming only bifurcations)
        cumulative_node_count += 1
        node_index += 1
        left_label = "nd" + str(cumulative_node_count)
        left_child = dp.Node(
            taxon=dp.Taxon(label=left_label),
            label=left_label,
            edge_length=0.0
        )
        left_child.is_sa = False
        left_child.is_sa_lineage = False
        left_child.is_sa_dummy_parent = False
        left_child.alive = True
        left_child.sampled = True
        left_child.state = left_arriving_state
        left_child.annotations.add_bound_attribute("state")
        left_child.index = node_index
        left_child.annotations.add_bound_attribute("index")
        state_representation_dict[left_child.state].add(
            left_child.label
        )

        cumulative_node_count += 1
        node_index += 1
        right_label = "nd" + str(cumulative_node_count)
        right_child = dp.Node(
            taxon=dp.Taxon(label=right_label),
            label=right_label,
            edge_length=0.0
        )
        right_child.is_sa = False
        right_child.is_sa_lineage = False
        right_child.is_sa_dummy_parent = False
        right_child.alive = True
        right_child.sampled = True
        right_child.state = right_arriving_state
        right_child.annotations.add_bound_attribute("state")
        right_child.index = node_index
        right_child.annotations.add_bound_attribute("index")
        state_representation_dict[right_child.state].add(
            right_child.label
        )

        ####################################
        # Adding state transition instance #
        # to class members that bookkeep   #
        # for either painting mid-branch,  #
        # or for stochastic map I/O        #
        ####################################
        for idx, child_node in enumerate((left_child, right_child)):
            arriving_state: int = -1

            if idx == 0:
                arriving_state = left_arriving_state

            elif idx == 1:
                arriving_state = right_arriving_state

            if arriving_state != chosen_node.state:
                state_trans = AttributeTransition(
                    "state",
                    child_node.label,
                    event_t,
                    chosen_node.state,
                    arriving_state,
                    at_speciation=True
                )

                if child_node.label in state_transition_dict:
                    state_transition_dict[child_node.label] \
                        .append(state_trans)

                else:
                    state_transition_dict[child_node.label] = [state_trans]

        if chosen_node.state != left_arriving_state or \
                chosen_node.state != right_arriving_state:
            clado_state_trans = AttributeTransition(
                "state",
                chosen_node.label,
                event_t,
                chosen_node.state,
                left_arriving_state,
                to_state2=right_arriving_state,
                at_speciation=True
            )

            clado_state_transition_dict[chosen_node.label] = clado_state_trans

        # updating parent node
        # (1) adding both children
        chosen_node.add_child(left_child)
        tr_namespace.add_taxon(left_child.taxon)
        chosen_node.add_child(right_child)
        tr_namespace.add_taxon(right_child.taxon)

        # (2) parent node does not represent its state
        # and cannot be targeted anymore
        try:
            # state of parent is now irrelevant
            state_representation_dict[chosen_node.state] \
                .remove(chosen_node.label)

        except Exception as e:
            # if chosen_node is root, it has never been added
            # to state_representation_dict, so we pass
            pass

        chosen_node.alive = False
        untargetable_node_set \
            .add(chosen_node.label)  # cannot pick parent anymore!

        # (3) mark parent node as the last node to undergo event
        #     will use this node to undo last event when number of species
        #     go beyond maximum
        last_node2speciate = chosen_node

        # (4) if chosen node was on a lineage with SAs, we update the SAs info
        if chosen_node.is_sa_lineage:
            self._update_sa_lineage_dict(
                event_t,
                sa_lineage_dict,
                [chosen_node.label],
                debug=debug)

        return last_node2speciate, cumulative_node_count, node_index

    def _execute_death(
            self,
            tr_namespace: dp.TaxonNamespace,
            chosen_node: dp.Node,
            state_representation_dict: ty.Dict[int, ty.Set[str]],
            sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]],
            untargetable_node_set: ty.Set[str],
            event_t: float,
            debug=False) -> dp.Node:
        """Execute lineage death.

        This method has both side-effects and a return. The
        side-effects include the updating of class members. These
        members are passed as arguments in the signature so that we
        can have a better idea of what to expect as behavior:
        
            (i)   self.tr.tr_namespace,
            (ii)  self.state_representation_dict
            (iii) self.sa_lineage_dict
            (iv)  self.untargetable_node_set

        Args:
            tr_namespace (dendropy.TaxonNamespace): Dendropy object
                recording taxa in the tree.
            chosen_node (dendropy.Node): Node that will undergo
                extinction.
            state_representation_dict (dict): Dictionary that keeps
                track of all states currently represented by living
                lineages. Keys are integer representing the states,
                values are sets of taxon names (strings).
            sa_lineage_dict (dict): Dictionary that keeps track of
                nodes that have direct (sampled) ancestor children.
                Keys are taxon names (strings), values are lists
                of SampledAncestor objects.
            untargetable_node_set (str): Set of Node labels that
                cannot be targeted for events anymore (went extinct).
            event_t (float): Time of death event taking place.
            debug (bool): Flag for printing debugging messages.
                If 'True', prints messages. Defaults to 'False'.

        Returns:
            (dendropy.Node) Last node to die.
        """

        if debug:
            print("EXTINCTION of node " + chosen_node.label)

        # with this node gone, this state
        # is not represented by it anymore
        state_representation_dict[chosen_node.state].remove(chosen_node.label)

        # we also cannot choose this extinct node
        # anymore to undergo an event
        chosen_node.alive = False
        chosen_node.sampled = False
        untargetable_node_set.add(chosen_node.label)

        # if chosen node was on a lineage with SAs,
        # we update the SAs info
        if chosen_node.is_sa_lineage:
            self._update_sa_lineage_dict(event_t,
                                         sa_lineage_dict,
                                         [chosen_node.label])

        ######################################
        # Special case: origin went extinct, #
        # we slap a brosc node               #
        ######################################
        if chosen_node.label == "origin":
            # at this point, the origin was chosen to die, but this node will
            # never be extended (the origin always has an
            # origin_edge_length = 0.0); for us to account for the evolution
            # (branch length) that has happened before this death -- between
            # the origin and the brosc_node being added -- we must add
            # 'event_t' as the brosc_node edge_length (other nodes are always
            # added with edge_length = 0.0, and have their edges extended
            # when a new event takes place)
            brosc_node = dp.Node(taxon=dp.Taxon(label="brosc"),
                                 label="brosc",
                                 edge_length=event_t)
            brosc_node.state = chosen_node.state
            brosc_node.alive = False
            brosc_node.sampled = False
            brosc_node.is_sa = False
            brosc_node.is_sa_dummy_parent = False
            brosc_node.is_sa_lineage = False
            untargetable_node_set.add(brosc_node.label)
            chosen_node.add_child(brosc_node)
            tr_namespace.add_taxon(brosc_node.taxon)

        last_node2die = chosen_node

        return last_node2die

    def _execute_anatrans(
            self,
            tr_namespace: dp.TaxonNamespace,
            chosen_node: dp.Node,
            state_representation_dict: ty.Dict[int, ty.Set[str]],
            state_transition_dict: ty.Dict[str, ty.List[AttributeTransition]],
            untargetable_node_set: ty.Set[str],
            sse_anatrans_rate_object: sseobj.DiscreteStateDependentRate,
            event_t: float,
            debug: bool = False) -> None:
        """
        Execute anagenetic trait-state transition on path to chosen
        node.

        This method has only side-effects. These include the updating
        of class members, which are passed as arguments in the
        signature so that we can have a better idea of what to expect
        as behavior:
        
            (i)   self.tr.tr_namespace,
            (ii)  self.state_representation_dict
            (iii) self.state_transition_dict
            (iv)  self.untargetable_node_set

        Args:
            tr_namespace (dendropy.TaxonNamespace): Dendropy object
                recording taxa in the tree.
            chosen_node (dendropy.Node): Node that will undergo
                anagenetic trait-state transition.
            state_representation_dict (dict): Dictionary that keeps
                track of all states currently represented by living
                lineages. Keys are integer representing the states,
                values are sets of taxon names (strings).
            state_transition_dict (dict): Dictionary that keeps track
                of nodes subtending which state transitions happen
                (used for plotting). Keys are taxon names (strings)
                and values are lists of AttributeTransition objects.
            untargetable_node_set (set): Set of Node labels that
                cannot be targeted for events anymore (went extinct).
            sse_anatrans_rate_object (DiscreteStateDependentRate): SSE
                rate parameter object holding the departing/arriving
                states.
            event_t (float): Time of anagenetic transition event taking
                place.
            debug (bool): If 'True', prints debugging messages.
                Defaults to False.
        """

        if debug:
            print("TRANSITION of node " + chosen_node.label \
                  + " from state " + str(chosen_node.state) \
                  + " to state " \
                  + str(sse_anatrans_rate_object.arriving_states[0]))

        # new state
        arriving_state = sse_anatrans_rate_object.arriving_states[0]

        # old state is not represented anymore
        state_representation_dict[chosen_node.state].remove(chosen_node.label)

        ########################################################
        # Special case: anagenetic trasition event before root #
        ########################################################
        if chosen_node.label == "origin":
            # at this point, the origin was chosen to undergo anagenetic
            # state transition, but this node will never become extended
            # (the origin always has an origin_edge_length = 0.0)
            #
            # for us to account for the evolution (branch length) that
            # has happened before this death -- between the origin and
            # the brosc_node being added -- we must add event_t as the
            # brosc_node edge_length (other nodes are always added with
            # edge_length = 0.0, and have their edges extended when a new
            # event takes place)
            brosc_node = dp.Node(taxon=dp.Taxon(label="brosc"),
                                 label="brosc",
                                 edge_length=event_t)
            brosc_node.state = chosen_node.state
            brosc_node.alive = True
            brosc_node.sampled = True
            brosc_node.is_sa = False
            brosc_node.is_sa_dummy_parent = False
            brosc_node.is_sa_lineage = False
            tr_namespace.add_taxon(brosc_node.taxon)

            # origin cannot be selected anymore
            chosen_node.alive = False  # origin node is no longer alive
            untargetable_node_set.add(chosen_node.label)  # cannot be targeted
            chosen_node.add_child(brosc_node)  # brosc is now child of origin

            # now make chosen_node become brosc, so the rest
            # of the anagenetic transition event can be executed
            chosen_node = brosc_node

        ####################################
        # Adding state transition instance #
        # to class member that bookkeeps   #
        # (for painting mid-branch)        #
        ####################################
        state_trans = AttributeTransition(
            "state",
            chosen_node.label,
            event_t,
            chosen_node.state,
            arriving_state)

        try:
            state_transition_dict[chosen_node.label].append(state_trans)

        except KeyError:
            state_transition_dict[chosen_node.label] = [state_trans]

        # last updates and bookkeping #
        # new state gets assigned
        chosen_node.state = arriving_state

        # new state is now represented
        state_representation_dict[chosen_node.state].add(chosen_node.label)

    def _execute_sample_ancestor(
            self,
            tr_namespace: dp.TaxonNamespace,
            chosen_node: dp.Node,
            state_representation_dict: ty.Dict[int, ty.Set[str]],
            sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]],
            untargetable_node_set: ty.Set[str],
            cumulative_sa_count: int,
            node_index: int,
            event_t: float,
            debug: bool = False) -> ty.Tuple[int]:
        """Execute sampling of direct lineage ancestor.

        This method has both side-effects and a return. The
        side-effects include the updating of class members. These
        members are passed as arguments in the signature so that we
        can have a better idea of what to expect as behavior:
        
            (i)   self.tr.tr_namespace,
            (ii)  self.state_representation_dict
            (iii) self.sa_lineage_dict
            (iv)  self.untargetable_node_set

        Args:
            tr_namespace (dendropy.TaxonNamespace): Dendropy object
                recording taxa in the tree.
            chosen_node (dendropy.Node): Node that will undergo the
                direct ancestor sampling.
            state_representation_dict (dict): Dictionary that keeps
                track of all states currently represented by living
                lineages. Keys are integer representing the states,
                values are sets of taxon names (strings).
            sa_lineage_dict (dict): Dictionary that keeps track of
                nodes that have direct (sampled) ancestor children.
                Keys are taxon names (strings) and values are lists
                of SampledAncestor objects.
            untargetable_node_set (str): Set of Node labels that
                cannot be targeted for events anymore (went extinct).
            cumulative_sa_count (int): Total number of direct ancestor
                nodes in the tree (to be used in labeling).
            node_index (int): Total node count to be used to number
                nodes (used in their 'index' attribute).
            event_t (float): Time of ancestor sampling event taking
                place.
            debug (bool): If 'True', prints debugging messages.
                Defaults to False.

        Returns:
            (int): Cumulative number of direct (sampled) ancestors.
        """

        if debug:
            print("> ANCESTOR-SAMPLING of node " + chosen_node.label \
                  + " , keeping state " + str(chosen_node.state))

        # modeling sampled ancestor as bifurcating event, one child is the same as parent lineage,
        # second child is a terminal node that will never be elongated (branch length = 0.0);
        # the original parent is made into a dummy node

        if chosen_node.label == "origin":
            # ancestor sampling event took place before root was born so we create a replacement
            # for the root, the 'brosc_node' before-root-origin-single-child; this is the node that
            # will undergo the ancestor sampling event
            #
            # the origin was chosen to undergo the ancestor sampling, but this node will never become
            # extended (the origin always has an origin_edge_length = 0.0); for us to account for the
            # evolution (branch length) that has happened before this death -- between the origin and
            # the brosc node being added -- we must add event_t as the brosc node edge_length (other
            # nodes are always added with edge_length = 0.0, and have their edges extended when a new
            # event takes place)
            node_index += 1
            brosc_node = dp.Node(taxon=dp.Taxon(label="brosc"),
                                 label="brosc",
                                 edge_length=event_t)
            brosc_node.state = chosen_node.state
            brosc_node.annotations.add_bound_attribute("state")
            brosc_node.index = node_index
            brosc_node.annotations.add_bound_attribute("index")
            brosc_node.alive = True
            brosc_node.sampled = True
            brosc_node.is_sa = False
            brosc_node.is_sa_dummy_parent = False
            brosc_node.is_sa_lineage = True
            # this node is alive and can be selected
            state_representation_dict[chosen_node.state].add(brosc_node.label)
            tr_namespace.add_taxon(brosc_node.taxon)

            # updating origin
            chosen_node.alive = False  # origin is no longer alive
            untargetable_node_set.add("origin")  # cannot be targeted
            state_representation_dict[chosen_node.state] \
                .remove(chosen_node.label)  # and does not represent state
            chosen_node.add_child(brosc_node)

            # we will execute the ancestor sampling on this new chosen_node
            chosen_node = brosc_node

            # print("\n\nMust have a sampled ancestor before the root!")

        parent_label = chosen_node.label  # (parent's label will change below)
        parent_state = chosen_node.state

        # doing lineage that remains alive
        left_child = dp.Node(taxon=dp.Taxon(label=parent_label),
                             label=parent_label,
                             edge_length=0.0)
        left_child.is_sa = False
        left_child.is_sa_lineage = True
        left_child.is_sa_dummy_parent = False
        left_child.alive = True
        left_child.sampled = True
        left_child.state = chosen_node.state  # gets parent's state
        left_child.annotations.add_bound_attribute("state")
        # we don't update state_representation_dict,
        # because parent's label is already there

        # doing sampled ancestor
        cumulative_sa_count += 1
        node_index += 1
        right_label = "sa" + str(cumulative_sa_count)
        right_child = dp.Node(taxon=dp.Taxon(label=right_label),
                              label=right_label,
                              edge_length=0.0)
        right_child.is_sa = True
        right_child.is_sa_lineage = False
        right_child.is_sa_dummy_parent = False
        right_child.alive = False
        right_child.sampled = False
        right_child.state = chosen_node.state  # gets parent's state
        right_child.annotations.add_bound_attribute("state")
        right_child.index = node_index
        right_child.annotations.add_bound_attribute("index")
        # state_representation_dict[right_child.state].add(right_child.label)
        # I don't think sampled ancestor represents state
        # the way we're implementing sampled ancestors here,
        # they cannot be picked for events; it is their sister
        # lineage that can be picked (their parent is becomes a dummy node)

        # cannot pick sampled ancestor to undergo
        # event the way its implemented
        untargetable_node_set.add(right_child.label)

        ####################################################
        # Adding sampled ancestor and its lineage node     #
        # to class member that stashes it, or initializing #
        # class member                                     #
        ####################################################
        sa = SampledAncestor(right_label, parent_label, event_t, parent_state)
        try:
            sa_lineage_dict[parent_label].append(sa)

        except KeyError:
            sa_lineage_dict[parent_label] = [sa]

        # updating parent node
        # (1) adding both children
        chosen_node.add_child(left_child)
        tr_namespace.add_taxon(left_child.taxon)
        chosen_node.add_child(right_child)
        tr_namespace.add_taxon(right_child.taxon)
        # (2) update parent's label
        chosen_node.label = "dummy" + str(cumulative_sa_count)
        # note that here, we do not remove parent's label from targeted nodes,
        # because left child will now have that label and should be targetable
        # (3) parent label cannot be targeted anymore
        chosen_node.alive = False
        chosen_node.sampled = False
        chosen_node.is_sa_dummy_parent = True
        # left child will always be the sa_lineage
        chosen_node.is_sa_lineage = False
        # cannot pick parent!
        untargetable_node_set.add(chosen_node.label)

        return cumulative_sa_count, node_index

    def _execute_event(
            self,
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
            debug: bool = False) -> ty.Tuple[dp.Node, int]:
        """Execute event on chosen node.

        This method calls the correct event executer depending on which
        discrete SSE rate was passed as argument. That executer will
        have side-effects and/or returns (the latter will be captured
        and returned).

        Args:
            tr_namespace (dendropy.TaxonNamespace): DendroPy object
                recording taxa in the tree.
            sse_rate_object (DiscreteStateDependentRate): SSE rate
                parameter object holding the departing/arriving states.
            chosen_node (dendropy.Node): Node that will undergo
                birth event.
            state_representation_dict (dict): Dictionary that keeps
                track of all states currently represented by living
                lineages. Keys are integer representing the states,
                values are sets of taxon names (strings).
            state_transition_dict (dict): Dictionary that keeps track
                of nodes subtending which state transitions happen
                (used for plotting). Keys are taxon names (strings)
                and values are lists of AttributeTransition objects.
            clado_state_transition_dict (dict): Dictionary that keeps
                track of speciating nodes whose children have different
                states (used for stochastic maps I/O). Keys are
                taxon names (strings) and values are
                AttributeTransition objects.
            sa_lineage_dict (dict): Dictionary that keeps track of
                nodes that have direct (sampled) ancestor children.
                Keys are taxon names (strings), values are lists
                of SampledAncestor objects.
            untargetable_node_set (set): Set of Node labels that
                cannot be targeted for events anymore (went extinct).
            cumulative_node_count (int): Total number of non-direct
                ancestor nodes in the tree (to be used in labeling).
            cumulative_sa_count (int): Total number of direct ancestor
                nodes in the tree (to be used in labeling).
            node_index (int): Total node count to be used to number
                nodes (used in their 'index' attribute).
            event_t (float): Time of event taking place.
            debug (bool): If 'True', prints debugging messages.
                Defaults to False.

        Returns:
            (tuple): Tuple with four objects, the last node to
            speciate (dendropy.Node), the tree's non-direct ancestor
            cumulative node count (int), cumulative direct (sampled)
            ancestor count (int), and total node count (int).
        """

        macroevol_event = sse_rate_object.event

        # unlike in C++ (which has 'block-level' scoping), we do not need to
        # declare/initialize variables in Python if used in if/for/while
        # blocks -- Python has 'function-level' scoping.
        #
        # however, the design in PJ is to pretend Python has block-level
        # scoping, so:
        last_chosen_node: dp.Node = None

        if macroevol_event == sseobj.MacroevolEvent.W_SPECIATION or \
           macroevol_event == sseobj.MacroevolEvent.BW_SPECIATION or \
                macroevol_event == sseobj.MacroevolEvent.ASYM_SPECIATION:

            last_chosen_node, cumulative_node_count, node_index = \
                self._execute_birth(
                    tr_namespace,
                    chosen_node,
                    state_representation_dict,
                    sa_lineage_dict,
                    state_transition_dict,
                    clado_state_transition_dict,
                    untargetable_node_set,
                    cumulative_node_count,
                    node_index,
                    sse_rate_object,
                    event_t,
                    debug=debug)

        elif macroevol_event == sseobj.MacroevolEvent.EXTINCTION:
            self._execute_death(
                tr_namespace,
                chosen_node,
                state_representation_dict,
                sa_lineage_dict,
                untargetable_node_set,
                event_t,
                debug=debug)

        elif macroevol_event == sseobj.MacroevolEvent.ANAGENETIC_TRANSITION:
            self._execute_anatrans(
                tr_namespace,
                chosen_node,
                state_representation_dict,
                state_transition_dict,
                untargetable_node_set,
                sse_rate_object,
                event_t,
                debug=debug)

        elif macroevol_event == sseobj.MacroevolEvent.ANCESTOR_SAMPLING:
            cumulative_sa_count, node_index = \
                self._execute_sample_ancestor(
                    tr_namespace,
                    chosen_node,
                    state_representation_dict,
                    sa_lineage_dict,
                    untargetable_node_set,
                    cumulative_sa_count,
                    node_index,
                    event_t,
                    debug=debug)

        return last_chosen_node, cumulative_node_count, cumulative_sa_count, node_index

    def _annotate_sampled(
            self,
            a_time: float,
            living_nodes: ty.List[dp.Node],
            sample_idx: int,
            single_node: ty.Optional[dp.Node] = None) -> None:
        """
        Annotate each living node as sampled or not.
        
        This method is called whenever a maximum-time stop condition
        has been met. (It is not called when the maximum taxon count is
        met because bookkeeping which taxa get sampled as the tree
        grows it more complicated; currently this class does not
        support it.)

        This method has the side-effect of setting node attribute
        .sampled to 'True' or 'False' depending on:
            (i)  the sampling probability of its associated state, and
            (ii) the time slice (epoch) that node lives in.

        By default, all events but 'death' set .sampled to True for
        terminal nodes, and to False for internal nodes.

        Args:
            a_time (float): The time at which the simulation stopped
                (i.e., the 'stop_value' provided to the distribution,
                or the tree height).
            living_nodes (dp.Node): List of DendroPy nodes.
            sample_idx (int): Index of current sample.
        """

        if single_node is None:
            for living_nd in living_nodes:
                st = living_nd.state
                is_sampled = self.sse_stash.get_prob_handler() \
                    .randomly_decide_taxon_sampling_at_time_at_state(
                        a_time, st, sample_idx)

                living_nd.sampled = is_sampled
        
        else:
            st = single_node.state
            is_sampled = self.sse_stash.get_prob_handler() \
                .randomly_decide_taxon_sampling_at_time_at_state(
                        a_time, st, sample_idx)
            single_node.sampled = is_sampled

    ##########################
    # Main simulation method #
    ##########################

    def simulate(self,
                 a_start_state: int,
                 a_stop_value: ty.Union[int, float],
                 sample_idx: int = 0) -> AnnotatedTree:
        """Sample (simulate) discrete-SSE tree.

        The simulated tree may or not meet the conditions specified by
        the user. Conditioning is handled outside of this method (see
        documentation for .generate()).

        Args:
            a_start_state (int): State at seed node (origin or root).
            a_stop_value (float): Value to stop simulation with (number
                of tips or tree height).
            sample_idx (int): Index of current sample.

        Returns:
            (AnnotatedTree): A tree annotated with discrete states.
        """

        def _get_next_event_time(total_rate: float) -> float:
            """Draw next exponentially distributed event time.

            Args:
                current_node_target_count (int): Number of nodes that
                    can experience an event.
                total_rate (float): Sum of all rates currently represented
                    by nodes that can experience an event.

            Returns:
                (float): Time to next event.
            """

            next_time = \
                float(
                    dnpar.DnExponential.draw_exp(1,
                                                total_rate,
                                                rate_parameterization=True)[0])

            return next_time
    
        def _extend_all_living_nodes(branch_length: float) -> None:
            """Extend branches subtending all nodes by some length.

            This method has only side-effects.

            Args:
                branch_length (float): The amount by which to extend
                    the branches subtending all nodes.
            """
            for nd in tr:
                # the root node is extended (root_edge is made > 0.0)
                # here too, when the root is picked
                if nd.label not in untargetable_node_set and \
                        nd.label != "origin" and nd.alive:
                    nd.edge_length += branch_length
        
        ############################################
        # START initializing values pre-simulation #
        #
        # start variables
        start_state = a_start_state

        # stop variables
        reached_stop_condition: bool = False
        t_stop: float = 0.0  # time, not age
        max_obs_nodes: int = 0
        if self.stop == "age" and isinstance(a_stop_value, float):
            t_stop = a_stop_value

        elif self.stop == "size" and isinstance(a_stop_value, int):
            max_obs_nodes = a_stop_value

        # variables for during simulation
        latest_t: float = 0.0  # will increase as simulation progresses
        last_chosen_node: dp.Node = None
        current_node_target_count: int = 1
        cumulative_node_count: int = 1  # not being used atm
        cumulative_sa_count: int = 0  # not being used atm
        untargetable_node_set = set()  # not being used atm (.alive is used)
        tree_invalid: bool = False
        tree_died: bool = False

        # bookkeeping members
        # { 0: ["origin", "root"], 1: ["nd1", "nd2"], 2:["nd3",...], ... }
        state_representation_dict: ty.Dict[int, ty.Set[str]] = \
            dict((i, set()) for i in range(self.events.state_count))
        # these last two are for plotting
        state_transition_dict: \
            ty.Dict[str, ty.List[AttributeTransition]] = dict()
        sa_lineage_dict: \
            ty.Dict[str, ty.List[SampledAncestor]] = dict()
        # this dict is for stochastic maps I/O
        clado_state_transition_dict: \
            ty.Dict[str, AttributeTransition] = dict()

        ######################
        # Germinate the tree #
        ######################
        tr: dp.Tree = None
        brosc_node: dp.Node = None # will be None if tree starts at root
        root_is_born: bool = False
        tr, brosc_node, cumulative_node_count, root_is_born = \
            self._germinate_tree(start_state,
                                 self.with_origin,
                                 state_representation_dict,
                                 untargetable_node_set)
        self.root_is_born = root_is_born
        node_index: int = cumulative_node_count

        # need these counts initialized for sampling first event
        living_nodes = [nd for nd in tr if nd.alive]
        current_node_target_count = len([nd for nd in tr if nd.alive])

        # info for printing if user asks for it
        size_set: set(int) = set([])

        # END initializing values pre-simulation #
        ##########################################

        # START simulation loop #
        time_slice_idx: int = 0
        while (time_slice_idx < self.n_time_slices and \
               not reached_stop_condition):

            # checking and printing tree size info #
            if self.info:
                tr_size = len(tr)
                
                if tr_size % 100 == 0 and tr_size not in size_set:
                    size_set.add(tr_size)
                    print("Tree size =", str(len(tr)))

            # (1) Calculate total overall rate
            # 
            # overall rate depends on the states of current targetable
            # nodes
            rate_for_exponential_distn, \
                state_unweighed_rate_sums = \
                self.events.total_rate(latest_t,
                                       state_representation_dict,
                                       value_idx=sample_idx,
                                       debug=self.debug)
            # [1] are all states with >= 1 representation

            # (2) get the time to the next event
            t_to_next_event = \
                _get_next_event_time(rate_for_exponential_distn)
            latest_t += t_to_next_event

            # (3) get end time of this time slice (epoch) #
            # 
            # user provides it as end ages, but self.events
            # (MacroevolEventHandler) converts it to time slice ends
            # upon initialization
            next_max_t: float
            t_end = self.slice_t_ends[time_slice_idx]
            if isinstance(t_end, float):
                next_max_t = t_end

            # (4) check if new event time cut through the end of a #
            # time slice 
            # 
            # NOTE: this step executes if there is more than 1 time
            # slice (epoch) and the stop condition is maximum age
            #
            # next_max_t will be None if self.slice_t_ends is empty
            excess_t = 0.0
            if self.stop == "age" and self.n_time_slices > 1 \
                    and latest_t > next_max_t:
                excess_t = latest_t - next_max_t
                latest_t = next_max_t
                time_slice_idx += 1

                # extend all lineages (could keep just the extend in the else-block
                # below) and keep everything up-to-date
                _extend_all_living_nodes(t_to_next_event - excess_t)

                # STOP CHECK: tree grew beyond its maximum time duration #
                #
                # this last event goes beyond _next_max_t, which, if t_stop
                # means we should stop this process
                if next_max_t == t_stop:
                    # updates SA info for plotting
                    sa_lineage_node_labels = \
                        [nd.label for nd in living_nodes if nd.is_sa_lineage]
                    self._update_sa_lineage_dict(
                        t_stop,
                        sa_lineage_dict,
                        sa_lineage_node_labels)

                    # if origin is the only node (root always has children),
                    # we tie off the process with a 'finish' node (the brosc node)
                    if self.with_origin and tr.seed_node.alive and \
                            len(tr.seed_node.child_nodes()) == 0:
                        # we make origin edge length the max age of the tree
                        brosc_node.edge_length = t_stop
                        brosc_node.state = tr.seed_node.state
                        brosc_node.alive = True  # .sampled=True upon init
                        self._annotate_sampled(t_stop,
                                               [],
                                               sample_idx,
                                               brosc_node)
                        tr.seed_node.add_child(brosc_node)
                        tr.taxon_namespace.add_taxon(brosc_node.taxon)

                        # rho-refactor
                        tr.seed_node.alive = False
                        tr.seed_node.sampled = False

                    # annotating .sampled (just brosc)
                    living_nodes = [nd for nd in tr if nd.alive]
                    self._annotate_sampled(t_stop, living_nodes, sample_idx)

                    reached_stop_condition = True

                    break

                # goes back to top of while-loop restarting at next time slice
                continue

            # (5) draw event and execute it #
            # 
            # new event time did not go over a time slice
            else:
                # STOP CHECK: tree grew beyond its maximum time duration #
                #
                # we check here as well as the if-block in (4) above because
                # that one only executes if there is more than one time slice
                # (epoch)
                if (self.stop == "age" and (latest_t > t_stop)) or \
                        latest_t < 0.0:
                    _extend_all_living_nodes(
                        t_stop - (latest_t - t_to_next_event))

                    # updates SA info for plotting
                    sa_lineage_node_labels = \
                        [nd.label for nd in living_nodes if nd.is_sa_lineage]
                    self._update_sa_lineage_dict(
                        t_stop,
                        sa_lineage_dict,
                        sa_lineage_node_labels)

                    # if origin is the only node (root always has children)
                    # we tie off the process with a 'finish' node (the brosc node)
                    if self.with_origin and tr.seed_node.alive and \
                            len(tr.seed_node.child_nodes()) == 0:
                        brosc_node.edge_length = t_stop
                        brosc_node.state = tr.seed_node.state
                        brosc_node.alive = True  # sampled=True upon init
                        self._annotate_sampled(t_stop,
                                               [],
                                               sample_idx,
                                               brosc_node)
                        tr.seed_node.add_child(brosc_node)
                        tr.taxon_namespace.add_taxon(brosc_node.taxon)

                        # rho-refactor
                        tr.seed_node.alive = False
                        tr.seed_node.sampled = False

                    # annotating .sampled (just brosc)
                    living_nodes = [nd for nd in tr if nd.alive]
                    self._annotate_sampled(t_stop, living_nodes, sample_idx)

                    reached_stop_condition = True

                    break

                # elongate all lineages
                # original implmn w/o slices
                _extend_all_living_nodes(t_to_next_event)

                # (6) draw a node we'll apply the event to #
                #
                # a lineage will be chosen in proportion to the total rate
                # of its state
                lineage_weights = \
                    [state_unweighed_rate_sums[nd.state] \
                     for nd in living_nodes]
                chosen_node = \
                    random.choices(living_nodes, weights=lineage_weights)[0]

                # debugging
                if self.debug:
                    print("\nChosen node: " + chosen_node.label +
                          + " in state " + str(chosen_node.state) + "\n")

                # (7) draw an event
                macroevol_event_state = chosen_node.state
                
                # examples/geosse_timehet.pj: if it prints, we have forbidden events
                # happening in the oldest epoch!
                # if latest_t < 1.0 and macroevol_event_state != 0:
                #     print("t = " + str(latest_t) + \
                #           ", state of node = " + str(macroevol_event_state))
                
                # denominator for sampling event (normalization factor is the
                # sum of rates departing from the chosen node's state)
                state_conditioned_rate_sum = \
                    state_unweighed_rate_sums[macroevol_event_state]
                
                sse_rate_param_in_list = \
                    self.events.sample_event_sse_rate_param(
                        state_conditioned_rate_sum,
                        latest_t,
                        [macroevol_event_state],
                        value_idx=sample_idx,
                        debug=self.debug)
                # total_rate and latest_t have been updated at this point
                sse_rate_param = sse_rate_param_in_list[0]

                # (8) execute event
                # TODO: later add the living_nodes list as a return
                # of execute_event, whenever execute_birth and execute_eath
                # are called, so we do not need to traverse the tree and
                # update living_nodes all the time (below, in step (9))
                # print("dn_discrete_sse.py: at (5.3)")
                (last_chosen_node,
                 cumulative_node_count,
                 cumulative_sa_count,
                 node_index) \
                    = self._execute_event(
                        tr.taxon_namespace,
                        sse_rate_param,
                        chosen_node,
                        state_representation_dict,
                        state_transition_dict,
                        clado_state_transition_dict,
                        sa_lineage_dict,
                        untargetable_node_set,
                        cumulative_node_count,
                        cumulative_sa_count,
                        node_index,
                        latest_t,
                        debug=self.debug
                    )

                # (9) update number of lineages after event
                #
                # TODO: at some point, make 'living_nodes' a list that is
                # passed to all execute_birth and execute_death functions
                # so that we don't need to traverse the tree every event
                # to list living nodes that can be  targeted
                #
                # what counts for targetable node here are living
                # terminal nodes
                living_nodes = [nd for nd in tr if nd.alive]
                current_node_target_count = len(living_nodes)

                # (10) check for stop conditions
                #
                # check if all lineages went extinct, stop!
                if current_node_target_count == 0:
                    tree_died = True

                    # updates SA info for plotting
                    sa_lineage_node_labels = \
                        [nd.label for nd in living_nodes
                         if nd.is_sa_lineage]
                    self._update_sa_lineage_dict(
                        latest_t,
                        sa_lineage_dict,
                        sa_lineage_node_labels)

                    reached_stop_condition = True

                    break

                # [Condition is "age", but rejection sampling
                # was requested for taxon count range]
                #
                # this helps stop trees growing out of control
                # within age limit
                if current_node_target_count == self.abort_at_alive_count:
                    tree_invalid = True
                    reached_stop_condition = True

                    break

                # STOP CHECK: tree grew beyond its max taxon count #
                #
                # it's '>' max_obs_nodes below, instead of '==' or '>='
                # because we want to have one extra species (that then gets removed),
                # i.e., we want the longest possible tree up to max_obs_nodes
                # by waiting until one-too-many speciation events
                if self.stop == "size" and \
                        current_node_target_count > max_obs_nodes:
                    # delete last created children
                    last_chosen_node.clear_child_nodes()

                    # parent is back alive
                    last_chosen_node.alive = True

                    # we must add parent back so it can be extended
                    untargetable_node_set.remove(last_chosen_node.label)

                    if max_obs_nodes == 1:
                        # special case: birth happens on just-added root, and
                        # we want just a single lineage; only makes sense if
                        # starting at origin, b/c if starting at root, two
                        # lineages already exist)
                        if self.with_origin:
                            # getting root attributes
                            brosc_node.state = last_chosen_node.state
                            brosc_node.edge_length = \
                                last_chosen_node.edge_length
                            # should be 'True'
                            brosc_node.alive = last_chosen_node.alive
                            # connecting brosc to origin
                            last_chosen_node.parent_node.add_child(brosc_node)
                            tr.taxon_namespace.add_taxon(brosc_node.taxon)
                            # removing root and replacing it with brosc
                            last_chosen_node.parent_node \
                                .remove_child(last_chosen_node)
                            last_chosen_node = brosc_node

                    # won't use it again, but to be safe
                    state_representation_dict[last_chosen_node.state] \
                        .add(last_chosen_node.label)

                    # updates SA info for plotting
                    sa_lineage_node_labels = \
                        [nd.label for nd in tr.leaf_node_iter() \
                         if nd.is_sa_lineage]
                    self._update_sa_lineage_dict(t_stop,
                                                 sa_lineage_dict,
                                                 sa_lineage_node_labels)

                    reached_stop_condition = True

                    break

        # (11) got out of while loop because met stop condition
        # 'at' is scoped to simulate_a_tree() function
        if self.stop == "age":
            # print("dn_discrete_sse.py: at (11)")
            at = AnnotatedTree(
                tr,
                self.events.state_count,
                start_at_origin=self.with_origin,
                max_age=a_stop_value,
                slice_t_ends=self.slice_t_ends,
                slice_age_ends=self.events.slice_age_ends,
                sa_lineage_dict=sa_lineage_dict,
                at_dict=state_transition_dict,
                clado_at_dict=clado_state_transition_dict,
                condition_on_obs_both_sides_root=\
                    self.condition_on_obs_both_sides_root,
                tree_invalid=tree_invalid,
                tree_died=tree_died)

        elif self.stop == "size":
            at = AnnotatedTree(
                tr,
                self.events.state_count,
                start_at_origin=self.with_origin,
                sa_lineage_dict=sa_lineage_dict,
                at_dict=state_transition_dict,
                clado_at_dict=clado_state_transition_dict,
                condition_on_obs_both_sides_root=\
                    self.condition_on_obs_both_sides_root,
                tree_invalid=tree_invalid)

        if self.debug:
            print(at.tree.as_string(schema="newick",
                                    suppress_annotations=False,
                                    suppress_internal_taxon_labels=True))

        return at
        # END simulation loop #

    def _is_tr_ok(
            self,
            ann_tr: AnnotatedTree,
            a_stop_value: ty.Union[int, float]) -> bool:
        """Check if simulated tree meets stop conditions.

        Args:
            a_tree (AnnotatedTree): A tree with discrete-state
                annotation.
            a_stop_value (Union[int, float]): Value specified for
                stopping simulation.

        Returns:
            (bool): Boolean specifying if tree meets stop condition.
        """

        ##############################
        # STOP: conditioning not met #
        ##############################
        # An aborted trajectory is never a completed draw, regardless of conditioning.
        if ann_tr.tree_invalid:
            if self.info:
                print("Rejected tree, it grew out of control!")
            return False

        if self.condition_on_survival and ann_tr.tree_died:
            return False

        if self.condition_on_speciation and not \
                isinstance(ann_tr.root_node, dp.Node):
            return False

        # need to worry about rejection sampling
        # when user conditions reconstructed tree
        if self.condition_on_obs_both_sides_root:
            # condition on speciation means we want a root node!
            # (one could still have internal nodes because of direct ancestor
            # sampling, which is why we test this condition here)
            if self.condition_on_speciation and not \
                    isinstance(ann_tr.root_node, dp.Node):
                return False

            # reconstructed tree is empty, i.e., maybe there was a root, but
            # everyone died
            rec_tree = ann_tr.extract_reconstructed_tree()
            if rec_tree.__str__() == "":
                return False

        ###################
        # STOP: Tree size #
        ###################
        if self.stop == "size":
            # we want a root node and that it has children (alive or not)
            if self.condition_on_speciation and \
                    isinstance(ann_tr.root_node, dp.Node) and \
                    len(ann_tr.root_node.child_nodes()) == 0:
                return False

            # we want a specific number of living extant taxa
            if self.condition_on_survival and \
                    ann_tr.n_extant_terminal_nodes != a_stop_value:
                return False

            # if conditions are not met, tree is fine
            return True

        #####################
        # STOP: Tree height #
        #####################
        elif self.stop == "age":
            # reconstructed tree is too small
            if ann_tr.n_extant_sampled_terminal_nodes < self.min_rec_taxa:
                return False

            # reconstructed tree is too large
            if ann_tr.n_extant_sampled_terminal_nodes > self.max_rec_taxa:
                return False

            # neither too large nor too small, and also not out of control
            #
            # from origin #
            seed_age: float = -1.0
            if ann_tr.with_origin and isinstance(ann_tr.origin_age, float):
                seed_age = ann_tr.origin_age
            # from root #
            elif not ann_tr.with_origin and isinstance(ann_tr.root_age, float):
                seed_age = ann_tr.root_age

            # tree is OK if either
            #     (i)  it reached the max age, or
            #     (ii) its age is smaller than max age (it died!)
            if abs(a_stop_value - seed_age) <= self.epsilon or \
                seed_age < a_stop_value:
                return True

            # somehow tree is taller than the maximum specified height;
            # this block has to be below the one above, otherwise trees
            # that are epsilon-taller than the maximum age get incorrectly
            # rejected
            elif seed_age > a_stop_value:
                return False

        return False

    def generate(self) -> ty.List[AnnotatedTree]:
        """
        Generate samples according to distribution.

        The number of samples will be the specified number of samples
        times the number of replicates (per sample).

        Returns:
            (AnnotatedTree): Valid simulated trees annotated with
                discrete traits
        """

        # output
        output: ty.List[AnnotatedTree] = list()

        # do something with self.runtime_limit
        start_time = time.time()
        ith_sim = 0
        n_failed_attempts = 0
        while len(output) < (self.n_sim * self.n_repl):
            ith_repl = 0

            while ith_repl < self.n_repl:
                # abort due to runtime limit
                ellapsed_time = \
                    pjh.get_ellapsed_time_in_seconds(start_time, time.time())
                if ellapsed_time >= self.runtime_limit:
                    raise ec.RunTimeLimit(self.runtime_limit)

                start_states_idx = (ith_sim * self.n_repl) + ith_repl

                # simulate!
                tr = self.simulate(
                    self.start_states[start_states_idx],
                    self.stop_val[ith_sim],
                    sample_idx=ith_sim)

                # check if tr has right specs
                if self._is_tr_ok(tr, self.stop_val[ith_sim]):
                    output.append(tr)
                    ith_repl += 1

                # tree not good, stay in while loop
                else:
                    # abort due to max number of attempts
                    n_failed_attempts += 1
                    if n_failed_attempts > self.max_n_failed_attempts:
                        raise ec.MaxNFailedAttemptsLimit(n_failed_attempts)

                    continue

            ith_sim += 1

        return output

    #########################
    # Output health methods #
    #########################

    def get_rev_inference_spec_info(self) -> ty.List[str]:
        rev_str_list: ty.List[str] = list()

        return rev_str_list