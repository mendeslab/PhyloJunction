import typing as ty
import numpy as np
import enum
import random
import math
from abc import ABC, abstractmethod
from copy import deepcopy

# pj imports
import phylojunction.utility.exception_classes as ec
import phylojunction.utility.helper_functions as pjh

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


# MacroevolEvent = enum.Enum("Macroevol. event",
#                            ["W_SPECIATION",
#                             "BW_SPECIATION",
#                             "ASYM_SPECIATION",
#                             "EXTINCTION",
#                             "ANAGENETIC_TRANSITION",
#                             "ANCESTOR_SAMPLING"],
#                             start=0)

class MacroevolEvent(enum.Enum):
    W_SPECIATION = 0
    BW_SPECIATION = 1
    ASYM_SPECIATION = 2
    EXTINCTION = 3
    ANAGENETIC_TRANSITION = 4
    ANCESTOR_SAMPLING = 5


class DiscreteStateDependentParameterType(enum.Enum):
    UNDEFINED = 0
    RATE = 1
    PROBABILITY = 2


class DiscreteStateDependentParameter():
    """
    Discrete state-dependent parameter.

    This is the parent class for either state-dependent rates or
    state-dependent parameters.

    Supports vectorization of values only.

    Parameters:
        value (Union[int, float, str]): Either an integer, float or
            string, or a list of any of those.
        name (str): Name of SSE parameter.
        state (int): Integer representing the state this SSE parameter
            is associated to.
        epoch_idx (int): Time slice (epoch) this SSE parameter is
            associated to.
    """

    value: ty.Union[int,
                    float,
                    str,
                    ty.List[ty.Union[int, float, str]]]
    name: str
    state: int
    epoch_idx: int

    def __init__(self,
                 val: ty.Union[int,
                               float,
                               str,
                               ty.List[ty.Union[int, float, str]]],
                 name: str = "",
                 state: int = 0,
                 epoch_idx: int = 1):

        self.name = name
        self.state = state
        self.epoch_idx = epoch_idx
        self.value = []

        # scalar
        if isinstance(val, (float, int, str)):
            # vectorizing if input wasn't in list or ndarray form
            self.value = [float(val)]

        # list
        elif isinstance(val, list):
            if isinstance(val[0], (int, float, str, np.float64)):
                self.value = [float(v) for v in val]

            else:
                raise ec.StateDependentParameterMisspec(
                    message=("Could not recognize type of value argument, "
                             " for parameter " + self.name + ". Cannot "
                             "initialize " + name + "."))

        else:
            raise ec.StateDependentParameterMisspec(
                message=("Argument to value parameter is not either in "
                         "scalar or vectorized form (it is likely an "
                         "object). Cannot initialize " + name + "."))

    @abstractmethod
    def _initialize_str_representation(self) -> None:
        pass

    def __len__(self) -> int:
        if isinstance(self.value, list):
            return len(self.value)

        else:
            return 1


class DiscreteStateDependentRate(DiscreteStateDependentParameter):
    """
    State-dependent rate main class.

    This class is derived from DiscreteStateDependentParameter. It
    holds a value for the rate, the type of rate, the name of
    variable holding it, the state associated to the rate, and
    the time slice (epoch) that rate applies to.

    Supports vectorization of values only.

    Parameters:
        state_tuple (int): Tuple storing all states associated to SSE
            parameter.
        departing_state (int): State event is departing from.
        arriving_state (int): State event is arriving at.
        event (MacroevolEvent): The type of macroevolutionary event
            (e.g., within-region speciation, extinction, etc.) this SSE
            parameter is associated to.
        str_representation (str): The string that is put together
            upon initialization and printed when __str__() is called.
    """

    departing_state: int
    arriving_state: int
    state_tuple: ty.Tuple[int]
    event: MacroevolEvent
    str_representation: str

    def __init__(
            self,
            val: ty.Union[int,
                          float,
                          str,
                          ty.List[ty.Union[int, float, str]]],
            event: MacroevolEvent,
            name: str = "",
            states: ty.List[int] = [],
            epoch_idx: int = 1) -> None:

        # side-effect: initializes
        # self.value, self.name, self.state, self.epoch
        self.states: ty.List[int] = []
        if len(states) == 0:
            self.states = [0, 0, 0]

        else:
            self.states = states

        state = int(self.states[0])
        super().__init__(val=val,
                         name=name,
                         state=state,
                         epoch_idx=epoch_idx)

        self.state_tuple = tuple(int(s) for s in self.states)
        self.departing_state = self.state_tuple[0]
        # does not include time_slice_index 1
        self.arriving_states = self.state_tuple[1:]

        # making sure inputs are ok
        self.event = event
        if self.event == MacroevolEvent.EXTINCTION and \
                self.arriving_states and (len(states) > 0):

            raise ec.StateDependentParameterMisspec(
                message=("Extinction only takes a departing state "
                         "(and no arriving states)."))

        if self.event == MacroevolEvent.W_SPECIATION and \
            len(set(self.arriving_states)) != 1 and \
                self.departing_state not in self.arriving_states:

            raise ec.StateDependentParameterMisspec(
                message=("Within-state speciation requires the "
                         "departing and arriving states to be all the "
                         "same."))

        if self.event == MacroevolEvent.BW_SPECIATION and \
            len(set(self.arriving_states)) != 2 and \
                self.departing_state in self.arriving_states:

            raise ec.StateDependentParameterMisspec(
                message=("Between-state speciation requires the departing "
                         "and the two arriving states to be all different."))

        if self.event == MacroevolEvent.ASYM_SPECIATION and \
            len(set(self.arriving_states)) != 2 and \
                self.departing_state not in self.arriving_states:

            raise ec.StateDependentParameterMisspec(
                message=("Asymmetric-state speciation requires the "
                         "departing and one of the arriving states to be "
                         "the same, and the other arriving state to be "
                         "different. Exiting..."))

        if self.event == MacroevolEvent.ANAGENETIC_TRANSITION and \
            len(set(self.arriving_states)) != 1 and \
                self.departing_state in self.arriving_states:

            raise ec.StateDependentParameterMisspec(
                message=("State transition requires the departing and "
                         "arriving states to be different."))

        if self.event == MacroevolEvent.ANCESTOR_SAMPLING \
                and self.arriving_states:

            raise ec.StateDependentParameterMisspec(
                message="Ancestor sampling "
                + "only takes a departing state (and no arriving "
                + "states).")

        self._initialize_str_representation()

    def _initialize_str_representation(self):
        self.str_representation = "Discrete state-dependent rate\n" \
            + "   Name:                " + self.name + "\n" \
            + "   Value:               " \
            + ", ".join(str(v) for v in self.value) + "\n" \
            + "   Departing state:     " + str(self.departing_state) + "\n" \
            + "   Arriving state(s):   " \
            + ", ".join(str(v) for v in self.arriving_states) + "\n" \
            + "   Associated event:    " + str(self.event) + "\n" \
            + "   Epoch:               " + str(self.epoch_idx) + "\n\n"

    def __str__(self) -> str:
        return self.str_representation

    def __repr__(self) -> str:
        return self.str_representation

    def __len__(self) -> int:
        return super().__len__()


class DiscreteStateDependentProbability(DiscreteStateDependentParameter):
    """
    State-dependent probability main class.

    This class is derived from DiscreteStateDependentParameter. It
    holds a value for the probability, the name of the variable
    holding it, the state associated to the probability, and
    the time slice (epoch) that probability applies to.

    Supports vectorization of values only.

    Parameters:
        str_representation (str): The string that is put together
            upon initialization and printed when __str__() is called.
    """

    state_representation: str

    def __init__(self,
                 val: ty.Union[int,
                               float,
                               str,
                               ty.List[ty.Union[int, float, str]]],
                 name: str = "",
                 state: int = 0,
                 epoch_idx: int = 1) -> None:

        # side-effect: populates self.values (list of floats)
        super().__init__(val=val,
                         name=name,
                         state=state,
                         epoch_idx=epoch_idx)

        for v in self.value:
            if v < 0.0:
                raise ec.NotBetweenZeroAndOneError(name, "negative")

            elif v > 1.0:
                raise ec.NotBetweenZeroAndOneError(name, "positive")

        self._initialize_str_representation()

    def _initialize_str_representation(self) -> None:
        self.str_representation = "Discrete state-dependent probability\n" \
            + "   Name:   " + self.name + "\n" \
            + "   Value:  " \
            + ", ".join(str(v) for v in self.value) + "\n" \
            + "   State:  " + str(self.state) + "\n" \
            + "   Epoch:  " + str(self.epoch_idx) + "\n\n"

    def __str__(self) -> str:
        return self.str_representation

    def __repr__(self) -> str:
        return self.str_representation

    def __len__(self) -> int:
        return super().__len__()


class DiscreteStateDependentParameterManager:
    """Class for checking, organizing, and storing SSE parameters.

    Upon initialization, this class verifies that user provided SSE
    parameters (and their time slice annotation) is acceptable. Then
    given a seed age for anchoring ages and converting them into times,
    and time slice age ends, this class produces computes time slice
    time ends.

    This class does not care about vectorization. It manipulates
    instances of DiscreteStateDependentParameter, which then in turn
    contain multiple values if vectors have been passed by user.

    Currently there is no support for vectorization of seed ages nor
    of time slices.

    Parameters:
        matrix_state_dep_params (DiscreteStateDependentParameter): 2D
            list, with first dimension being time slices, the second
            being SSE parameters. This is just a raw container of SSE
            parameters, and is not used for computations.
        state_dep_params_dict (dict): Dictionary with state indices
            (int) as keys, and a list of SSE parameters as values.
            This is the container that is used by MacroEvolEventHandler
            when calculating rates.
        seed_age (float, optional): Age of origin or root used to
            anchor time slice age ends so as to convert them into time
            slice time ends.
        slice_age_ends (float): List of time slice age ends.
        slice_t_ends (float, optional): List of time slice time ends.
            It is automatically filled upon initialization.
        n_slices_to_ignore (int): How many time slices are in their
            totality older than the whole tree (the user may have
            specified more time slice age ends than necessary given
            an origin or root time for the tree). Defaults to 0.
        param_type (DiscreteStateDependentParameterType): Attribute
            holding what type of SSE parameter this is (rate or
            probability).
        epsilon (float): Threshold for considering a difference equal
            to zero.
    """

    matrix_state_dep_params: \
        ty.List[ty.List[DiscreteStateDependentParameter]]
    state_dep_params_dict: \
        ty.Dict[int,
                ty.List[ty.List[DiscreteStateDependentParameter]]]
    seed_age: ty.Optional[float]
    slice_age_ends: ty.List[float]
    slice_t_ends: ty.Optional[ty.List[float]]
    n_slices_to_ignore: int
    param_type: DiscreteStateDependentParameterType
    epsilon: float

    # NOTE: This class is flexible in that it allows different parameter
    # numbers per time slice, for whatever that is worth. However, the
    # user interface has a check inside make_SSEStash()
    # that forces the user to specify the same number of parameters in
    # all time slices

    def __init__(self,
                 matrix_state_dep_params: \
                    ty.List[ty.List[DiscreteStateDependentParameter]],
                 total_state_count: int,
                 seed_age_for_time_slicing: ty.Optional[float] = None,
                 list_time_slice_age_ends: ty.Optional[ty.List[float]] = None,
                 epsilon: float = 1e-12) -> None:
        
        self.epsilon = epsilon

        self.n_slices_to_ignore = 0

        # 1D: time slices
        # 2D: list of SSE params
        self.matrix_state_dep_params = matrix_state_dep_params

        # make sure it's all rates or all probabilities and
        # initialize self.param_type
        self._init_check_single_and_init_param_type()

        self.state_count = total_state_count

        # this is the origin or root age, and it is used to anchor
        # the user-specified ages to convert it to time
        self.seed_age = seed_age_for_time_slicing

        # default is one slice, but this gets updated below
        self.n_time_slices = 1

        # age ends (larger in the past, 0.0 in the present)
        self.slice_age_ends = [0.0]

        # default slice t's at present
        # (t = seed_age = stop_value with "age" condition)
        self.slice_t_ends = [self.seed_age]

        # if user provides age ends for time slices
        if list_time_slice_age_ends:
            # updates self.n_slices_to_ignore
            #         self.slice_age_ends
            #         self.slice_t_ends
            self._init_update_slice_age_and_t_ends(
                list_time_slice_age_ends)

        # if rate, must have at least one per time slice;
        # if probability, muast have one per state for all time slices
        self._init_check_correct_number_params_per_time_slice()

        # repeated SSE parameters per time slice are not allowed
        self._init_check_repeated_parameter_in_time_slice()

        # side effect: initializes and populates self.state_dep_params_dict
        self._init_matrix_state_dep_params_dict()

    ######################
    # populate functions #
    ######################
    def _init_matrix_state_dep_params_dict(self) -> None:
        """Populate class member dictionary with SSE parameters.
        
        Raises:
            ObjInitInvalidArgError: Is raised if the state (int)
                associated to the SSE parameter is not an integer
                between 0 and the total number of states.
        """

        # key (int): state
        # value (2D list): 1D are time slices, 2D are SSE parameters
        self.state_dep_params_dict: \
            ty.Dict[int,
                    ty.List[ty.List[DiscreteStateDependentParameter]]] = \
                        dict((s, [[] for j in range(self.n_time_slices)])
                             for s in range(self.state_count))

        # t-th time slice
        for t, list_params in enumerate(
                self.matrix_state_dep_params[
                self.n_slices_to_ignore:(self.n_time_slices + self.n_slices_to_ignore)
                ]):
            for param in list_params:
                try:
                    self.state_dep_params_dict[param.state][t].append(param)

                except KeyError:
                    arg_name: str = "flat_rate_mat"
                    if self.param_type == DiscreteStateDependentParameterType\
                            .PROBABILITY:
                        arg_name = "flat_rate_prob"

                    raise(ec.ObjInitInvalidArgError(
                        ("'sse_stash' ('DiscreteStateDependentParameter'"
                         "Manager' in the backend)"),
                        arg_name + (" ('matrix_state_dep_params' in the "
                         "backend)"),
                        message="Parameter " + param.name + "'s associated " \
                            + "(possibly departing) state is " \
                            + str(param.state) \
                            + ", but this state is not represented between " \
                            + "0 and " + str(self.state_count) \
                            + (". Likely the total number of states was "
                               "misspecified.")))

    def _init_update_slice_age_and_t_ends(self,
                                          list_time_slice_age_ends) -> None:
        """Update class members related to time slices.

        This method will update:
            (i)   self.n_time_slices
            (ii)  self.slice_age_ends
            (iii) self.slice_t_ends
        """

        # we do not want to change this object outside
        # of this class instance, because other state-dep
        # param managers may be created (and will be taking
        # 'list_time_slice_age_ends' as arguments)
        self.slice_age_ends = deepcopy(list_time_slice_age_ends)

        # appends age end of time slice ending in present
        self.slice_age_ends.append(0.0)

        # convert into time ends
        # (0.0 at origin, seed_age = stop_val for "age" stop
        # condition in the present)
        if self.seed_age:
            # reset slice_t_ends member
            self.slice_t_ends: ty.List[float] = list()

            # old first, young last
            n_slices_to_ignore = 0
            for age_end in self.slice_age_ends:
                # we ignore user-specified age ends that for some reason
                # older than the seed (origin/root) age; must have been
                # a user oversight...
                if age_end > self.seed_age:
                    n_slices_to_ignore += 1
                    continue

                else:
                    self.slice_t_ends.append(self.seed_age - age_end)

            # these slices are older than the tree, so we will ignore them!
            self.n_slices_to_ignore = n_slices_to_ignore

            # reflects n_slices_to_ignore!
            self.n_time_slices = len(self.slice_t_ends)

    ##########################
    # health check functions #
    ##########################
    def _init_check_single_and_init_param_type(self) -> None:
        """Check discrete SSE paramaters are all of the same kind.
        
        This initialization method verifies that the manager was passed
        either all SSE rates, or SSE probabilities. It will raise an
        error if not.

        The method has a side-effect, the initialization of
        self.param_type to either:
            (i)  DiscreteStateDependentParameterType.RATE
            (ii) DiscreteStateDependentParameterType.PROBABILITY.

        Raises:
            ObjInitRequireSameParameterTypeError: Is raised if
                parameters passed to manager by the user are not
                all of the same type.
        """

        different_par_type_set: \
            ty.Set(DiscreteStateDependentParameter) = set()

        param_type: DiscreteStateDependentParameterType = None

        # iterating over time slices (epochs)
        for list_params in self.matrix_state_dep_params:
            for param in list_params:
                this_param_type = type(param)
                param_type = this_param_type

                if this_param_type not in different_par_type_set:
                    different_par_type_set.add(this_param_type)

        if len(different_par_type_set) >= 2:
            raise ec.ObjInitRequireSameParameterTypeError(
                "DiscreteStateDependentParameterManager",
                len(different_par_type_set))
        
        if param_type is DiscreteStateDependentRate:
            self.param_type = \
                DiscreteStateDependentParameterType.RATE

        elif param_type is DiscreteStateDependentProbability:
            self.param_type = \
                DiscreteStateDependentParameterType.PROBABILITY

    def _init_check_correct_number_params_per_time_slice(self) -> None:
        """Check if number of SSE paramaters per time slice is valid.

        This initialization method verifies that the manager was passed
        valid numbers of discrete SSE parameters in all time slices. If
        parameters are probabilities, then we need one per state per
        time slice. If rates, we need at least one non-zero rate per
        time slice.

        Raises:
            ObjInitIncorrectDimensionError: Is raised if at least one
                state-dependent probability is not specified for each
                state for each time slice. Is also raised if not at
                least one state-dependent rate is not specified for
                each time slice.
            ObjInitRequireNonZeroStateDependentParameterError: Is
                raised if at least one non-zero state-dependent rate
                is not provided for all time slices.
        """
        
        # one SSE probability per state per time slice (at the moment SSE
        # probabilities cannot be annotated with a time slice they belong)
        if self.param_type == \
                DiscreteStateDependentParameterType.PROBABILITY:
            # time slice t
            for t, sse_param_list in enumerate(self.matrix_state_dep_params):
                n_sse_params: int = len(sse_param_list)
                if n_sse_params != self.state_count:
                    raise ec.ObjInitIncorrectDimensionError(
                        ("'sse_stash' ('DiscreteStateDependentParameter"
                         "Manager' in the backend)"),
                        ("'flat_prob_mat' ('matrix_state_dep_params[" + str(t) \
                         + "]' in the backend)"),
                        n_sse_params,
                        exp_len=self.state_count
                    )
            
            # min_param_count = \
            #     min(len(t) for t in self.matrix_state_dep_params)
            
            # max_param_count = \
            #     max(len(t) for t in self.matrix_state_dep_params)
            
            # if min_param_count != max_param_count:
            #     raise ec.ObjInitIncorrectDimensionError(
            #         "DiscreteStateDependentParameterManager",
            #         "matrix_state_dep_params",
            #         max_param_count,
            #         exp_len=min_param_count
            #     )
            
            # self.n_time_slices = len(self.matrix_state_dep_params)

        # at least one non-zero SSE rate per time slice
        elif self.param_type == \
                DiscreteStateDependentParameterType.RATE:
            # t-th time slice
            for t in range(self.n_time_slices):
                tth_param_list = self.matrix_state_dep_params[t]

                # first we check that at least one rate was provided
                if len(tth_param_list) < 1:
                    raise ec.ObjInitIncorrectDimensionError(
                        ("'sse_stash' ('DiscreteStateDependentParameter"
                         "Manager' in the backend)"),
                        "'flat_rate_mat' ('matrix_state_dep_params[" + str(t+1) \
                        + "]' in the backend)",
                        0,
                        1,
                        at_least=True)
                
                # second, we check that of the provided rates, at 
                # least one is non-zero
                else:
                    # each param.value will contain values for all samples
                    # of a specific type of rate, so if ANY of the values of
                    # param_vals_list is zero, it means that specific rate is 0.0
                    all_rate_types_are_zero = True
                    for param in tth_param_list:
                        if any(v != 0.0 for v in param.value):
                            all_rate_types_are_zero = False

                    if all_rate_types_are_zero:
                        raise ec.ObjInitRequireNonZeroStateDependentParameterError(
                                "DiscreteStateDependentParameterManager",
                                t)

    def _init_check_repeated_parameter_in_time_slice(self) -> None:
        """Make sure parameters are not repeated in epoch.
        
        Raises:
            ObjInitRepeatedStateDependentParameterError: Is raised if
                there is more than one SSE rate (or probability)
                representing the same state.
        """

        states_per_epoch_dict = \
            dict((t, set()) for t in range(self.n_time_slices))

        # t-th time slice (note that we slice matrix_state_dep_params
        # because we want to only use the time slices that are younger
        # than the tree seed's age)
        for t, list_params in enumerate(
                self.matrix_state_dep_params[:self.n_time_slices]):
            for param in list_params:
                st_event_tup = tuple()
                sts = tuple()

                if self.param_type is \
                        DiscreteStateDependentParameterType.RATE:
                    sts = param.state_tuple
                    st_event_tup = tuple([sts, param.event])

                if self.param_type is \
                        DiscreteStateDependentParameterType.PROBABILITY:
                    sts = param.state
                    st_event_tup = tuple([sts])                    

                # if this is raising when you would not expect it to,
                # double check that the rates are being passed
                # by epoch first, then by type of rate, e.g.,
                # birth0, death0, birth1, death1, etc.
                if st_event_tup not in states_per_epoch_dict[t]:
                    states_per_epoch_dict[t].add(st_event_tup)
                
                else:
                    raise ec.ObjInitRepeatedStateDependentParameterError(
                        t, st_event_tup
                    )
                
                states_per_epoch_dict[t].add(sts)

    ###########
    # getters #
    ###########
    def state_dep_params_at_time(
            self,
            a_time: float,
            params_matrix: ty.Optional[
                ty.List[ty.List[DiscreteStateDependentParameter]]] = None) \
            -> ty.List[DiscreteStateDependentParameter]:
        """Return list of SSE parameters given a specific time.

        This method finds the time slice a specific time belongs to,
        and retrieves the list of SSE parameters pertaining to that
        time slice. This list is retrieved from a user-specified
        matrix, or from the matrix member of the enclosing class,
        self.matrix_state_dep_params.

        Args:
            a_time (float): A time at which we want to grab SSE
                parameters.
            param_matrix (DiscreteStateDependentParameter): A 2D list
                where the first dimension are time slices, and the
                second dimension are SSE parameters.
        
        Returns:
            (DiscreteStateDependentParameter): a list of SSE parameters
                from the appropriate time slice.
        """

        # find out which time slice (by getting its index) to get
        # SSE parameter value from
        time_slice_index: int = 0
        # (if n_time_slices is 0, for loop has no effect)

        # number of time slices is the number of age ends provided by
        # user, but only those that are as old or younger than the tree!!!
        for time_slice_index in range(0, self.n_time_slices):
            # try:
            if isinstance(self.slice_t_ends, list) and \
                    isinstance(self.slice_t_ends[time_slice_index], float):
                # time end of this time slice
                time_slice_t_end = \
                    ty.cast(float, self.slice_t_ends[time_slice_index])
                
                # if time is larger than this time slice time end OR
                # if time is the same as the time slice, we continue to
                # the next time slice and check again
                if a_time > time_slice_t_end or \
                        (abs(a_time - time_slice_t_end) <= self.epsilon):
                    continue

            break
        
        # now that we know the time slice, we get the SSE parameter
        #
        # params_matrix will have been provided as a matrix conditioned on a
        # departing state
        if params_matrix is not None:
            return params_matrix[time_slice_index]

        # when no params_matrix is specified, just grab
        # all parameters at time t, but warning:
        # no control for the order (in different time slices)
        # the parameters have been passed by user!!!!
        else:
            return self.matrix_state_dep_params[time_slice_index]

    def __len__(self) -> int:
        if isinstance(self.matrix_state_dep_params, list):
            return sum([len(list_atomic_rates_in_time_slice)
                        for list_atomic_rates_in_time_slice
                        in self.matrix_state_dep_params])

        else:
            return 0


# class DiscreteStateDependentParameterManager_old_working:
#     """Stash for discrete state-dependent parameters and time slices

#     At the moment, this class does not care about vectorization. It manipulates
#     whole instances of DiscreteStateDependentParameter, which then in turn
#     contain multiple values if vectors have been passed by user.

#     Later, vectorization of seed ages for time slicing and time slices themselves
#     might become vectorized.
#     """

#     seed_age: ty.Optional[float]
#     slice_age_ends: ty.List[float]
#     slice_t_ends: ty.List[ty.Optional[float]]

#     # NOTE: This class is flexible in that it allows different parameter
#     # numbers per time slice, for whatever that is worth. However, the
#     # user interface has a check inside make_SSEStash()
#     # that forces the user to specify the same number of parameters in
#     # all time slices

#     def __init__(self,
#                  matrix_atomic_rate_params:
#                  ty.List[ty.List[DiscreteStateDependentRate]],
#                  total_state_count: int,
#                  seed_age_for_time_slicing: ty.Optional[float] = None,
#                  list_time_slice_age_ends: ty.Optional[ty.List[float]] = None,
#                  epsilon: float = 1e-12):
#         # how do we need to access parameters and values during simulation?
#         # we probably need two types of parameter containers.
#         # type 1: records "high-level" governing parameters, rho, sigma, phi
#         # type 2: records "low-level" event-specific rates, which are used
#         #         to simulate next event times, classes, and outcomes

#         # 1D: time slices, 2D: list of atomic rate params
#         self.atomic_rate_params_matrix = matrix_atomic_rate_params
#         self.state_count = total_state_count
#         # TODO: make this vectorized next, so a prior can be put on seed age

#         # this is the origin or root age, and it is used to anchor the
#         # user-specified ages to convert it to time
#         self.seed_age = seed_age_for_time_slicing
#         self.n_time_slices = 1  # default is one slice
#         # default slice ends at present (age = 0.0)
#         self.slice_age_ends = [0.0]
#         # default slice t's at present (t = seed_age = stop_value
#         # with "age" condition)
#         self.slice_t_ends = [self.seed_age]

#         # age ends (larger in the past, 0.0 in the present)
#         if list_time_slice_age_ends:
#             self.n_time_slices += len(list_time_slice_age_ends)
#             self.slice_age_ends = list_time_slice_age_ends
#             self.slice_age_ends.append(0.0)
#             # appends age end of time slice
#             # ending in present

#             # convert into time ends (0.0 at origin, seed_age = stop_val for
#             # "age" stop condition in the present)
#             if self.seed_age:
#                 self.slice_t_ends = \
#                     [self.seed_age - age_end
#                      for age_end
#                      in self.slice_age_ends]
#                 # no need to append seed_age, because self.slice_age_ends
#                 # already has 0.0 in it

#         self.atomic_rate_params_dict: \
#             ty.Dict[int,
#                     ty.List[ty.List[DiscreteStateDependentRate]]] = \
#             dict((s, [[] for j in range(self.n_time_slices)])
#                  for s in range(self.state_count))

#         # side effect: initializes self.atomic_rate_params_dict
#         self.init_atomic_rate_param_dict(matrix_atomic_rate_params)
#         # self.atomic_rate_params_dict =
#         # { state0:
#         #           [ [ #slice1#; atomic_param1, atomic_param2, ...] [ #slice2#; atomic_param1, atomic_param2 ] ], ... ]
#         #   state1:
#         #           [ [ #slice1#; atomic_param1, atomic_param2, ...] [ #slice2#; atomic_param1, atomic_param2 ] ], ... ]
#         # }

#         self.epsilon = epsilon

#     def __len__(self) -> int:
#         if isinstance(self.atomic_rate_params_matrix, list):
#             return sum(
#                 [len(list_atomic_rates_in_time_slice)
#                  for list_atomic_rates_in_time_slice
#                  in self.atomic_rate_params_matrix])

#         else:
#             return 0

#     def init_atomic_rate_param_dict(
#             self,
#             matrix_atomic_rate_params:
#             ty.List[ty.List[DiscreteStateDependentRate]]):
#         # original implementation (2nd dimension were states, rather than all rates from all states together)
#         # for k, s_state_list in enumerate(matrix_atomic_rate_params):

#         # k-th time slice
#         for k, list_atomic_rate_params in enumerate(matrix_atomic_rate_params):

#             # original implementation had this extra loop
#             # for s, list_atomic_rate_params in enumerate(s_state_list):

#             for atomic_rate_param in list_atomic_rate_params:
#                 try:
#                     self.atomic_rate_params_dict[atomic_rate_param.departing_state][k].append(atomic_rate_param)

#                 except Exception as e:
#                     exit("Parameter " + atomic_rate_param.name
#                          + "'s associated departing state is "
#                          + str(atomic_rate_param.departing_state)
#                          + ", but this departing state is not represented "
#                          + "between 0 and " + str(self.state_count)
#                          + ". Likely the total number of states was "
#                          + "misspecified.")

#     def atomic_rate_params_at_time(self, atomic_rate_params_matrix, a_time: float):
#         # see where a_time falls within (get time_slice_index)
#         # grab list of atomic rate params at that time_slice_index
#         # print("self.n_time_slices = " + str(self.n_time_slices))

#         time_slice_index = -1
#         # try:
#         while time_slice_index < self.n_time_slices:
#             time_slice_index += 1

#             try:
#                 if isinstance(self.slice_t_ends, list) and \
#                         isinstance(self.slice_t_ends[time_slice_index],
#                                    float):
#                     time_slice_t_end = \
#                         ty.cast(float, self.slice_t_ends[time_slice_index])

#                     if a_time > time_slice_t_end or \
#                             (abs(a_time - time_slice_t_end) <= self.epsilon):
#                         continue

#             # self.slice_t_ends will be None if no seed age or time slice
#             # age ends are provided
#             except Exception as e:
#                 print("Exception 3 inside discrete_sse.py: ",
#                       type(e).__name__, " - ", e)
#                 # in which case we want the index to just be 0
#                 time_slice_index = 0

#             break

#         # adjusting
#         # if time_slice_index > (self.n_time_slices - 1):
#         #     time_slice_index = self.n_time_slices - 1

#         return atomic_rate_params_matrix[time_slice_index]


class MacroevolEventHandler():
    """Class for handling state-dependent rates.

    Some of the methods in this class are vector-aware (through
    parameter value_idx), because instances of this class will
    be directly called by dn_sse.simulate() method, which in turn
    will try to access parameter values stored in a
    #-simulations-sized list.

    Parameters:
        sse_rate_manager (DiscreteStateDependentParameterManager):
            Object that checks, organizes and stores SSE rates.
        state_count (int): Number of states characterizing the SSE
            process.
        n_time_slices (int): Number of time slices.
        slice_age_ends (float): List of time slice age ends.
        slice_t_ends (float, optional): List of time slice time ends.
        seed_age (float, optional): Age of seed (origin or root),
            stored in 'state_dep_rate_manager'.
        str_representation (str): The string that is put together
            upon initialization and printed when __str__() is called.
    """

    sse_rate_manager: DiscreteStateDependentParameterManager
    state_count: int
    n_time_slices: int
    slice_age_ends: ty.List[float]
    slice_t_ends: ty.Optional[ty.List[float]]
    seed_age: ty.Optional[float]
    str_representation: str

    def __init__(self,
                 sse_rate_manager: \
                    DiscreteStateDependentParameterManager) -> None:

        self.sse_rate_manager = sse_rate_manager
        self.state_count = self.sse_rate_manager.state_count
        self.n_time_slices = self.sse_rate_manager.n_time_slices
        self.slice_age_ends = self.sse_rate_manager.slice_age_ends
        self.slice_t_ends = self.sse_rate_manager.slice_t_ends
        # if there are more time slice age ends than there are time ends,
        # for each extra age end we add a 'Older than tree' string to this
        # list to use it when printing an instance of this class as a string
        filler_t_ends_list = \
            ["Older than tree" for i in range(
             len(self.slice_age_ends) - len(self.slice_t_ends))]
        # and here we prepend those filler strings, so they can be printed
        # as the time end of time slices that are older than the seed age
        # of the tree
        self.slice_t_ends_for_printing = \
            filler_t_ends_list + self.slice_t_ends

        self.seed_age = self.sse_rate_manager.seed_age

        # side-effect: populates self.str_representation
        self._initialize_str_representation()

    # init methods
    def _initialize_str_representation(self) -> None:
        """Populate str_representation class member."""

        self.str_representation = "State-dependent rates"

        # state s
        for s, sse_rate_state_mat in \
                self.sse_rate_manager.state_dep_params_dict.items():
            self.str_representation += "\n  State " + str(s) + ":\n"

            # time slice t
            for t, list_sse_rate_slice in enumerate(sse_rate_state_mat):
                if self.seed_age and isinstance(self.slice_t_ends_for_printing, list):
                    time_slice_t_end = self.slice_t_ends_for_printing[t]

                    self.str_representation += \
                        "    Time slice " + str(t + 1) + " (time = " \
                        + (str(round(time_slice_t_end, 4)) \
                           if isinstance(time_slice_t_end, float) \
                           else time_slice_t_end) + ", age = " \
                        + str(round(self.slice_age_ends[t], 4)) + ")\n"

                else:
                    self.str_representation += "    Time slice " \
                        + str(t + 1) + " (age = " \
                        + str(round(self.slice_age_ends[t], 4)) + ")\n"

                for sse_rate in list_sse_rate_slice:
                    self.str_representation += "      " + sse_rate.name \
                        + " = " + ", ".join(str(v) for v in sse_rate.value) \
                        + "\n"
                    
    # this function deals with vectorization
    def total_rate(
            self,
            a_time: float,
            state_representation_dict: ty.Dict[int, ty.Set[str]],
            value_idx: int = 0,
            departing_state: ty.Optional[int] = None,
            debug: ty.Optional[bool] = False) -> \
                ty.Union[float, ty.Tuple[float, ty.List[float]]]:
        """Get total rate given possible events at time slice and tree.

        This method is called by DnSSE (dn_discrete_sse.py) when it is
        time to draw the time to the next event, which is exponentially
        distributed with rate equal to the total sum of rates being
        computed by this method (the first return in the tuple).

        For the first return in the tuple, this method gets all event
        rates for a given time slice, and adds them up. Different rates
        are weighed depending on their departing states, with the weight
        being the number of lineages representing that state.
        
        For the second return in the tuple, this method adds up rates
        by their departing states, but there is no weighing. The
        grouped rates are put into a list, at index = state index.
        
        The 'departing_state' argument remains as legacy; this method
        never actually gets called with departing_state. This also means
        that a pure float return also never happens (only the tuple
        return).

        Args:
            a_time (float): Forward time (not age!) we want to recover
                rates with.
            state_representation_dict (dict): Dictionary with (int) states
                as keys and list of (str) node labels at that state.
            value_idx (int, optional): Index specifying which parameter
                value (within a vector of #-of-simulations-size) we care about.
                Defaults to 0 (first element).
            departing_state (int, optional): State we can condition rates on.
                Defaults to None.
            debug (bool, optional): Flag for printing debugging messages.
                Defaults to False.

        Returns:
            (float): If no departing_state is provided, returns the global
                rate. Otherwise, returns a state-conditioned total rate and
                a list with each state total rate.
        """

        if debug:
            if departing_state is None:
                print("\nCalculating rate for exponential:")

            else:
                print(("\nCalculating denominator for randomly choosing an "
                       "event to take place:"))

        # TODO: the total_rate (for each time slice) when there is a
        # departing_state can calculated just once and cached; this should
        # be done in sse_rate_manager (need to implement a method in that
        # class's _init_ for calculation, and then implement a getter)
        #
        # we only need to calculate total_rate again (as below) when not
        # conditioning on a departing_state (because of weighing, which
        # depends on the stochastically growing tree)

        # returns
        total_rate: float = 0.0
        state_rates: ty.List[float] = \
            [0.0 for i in range(len(state_representation_dict))]

        for state_idx, nd_labels in state_representation_dict.items():
            # if departing state is provided, then we do not care about rates
            # not departing from it
            if departing_state is not None and departing_state != state_idx:
                continue

            n_lineages_in_state = len(nd_labels)

            if debug:
                if departing_state is None:
                    print("  state " + str(state_idx) + " is represented by " \
                          + str(n_lineages_in_state) + " lineages (weight).")

            # conditioning on state (variable scoped to function!)
            sse_rate_params_from_state_mat = \
                self.sse_rate_manager.state_dep_params_dict[state_idx]

            sse_rate_params_from_state_at_time = \
                self.sse_rate_manager.state_dep_params_at_time(
                    a_time,
                    params_matrix=sse_rate_params_from_state_mat)

            for sse_rate_param_from_state in \
                    sse_rate_params_from_state_at_time:
                w: float = 1.0  # weight

                # we are getting the total rate for the whole tree,
                # before picking a node; the weight is then the number of
                # targetable nodes at a given state
                if departing_state is None:
                    w = float(n_lineages_in_state)

                if debug:
                    print("    This rate's value is "
                          + str(sse_rate_param_from_state.value[value_idx]))

                # value_idx is what introduces vectorization here
                # if we have 100 lambda values (b/c we are doing 100 simulations)
                # then we will be only looking at the i-th (0 < i < 99) lambda at a time
                state_rates[state_idx] += sse_rate_param_from_state.value[value_idx]
                total_rate += w * sse_rate_param_from_state.value[value_idx]

        if debug:
            print("Total rate =", total_rate)
        
        # total_rate will be simply the sum of all possible event rates,
        # normalized by the number of lineages representing each and all
        # possible departing states
        #
        # state_rates just groups all rates by their departing states,
        # no weighing
        #
        # this is the only return we actually use in dn_discrete_sse (DnSSE)
        # when determining the rate of (any) events for the Poisson process
        # (the rate of the exponential)
        if departing_state is None:
            return total_rate, state_rates

        # total_rate will be the sum of event rates that start at the
        # specified departing_state, with all events weighed equally
        #
        # we currently do not use this return
        else:
            return total_rate

    # this function deals with vectorization
    def sample_event_sse_rate_param(
            self,
            denominator: float,
            a_time: float,
            state_indices: ty.List[int],
            value_idx: int = 0,
            a_seed: ty.Optional[float] = None,
            debug: ty.Optional[bool] = False):
        """Return one-sized list with a random SSE rate.

        This method accesses the sse_rate_manager member to grab all
        SSE rates departing from one or more focal states. It then
        weighs each SSE rate by their value divided by the sum of the
        values of all SSE rates departing from the focal state. Given
        these weighed SSE rates, it randomly picks one and returns.

        This method is called by DnSSE (dn_discrete_sse.py) when
        randomly choosing an event, after an event time and a node have
        already been chosen.

        Args:
            denominator (float): Normalization term for computing weights
                of different events.
            a_time (float): When the event will take place.
            state_indices (int): Event is conditioned on departing from
                this state(s).
            value_idx (int, optional): Index specifying which parameter
                value (within a vector of #-of-simulations-size) we care
                about. Defaults to 0 (first element).
            a_seed (int, optional): Random seed for simulation (not
                working at the moment). Defaults to None.
            debug (bool, optional): Flag for printing debugging messages.
                Defaults to False.

        Returns:
            DiscreteStateDependentRate: List with a single SSE rate.
        """

        if debug:
            print("    Will sample event from state_indices = "
                  + ", ".join(str(i) for i in state_indices))

        if a_seed:
            random.seed(a_seed)

        all_states_sse_rate_params_at_time: \
            ty.List[DiscreteStateDependentRate] = list()
        # weights for sampling proportional to rate value
        ws: ty.List[float] = list()

        # produce weight list (one weight per SSE rate)
        # so we can sample SSE rates proportionally
        for state_idx in state_indices:
            this_state_sse_rate_params_mat = \
                self.sse_rate_manager.state_dep_params_dict[state_idx]
            this_state_sse_rate_params_at_time_list = \
                self.sse_rate_manager.state_dep_params_at_time(
                    a_time,
                    params_matrix=this_state_sse_rate_params_mat)
            all_states_sse_rate_params_at_time += \
                this_state_sse_rate_params_at_time_list

            # total rate of outcomes must depend
            # on "adjacent" states across events
            for sse_rate_param in \
                    this_state_sse_rate_params_at_time_list:
                ws.append((sse_rate_param.value[value_idx] / denominator))

        if debug:
            print("    \nThe following events are allowed:")
            print("   " + "   ".join(
                str(ap) for ap in all_states_sse_rate_params_at_time))

        return random.choices(
            [sse_rate_param for sse_rate_param
             in all_states_sse_rate_params_at_time], weights=ws)

    def __len__(self) -> int:
        if self.sse_rate_manager:
            return len(self.sse_rate_manager)

        else:
            return 0

    def __str__(self) -> str:
        return self.str_representation

    def __repr__(self) -> str:
        return self.str_representation

    # def get_gcf(self):
    #     pass

    # def get_length(self):
    #     pass


class DiscreteStateDependentProbabilityHandler():
    """Class for handling state-dependent taxon sampling.

    This is the sister class to MacroevolEventHandler. Note that time
    slices do not have to be the same as those for SSE rates.

    Parameters:
        state_dep_prob_manager (DiscreteStateDependentParameterManager):
            Object that checks, organizes and stores SSE probabilities.
        state_count (int): Number of states characterizing the SSE
            process.
        n_time_slices (int): Number of time slices.
        slice_age_ends (float): List of time slice age ends.
        slice_t_ends (float, optional): List of time slice time ends.
        str_representation (str): The string that is put together
            upon initialization and printed when __str__() is called.
    """

    state_dep_prob_manager: DiscreteStateDependentParameterManager
    state_count: int
    n_time_slices: int
    seed_age: ty.Optional[float]
    slice_age_ends: ty.List[float]
    slice_t_ends: ty.Optional[ty.List[float]]
    str_representation: str

    def __init__(
            self,
            state_dep_prob_manager: DiscreteStateDependentParameterManager) \
            -> None:

        self.state_dep_prob_manager = state_dep_prob_manager
        self.state_count = self.state_dep_prob_manager.state_count
        self.n_time_slices = self.state_dep_prob_manager.n_time_slices
        self.seed_age = self.state_dep_prob_manager.seed_age
        self.slice_t_ends = self.state_dep_prob_manager.slice_t_ends
        self.slice_age_ends = self.state_dep_prob_manager.slice_age_ends
        filler_t_ends_list = \
            ["Older than tree" for i in range(
             len(self.slice_age_ends) - len(self.slice_t_ends))]
        self.slice_t_ends_for_printing = \
            filler_t_ends_list + self.slice_t_ends

        # side-effect: initializes self.str_representation
        self._initialize_str_representation()

    def _initialize_str_representation(self) -> None:
        """Populate str_representation class member."""

        self.str_representation = "State-dependent probabilities"

        # state s
        for s, atomic_rates_state_mat in \
                self.state_dep_prob_manager.state_dep_params_dict.items():
            self.str_representation += "\n  State " + str(s) + ":\n"

            # time slice t
            for t, list_state_dep_params_slice in \
                    enumerate(atomic_rates_state_mat):
                if self.seed_age and \
                        isinstance(self.slice_t_ends_for_printing, list):
                    time_slice_t_end = self.slice_t_ends_for_printing[t]

                    self.str_representation += "    Time slice " \
                        + str(t + 1) + " (time = " \
                            + (str(round(time_slice_t_end, 4)) if \
                               isinstance(time_slice_t_end, float) else \
                                time_slice_t_end) \
                        + ", age = " \
                        + str(round(self.slice_age_ends[t], 4)) + ")\n"

                else:
                    self.str_representation += "    Time slice " \
                        + str(t + 1) + " (age = " \
                        + str(round(self.slice_age_ends[t], 4)) + ")\n"

                for state_dep_param in list_state_dep_params_slice:
                    self.str_representation += "      " \
                        + state_dep_param.name + " = "
                    self.str_representation += ", ".join(
                        str(v) for v in state_dep_param.value) + "\n"

    def _state_dep_prob_at_time(self, a_time, state_idx) \
            -> ty.List[DiscreteStateDependentProbability]:
        """Get the SSE probability of a state at a specific time.
        
        Args:
            a_time (float): Time at which we want to know the
                probability of sampling a taxon at the state
                represented by 'state_idx'.
            state_idx (int): Integer representing a state.

        Returns:
            DiscreteStateDependentProbability: A list with a single
                SSE probability inside.
        """

        # scoped to total_rate
        state_cond_prob_matrix = \
            self.state_dep_prob_manager \
                .state_dep_params_dict[state_idx]  # conditioning


        # in the context of rates (multiple rates per departing state)
        # state_dep_params_at_time returns a list of rates,
        # but here, we have a single probability inside a list
        state_cond_prob_at_time = \
            self.state_dep_prob_manager \
                .state_dep_params_at_time(
                    a_time,
                    params_matrix=state_cond_prob_matrix)

        return state_cond_prob_at_time

    def randomly_decide_taxon_sampling_at_time_at_state(
            self, a_time, state_idx, sample_idx) \
            -> bool:
        """

        Raises:
        """

        prob_at_state_in_slice_list = \
            self._state_dep_prob_at_time(a_time, state_idx)

        if len(prob_at_state_in_slice_list) > 1:
            exit(("Should only have one probability deterministic node "
                  " per state per slice. Exiting..."))

        prob_at_state_in_slice_ith_sample = \
            prob_at_state_in_slice_list[0].value[sample_idx]
        
        if prob_at_state_in_slice_ith_sample == 1.0:
            return True

        elif prob_at_state_in_slice_ith_sample == 0.0:
            return False

        elif np.random.binomial(1, prob_at_state_in_slice_ith_sample):
            return True

        else:
            return False

    def __len__(self) -> int:
        if self.state_dep_prob_manager:
            return len(self.state_dep_prob_manager)

        else:
            return 0

    def __str__(self) -> str:
        return self.str_representation

    def __repl__(self) -> str:
        return self.str_representation


class SSEStash():
    """
    Class for stashing state-dependent rate and probability handlers.

    Parameters:
        meh (MacroevolEventHandler): Instance of state-dependent rate
            handler class.
        prob_handler (DiscreteStateDependentProbabilityHandler):
            Instance of state-dependent probability handler class.
    """

    meh: MacroevolEventHandler
    prob_handler: DiscreteStateDependentProbabilityHandler
    str_representation: str

    def __init__(
            self,
            macroevol_event_handler: MacroevolEventHandler,
            state_dep_prob_handler:
            ty.Optional[DiscreteStateDependentProbabilityHandler] = None) \
            -> None:

        self.meh = macroevol_event_handler
        self.prob_handler = state_dep_prob_handler

        if state_dep_prob_handler is None:
            self._initialize_missing_prob_handler()

        self._initialize_str_representation()

    # side-effect: creates and populates self.prob_handler if user did
    # not provide it (one probability per state, per time slice)
    def _initialize_missing_prob_handler(self) -> None:
        """Populate class member self.prob_handler.
        
        This method is called when the user only specified
        state-dependent rates for an SSE process, but no
        state-dependent probabilities. We still need the probabilities,
        so we create them here from scratch, assuming they are all
        1.0.
        
        The class member self.prob_handler is initialized as a
        side-effect.
        """

        matrix_state_dep_probs: \
            ty.List[ty.List[DiscreteStateDependentProbability]] = list()

        # we need a full size matrix of SSE probs, with one list of
        # SSE probs for each epoch with an age end provided by the
        # user (as opposed to the 'actual' number of epochs, which
        # would ignore epochs older than the root or origin)
        #
        # this matrix is later subset in
        # DiscreteStateDependentParameterManager's
        # _init_matrix_state_dep_params_dict
        for t in range(0, len(self.meh.slice_age_ends)):
            state_dep_probs: \
                ty.List[DiscreteStateDependentProbability] = list()

            for st in range(self.meh.state_count):
                # st = j % n_prob_per_slice
                state_dep_prob = \
                    DiscreteStateDependentProbability(
                        name="rho" + str(st) + "_t" + str(t),
                        val=1.0,
                        state=st)

                # so mypy won't complain
                if isinstance(state_dep_prob,
                              DiscreteStateDependentProbability):
                    state_dep_probs.append(state_dep_prob)

            # appending a list of DiscreteStateDependentProbability
            # 1D: time slices, 2D: state-dep prob list
            matrix_state_dep_probs.append(state_dep_probs)

        # we have to throw away the last element in slice_age_ends
        # because it has already been added a 0.0 for the youngest
        # time slice when 'meh' was initialized; if we don't do this
        # the number of time slices inside 'state_dep_probs_manager'
        # will be 1 too large
        state_dep_probs_manager = \
            DiscreteStateDependentParameterManager(
                matrix_state_dep_probs,
                self.meh.state_count,
                seed_age_for_time_slicing=self.meh.seed_age,
                list_time_slice_age_ends=self.meh.slice_age_ends[:-1]
            )

        self.prob_handler = DiscreteStateDependentProbabilityHandler(
            state_dep_probs_manager
        )

    def _initialize_str_representation(self) -> None:
        """Populate str_representation class member."""
        self.str_representation = str(self.meh)

        if self.prob_handler:
            self.str_representation += "\n" + str(self.prob_handler)

    def get_meh(self) -> MacroevolEventHandler:
        return self.meh

    def get_prob_handler(self) -> DiscreteStateDependentProbabilityHandler:
        return self.prob_handler

    def __str__(self) -> str:
        return self.str_representation

    def __repr__(self) -> str:
        return self.str_representation


class StateIntoPatternConverter:
    """
    Stash and machinery for checking and converting character
    compound-states into bit patterns, and vice-versa.

    For example, if we are thinking about regions as characters,
    then compound-state means different ranges (e.g., "A", "AB"),
    while 'state' is the number of different values a single
    character can take. In the case of 1-character-1-region, then
    the number of states is 2, because a species is either present
    or not in a region.

    Parameters:
        n_char (int): Number of characters.
        n_states_per_char (int): Number of states per character.
        n_states (int): Total number of states.
        int2set_dict (dict of str: str): Dictionary for converting a
            state coded as an integer into a bit set.
        set2int_dict (dict of str: str): Dictionary for converting a
            state coded as a bit set into an integer.
    """

    # e.g., for 2 states
    # Region | compound-state  | Bit pattern
    # --------------------------------------
    # Dead   |        0        | 000
    # A      |        1        | 100
    # B      |        2        | 010
    # C      |        3        | 001
    # AB     |        4        | 110
    # AC     |        5        | 101
    # BC     |        6        | 011
    # ABC    |        7        | 111

    n_char: int
    n_states_per_char: int
    n_states: int
    int2set_dict: ty.Dict[str, str]
    set2int_dict: ty.Dict[str, str]

    def __init__(self, n_characters: int, n_states_per_char: int) -> None:
        self.n_char = n_characters
        self.n_states_per_char = n_states_per_char
        self.n_states = int(n_states_per_char ** n_characters)

        self.int2set_dict = dict()
        self.set2int_dict = dict()

        # side-effect: initializes
        #     (i)  self.int2set_dict
        #     (ii) self.set2int_dict
        self._initialize_dicts()

    def _initialize_dicts(self) -> None:
        """Initialize and populate conversion dictionary members.
        
        This method non-recursively generates all bit patterns for a
        given number of characters and states per characters.

        Note that if we are thinking about regions as characters, then
        "compound state" means the range (e.g., "AB"), and state is the
        number of options per character (e.g., binary if present/absent
        in a region).

        Dictionary members are initialized as a side-effect.
        """

        # non-recursive method
        list_for_sorting: ty.List[ty.List[str]] = \
            [[] for i in range(self.n_char + 1)]
        for compound_state_idx in range(self.n_states):
            bit_pattern = [0 for i in range(self.n_char)]
            z = compound_state_idx

            for bit_idx in range(self.n_char - 1, -1, -1):
                # n_states_per_char here is present/absent for regions
                # (i.e., binary state)
                v = self.n_states_per_char ** bit_idx

                if z >= v:
                    bit_pattern[bit_idx] = 1
                    z -= v

                else:
                    bit_pattern[bit_idx] = 0

            _bit_pattern_str = "".join(str(b) for b in bit_pattern)
            _n_bits_on = sum(bit_pattern)
            list_for_sorting[_n_bits_on].append(_bit_pattern_str)

        # will store all bit patterns after sorting first by number of bits
        # that are on, and then by the decimal value underlying the bit pattern
        sorted_list_patterns: ty.List[str] = []
        for list_of_patterns in list_for_sorting:
            sorted_list_patterns.extend(list_of_patterns)

        # getting the other dict by reversing this one
        self.int2set_dict = \
            dict((str(idx), bit_pattern)
                 for idx, bit_pattern
                 in enumerate(sorted_list_patterns))
        for k, v in self.int2set_dict.items():
            self.set2int_dict[v] = k


if __name__ == "__main__":
    # can be called from calculation/
    # $ python3 discrete_sse.py
    #
    # can also be called from phylojunction/
    # $ python3 calculation/discrete_sse.py
    # or
    # $ python3 -m calculation.discrete_sse
    #
    # can also be called from VS Code, if open folder is phylojuction/

    # yule, one state, one epoch
    # arp = DiscreteStateDependentRate(1.0, MacroevolEvent.W_SPECIATION, name="lambda")
    # print(arp)
    # rates_t0_s0 = [ arp ]
    # matrix_atomic_rate_params = [ rates_t0_s0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice
    # state_dep_par_manager = DiscreteStateDependentParameterManager(matrix_atomic_rate_params, 1)

    # n_characters = 3 # regions A and B
    # n_states_per_char = 2 # presence/absence
    # svc = StateIntoPatternConverter(n_characters, n_states_per_char)

    # print(svc.int2set_dict)
    # print(svc.set2int_dict)

    n_characters = 8  # regions A and B
    n_states_per_char = 2  # presence/absence
    svc = StateIntoPatternConverter(n_characters, n_states_per_char)

    # print(svc.int2set_dict)
    for k, v in svc.int2set_dict.items():
        print(k + "\t" + v)
    # print(svc.set2int_dict)
