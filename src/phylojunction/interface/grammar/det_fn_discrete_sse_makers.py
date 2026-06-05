import typing as ty

# pj imports
import phylojunction.utility.exception_classes as ec
import phylojunction.utility.helper_functions as pjh
import phylojunction.pgm.pgm as pgm
import phylojunction.calculation.discrete_sse as sseobj  # type: ignore

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


def make_DiscreteStateDependentRate(
    det_fn_name: str,
    det_fn_param_dict:
        ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
        -> sseobj.DiscreteStateDependentRate:
    """
    Create and return DiscreteStateDependentRate as prompted by
    deterministic function call.

    Args:
        det_fn_name (str): Name of the function being called
        det_fn_param_dict (dict): dictionary containing parameter
            strings as keys, and lists of either strings or NodeDAGs
            as value(s)

    Returns:
        Object holding the info about a discrete state-dependent rate
    """

    if not det_fn_param_dict:
        raise ec.ParseMissingSpecificationError(det_fn_name)

    event_type = None
    value: ty.List[float] = []
    the_states: ty.List[int] = []
    epoch_idx: int = 1

    ######################################
    # Checking for input correctness,    #
    # and preparing for instantiation of #
    # parameter objects                  #
    ######################################

    # val is a list
    for arg, val in det_fn_param_dict.items():
        if not val:
            raise ec.ParseMissingArgumentError(arg)

        if arg == "value":
            # val is a list of random variable objects
            # if type(val[0]) != str:
            if isinstance(val[0], pgm.NodeDAG):
                # need to declare cast_val separately so mypy won't complain
                cast_val1: ty.List[pgm.NodeDAG] = \
                    ty.cast(ty.List[pgm.NodeDAG], val)
                value = pgm.extract_value_from_dagnodes(cast_val1)

            # val is a list of strings
            elif isinstance(val[0], str):
                cast_val2: ty.List[str] = ty.cast(ty.List[str], val)
                value = [float(v) for v in cast_val2]

        elif arg == "event":
            if not (val[0].startswith("\"") or val[0].startswith("\'")) or \
                    not (val[0].endswith("\"") or val[0].endswith("\'")):
                raise ec.ParseFunctionArgError(
                    arg, ("Event name must start and end with either single "
                          "(\') or double (\") quotes."))

            if val[0] == "\"w_speciation\"" or val[0] == "\"speciation\"":
                event_type = sseobj.MacroevolEvent.W_SPECIATION

            elif val[0] == "\"bw_speciation\"":
                event_type = sseobj.MacroevolEvent.BW_SPECIATION

            elif val[0] == "\"asym_speciation\"":
                event_type = sseobj.MacroevolEvent.ASYM_SPECIATION

            elif val[0] == "\"extinction\"":
                event_type = sseobj.MacroevolEvent.EXTINCTION

            elif val[0] == "\"transition\"":
                event_type = sseobj.MacroevolEvent.ANAGENETIC_TRANSITION

            elif val[0] == "\"anc_sampling\"":
                event_type = sseobj.MacroevolEvent.ANCESTOR_SAMPLING

            else:
                raise ec.ParseFunctionArgError(
                    arg, ("Did not recognize event name. Please check "
                          "documentation for allowed events"))

        elif arg == "states":
            # need to declare cast_val separately
            # so mypy won't complain
            the_states = ty.cast(ty.List[int], val)

        elif arg == "epoch":
            if len(val) > 1:
                raise ec.ParseRequireSingleValueError(det_fn_name,
                                                      arg)

            try:
                # default is 1
                epoch_idx = int(val[0])  # start at 1

                if epoch_idx < 1:
                    raise ec.ParseRequirePositiveIntegerError(det_fn_name,
                                                              arg)

            except ValueError:
                raise ec.ParseRequireIntegerError(det_fn_name, arg)


    ###################################
    # Instantiating parameter objects #
    ###################################

    sse_rate_name = det_fn_param_dict["name"][0]

    if the_states:
        # TODO: deal with vectorization later
        if isinstance(sse_rate_name, str):
            # see 'invariance vs covariance' in mypy's doc for why list() is called here
            return sseobj.DiscreteStateDependentRate(
                name=sse_rate_name,
                val=list(value),
                event=event_type,
                states=the_states,
                epoch_idx=epoch_idx)

    else:
        # TODO: deal with vectorization later
        if isinstance(sse_rate_name, str):
            # see 'invariance vs covariance' in mypy's doc for why list() is called here
            return sseobj.DiscreteStateDependentRate(
                name=sse_rate_name,
                val=list(value),
                event=event_type,
                epoch_idx=epoch_idx)


def make_DiscreteStateDependentProbability(
    det_fn_name: str,
    det_fn_param_dict:
        ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
        -> sseobj.DiscreteStateDependentProbability:
    """
    Return SSE probability as prompted by deterministic function call.

    Args:
        det_fn_name (str): Name of the function being called.
        det_fn_param_dict (dict): Dictionary containing parameter
            strings as keys, and lists of either strings or NodeDAGs
            as value(s).

    Returns:
        DiscreteStateDependentProbability: Object holding the info
            about a discrete state-dependent probability.
    """

    if not det_fn_param_dict:
        raise ec.ParseMissingSpecificationError(det_fn_name)

    value: ty.List[float] = []
    state: ty.Optional[int] = 0
    epoch_idx: int = 1

    ######################################
    # Checking for input correctness,    #
    # and preparing for instantiation of #
    # parameter objects                  #
    ######################################

    # val is a list
    for arg, val in det_fn_param_dict.items():
        if not val:
            raise ec.ParseMissingArgumentError(arg)

        if arg == "value":
            # val is a list of random variable objects
            # if type(val[0]) != str:
            if isinstance(val[0], pgm.NodeDAG):
                # need to declare cast_val separately so mypy won't complain
                cast_val1: ty.List[pgm.NodeDAG] = \
                    ty.cast(ty.List[pgm.NodeDAG], val)
                value = pgm.extract_value_from_dagnodes(cast_val1)

            # val is a list of strings
            elif isinstance(val[0], str):
                cast_val2: ty.List[str] = ty.cast(ty.List[str], val)
                value = [float(v) for v in cast_val2]

        elif arg == "state":
            # need to declare cast_val separately
            # so mypy won't complain
            if len(val) > 1:
                raise ec.ParseRequireSingleValueError(
                    det_fn_name,
                    arg
                )

            state = int(val[0])

        elif arg == "epoch":
            if len(val[0]) > 1:
                raise ec.ParseRequireSingleValueError(det_fn_name,
                                                      arg)
            try:
                # default is 1
                epoch_idx = int(val[0])  # start at 1

                if epoch_idx < 1:
                    raise ec.ParseRequirePositiveIntegerError(det_fn_name,
                                                              arg)

            except ValueError:
                raise ec.ParseRequireIntegerError(det_fn_name, arg)

    ##################################
    # Istantiating parameter objects #
    ##################################

    state_dep_prob_name = det_fn_param_dict["name"][0]

    if isinstance(state_dep_prob_name, str):
        if state is None:
            return sseobj.DiscreteStateDependentProbability(
                name=state_dep_prob_name,
                val=list(value))

        else:
            return sseobj.DiscreteStateDependentProbability(
                name=state_dep_prob_name,
                val=list(value),
                state=state,
                epoch_idx=epoch_idx)


def make_SSEStash(
    det_fn_name: str,
    det_fn_param_dict:
        ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
        -> sseobj.SSEStash:
    """
    Create SSEStash as prompted by deterministic function call

    Args:
        det_fn_name (str): Name of the function being called
        det_fn_param_dict (dict): dictionary containing parameter
            strings as keys, and lists of either strings or NodeDAGs

    Return:
        Object holding all discrete state-dependent rates and
            probabilities
    """

    n_states: int = 1
    n_time_slices: int = 1
    time_slice_age_ends: ty.List[float] = list()
    seed_age_for_time_slicing: ty.Optional[float] = None
    flat_state_dep_rate_mat: ty.List[pgm.DeterministicNodeDAG] = []
    flat_state_dep_prob_mat: ty.List[pgm.DeterministicNodeDAG] = []

    #############################################
    # Reading all arguments and checking health #
    #############################################

    # val is a list of strings or nodes
    for arg, val in det_fn_param_dict.items():
        first_element = val[0]
        extracted_value = first_element  # can be scalar or container

        if isinstance(first_element, pgm.NodeDAG) \
                and len(val) == 1:
            extracted_value = first_element.value

        # so mypy won't complain
        # if isinstance(val[0], str):
        if arg in ("n_states", "n_epochs", "seed_age"):
            # none of the above can be more than one value
            if isinstance(extracted_value, list) and \
                    len(extracted_value) > 1:
                # should provide only one number of states
                raise ec.ParseRequireSingleValueError(det_fn_name, arg)

        if arg == "n_states":
            try:
                n_states = int(extracted_value)

            except ValueError:
                raise ec.ParseRequireIntegerError(det_fn_name, arg)

        elif arg == "n_epochs":
            try:
                n_time_slices = int(extracted_value)

            except ValueError:
                raise ec.ParseRequireIntegerError(det_fn_name, arg)

        elif arg == "seed_age":
            try:
                seed_age_for_time_slicing = float(extracted_value)

            except ValueError:
                raise ec.ParseRequireNumericError(det_fn_name, arg)

        elif arg == "epoch_age_ends":
            if len(val) != (n_time_slices - 1):
                raise ec.IncorrectDimensionError(
                    "epoch_age_ends",
                    len(val),
                    exp_len=(n_time_slices - 1))

            try:
                time_slice_age_ends = \
                    [float(v) for v in val if isinstance(v, str)]

            except ValueError:
                raise ec.ParseRequireNumericError(det_fn_name, arg)

        elif arg == "flat_prob_mat":
            if det_fn_param_dict["flat_prob_mat"]:
                flat_state_dep_prob_mat = \
                    [v for v in det_fn_param_dict["flat_prob_mat"]
                        if isinstance(v, pgm.DeterministicNodeDAG)]

            # total number of rates has to be divisible by number of slices
            # if len(flat_state_dep_prob_mat) % n_time_slices != 0:
            #     raise ec.IncorrectDimensionError(
            #         arg, len(flat_state_dep_prob_mat))

        elif arg == "flat_rate_mat":
            if det_fn_param_dict["flat_rate_mat"]:
                # list of NodeDAG's
                flat_state_dep_rate_mat = \
                    [v for v in det_fn_param_dict["flat_rate_mat"]
                        if isinstance(v, pgm.DeterministicNodeDAG)]

                # total number of rates has to be divisible by number of slices
                # if len(flat_state_dep_rate_mat) % n_time_slices != 0:
                #     raise ec.IncorrectDimensionError(
                #         arg, len(flat_state_dep_rate_mat))

            else:
                raise ec.ParseMissingParameterError(arg)

    ##################################
    # Putting rates in their manager #
    ##################################

    matrix_state_dep_rates: \
        ty.List[ty.List[sseobj.DiscreteStateDependentRate]] = \
            [[] for i in range(n_time_slices)] # []

    # print("running new code for initializing matrix_state_dep_rates")
    for state_dep_rate_det_nd in flat_state_dep_rate_mat:
        if isinstance(state_dep_rate_det_nd.value,
                      sseobj.DiscreteStateDependentRate):
            if state_dep_rate_det_nd.value.epoch_idx > n_time_slices:
                raise ec.ParseInvalidArgumentError(
                    "epoch",
                    "\'" + str(state_dep_rate_det_nd.value.epoch_idx) + "\'",
                    ("One of the provided rates was placed in an "
                     "epoch that exceeds the specified number of "
                     "epochs."))

            # -1 for offset
            matrix_state_dep_rates[state_dep_rate_det_nd.value.epoch_idx - 1] \
                .append(state_dep_rate_det_nd.value)
            
    # print(matrix_state_dep_rates)

    # DEPRECATED: initial method expected user to specify the same
    # number of rates for each and every epoch
    #
    # populating state-dependent rate matrix #
    # total_n_rate = len(flat_state_dep_rate_mat)
    # n_rates_per_slice = int(total_n_rate / n_time_slices)
    # for i in range(0, total_n_rate, n_rates_per_slice):
    #     state_dep_rates: ty.List[sseobj.DiscreteStateDependentRate] = []
    #     for state_dep_rate_det_nd in \
    #             flat_state_dep_rate_mat[i:(i + n_rates_per_slice)]:

    #         # so mypy won't complain
    #         if isinstance(state_dep_rate_det_nd.value,
    #                       sseobj.DiscreteStateDependentRate):
    #             state_dep_rates.append(state_dep_rate_det_nd.value)

    #     # appending a list of DiscreteStateDependentRate
    #     # 1D: time slices, 2D: state-dep rate list
    #     matrix_state_dep_rates.append(state_dep_rates)

    state_dep_rates_manager = \
        sseobj.DiscreteStateDependentParameterManager(
            matrix_state_dep_rates,
            n_states,
            seed_age_for_time_slicing=seed_age_for_time_slicing,
            list_time_slice_age_ends=time_slice_age_ends
        )

    ##################################
    # Putting probs in their manager #
    ##################################

    matrix_state_dep_probs: \
        ty.List[ty.List[sseobj.DiscreteStateDependentProbability]] = []
    expected_n_prob_per_slice = n_states * n_time_slices
    total_n_prob = n_states * n_time_slices  # initialize to expected
    n_prob_per_slice = int(total_n_prob / n_time_slices)

    # must provide one probability parameter per state
    # per time slice!
    probs_were_provided = False
    if len(flat_state_dep_prob_mat) > 0:
        if n_prob_per_slice != expected_n_prob_per_slice:
            # TODO: write exception
            pass

        probs_were_provided = True

    state_dep_probs_manager: \
        sseobj.DiscreteStateDependentParameterManager = None

    # populating state-dependent prob matrix #
    #
    # this iteration here is effectively over
    # time-slices if dimensions are indeed all
    # correct
    if probs_were_provided:
        for i in range(0, total_n_prob, n_prob_per_slice):
            state_dep_probs: \
                ty.List[sseobj.DiscreteStateDependentProbability] = []

            for state_dep_prob_det_nd in \
                    flat_state_dep_prob_mat[i:(i + n_prob_per_slice)]:
                # so mypy won't complain
                if isinstance(state_dep_prob_det_nd.value,
                              sseobj.DiscreteStateDependentProbability):
                    state_dep_probs.append(state_dep_prob_det_nd.value)

            # appending a list of DiscreteStateDependentRate
            # 1D: time slices, 2D: state-dep rate list
            matrix_state_dep_probs.append(state_dep_probs)

        state_dep_probs_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                matrix_state_dep_probs,
                n_states,
                seed_age_for_time_slicing=seed_age_for_time_slicing,
                list_time_slice_age_ends=time_slice_age_ends
            )

        return \
            sseobj.SSEStash(
                sseobj.MacroevolEventHandler(state_dep_rates_manager),
                sseobj.DiscreteStateDependentProbabilityHandler(
                    state_dep_probs_manager)
            )

        # # else here means we will have i being
        # # one integer per probability (per state
        # # per slice)
        # #
        # # we will make all probabilities 1.0 by
        # # default
        # else:
        #     st = i % n_prob_per_slice
        #     state_dep_prob = sseobj.DiscreteStateDependentProbability(
        #         name="rho" + str(st) + "_t" + str(i // n_prob_per_slice),
        #         val=1.0, state=st)

        #     # so mypy won't complain
        #     if isinstance(state_dep_prob,
        #         sseobj.DiscreteStateDependentProbability):
        #         state_dep_probs.append(state_dep_prob)

        #     # appending a list of DiscreteStateDependentProbability
        #     # 1D: time slices, 2D: state-dep prob list
        #     matrix_state_dep_probs.append(state_dep_probs)

        #     state_dep_probs_manager = \
        #         sseobj.DiscreteStateDependentParameterManager(
        #             matrix_state_dep_probs,
        #             n_states,
        #             seed_age_for_time_slicing=seed_age_for_time_slicing,
        #             list_time_slice_age_ends=time_slice_age_ends
        #         )

    # state_dep_probs_manager = \
    #     sseobj.DiscreteStateDependentParameterManager(
    #         # TODO
    #         n_states,
    #         seed_age_for_time_slicing=seed_age_for_time_slicing,
    #         list_time_slice_age_ends=time_slice_age_ends
    #     )

    # return sseobj.MacroevolEventHandler(
    #     state_dep_rates_manager,
    #     state_dep_probs_manager,

    # )

    return \
        sseobj.SSEStash(
            sseobj.MacroevolEventHandler(state_dep_rates_manager)
        )
