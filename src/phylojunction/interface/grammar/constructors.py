"""Construct domain objects from evaluated values, without temporary DAG nodes."""
import math
from phylojunction.calculation import discrete_sse as sse
from phylojunction.distribution.dn_discrete_sse import DnSSE
from phylojunction.constant_function.ct_function_read_tree import CtFnTreeReader
from phylojunction.interface.grammar.det_fn_map_attribute import map_attr
from phylojunction.distribution import dn_parametric as dn
from phylojunction.distribution.dn_iid import DnIID
from phylojunction.distribution.dn_quasse import DnQuaSSE
from phylojunction.interface.grammar.det_fn_quasse_makers import RATE_CONSTRUCTORS


# Domain factories share the scalar/flag validation used by distribution constructors.
def scalar(value):
    return dn.finite(value, 'parameter')


def flag(value):
    return dn.boolean(value)


# Integer controls are distinct from real parameters; do not silently round them.
def integer(value, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(f"Expected an integer at least {minimum}.")
    return value


# Rates hold one scalar value per execution; their states remain genuine vector arguments.
def sse_rate(name, value, event, states=None, epoch=1):
    events = {'w_speciation': sse.MacroevolEvent.W_SPECIATION,
              'speciation': sse.MacroevolEvent.W_SPECIATION,
              'bw_speciation': sse.MacroevolEvent.BW_SPECIATION,
              'asym_speciation': sse.MacroevolEvent.ASYM_SPECIATION,
              'extinction': sse.MacroevolEvent.EXTINCTION,
              'transition': sse.MacroevolEvent.ANAGENETIC_TRANSITION,
              'anc_sampling': sse.MacroevolEvent.ANCESTOR_SAMPLING}
    if event not in events:
        raise ValueError(f"Unknown SSE event: {event}")
    states = [] if states is None else [integer(state) for state in states]
    return sse.DiscreteStateDependentRate(scalar(value), events[event], name, states, integer(epoch, 1))


def sse_prob(name, value, state=0, epoch=1):
    return sse.DiscreteStateDependentProbability(scalar(value), name, integer(state), integer(epoch, 1))


# Put each parameter in its declared epoch, retaining the managers' domain validation.
# Epoch counts describe the biological process, not IID vectors or model executions.
def sse_stash(flat_rate_mat, n_states=1, n_epochs=1, seed_age=None,
              epoch_age_ends=None, flat_prob_mat=None):
    n_epochs, n_states = integer(n_epochs, 1), integer(n_states, 1)
    ends = [] if epoch_age_ends is None else [scalar(age) for age in epoch_age_ends]
    if len(ends) != n_epochs - 1:
        raise ValueError("epoch_age_ends must contain n_epochs - 1 ages.")
    managers = []
    for parameters, expected in ((flat_rate_mat, sse.DiscreteStateDependentRate),
                                  (flat_prob_mat, sse.DiscreteStateDependentProbability)):
        if parameters is None:
            managers.append(None)
            continue
        matrix = [[] for _ in range(n_epochs)]
        for parameter in parameters:
            if not isinstance(parameter, expected):
                raise TypeError(f"Expected {expected.__name__} objects.")
            if not 1 <= parameter.epoch_idx <= n_epochs:
                raise ValueError("Parameter epoch exceeds the declared number of epochs.")
            matrix[parameter.epoch_idx - 1].append(parameter)
        managers.append(sse.DiscreteStateDependentParameterManager(
            matrix, n_states, seed_age_for_time_slicing=seed_age, list_time_slice_age_ends=ends))
    return sse.SSEStash(sse.MacroevolEventHandler(managers[0]),
                       None if managers[1] is None else sse.DiscreteStateDependentProbabilityHandler(managers[1]))


# Adapt language parameter names to the established simulator without batch dimensions.
def discrete_sse(stash, stop, stop_value, start_state=0, origin=True, cond_spn=False,
                 cond_surv=True, cond_obs_both_sides=False, min_rec_taxa=0,
                 max_rec_taxa=int(1e12), abort_at_alive_count=int(1e12), eps=1e-12,
                 runtime_limit=300, max_n_attempts=200):
    return DnSSE(stash, origin=flag(origin), start_state=integer(start_state),
                 stop=stop, stop_value=scalar(stop_value), condition_on_speciation=flag(cond_spn),
                 condition_on_survival=flag(cond_surv),
                 condition_on_obs_both_sides_root=flag(cond_obs_both_sides),
                 min_rec_taxa=integer(min_rec_taxa), max_rec_taxa=integer(max_rec_taxa),
                 abort_at_alive_count=integer(abort_at_alive_count, 1), epsilon=scalar(eps),
                 runtime_limit=scalar(runtime_limit), max_n_failed_attempts=integer(max_n_attempts, 1))


# Reading one or several trees always returns a list; [1] explicitly selects one tree.
def read_tree(file_path=None, string=None, node_name_attr='', epoch_age_ends=None):
    if (file_path is None) == (string is None):
        raise ValueError("Specify exactly one of file_path or string.")
    return CtFnTreeReader(file_path if file_path is not None else string,
                          node_name_attr, file_path is not None,
                          [] if epoch_age_ends is None else epoch_age_ends).generate()


CONSTRUCTORS = {
    'normal': dn.DnNormal, 'lognormal': dn.DnLogNormal, 'exponential': dn.DnExponential,
    'gamma': dn.DnGamma, 'unif': dn.DnUnif, 'categorical': dn.DnCategorical,
    'beta': dn.DnBeta, 'iid': DnIID, 'quasse': DnQuaSSE,
    'sse_rate': sse_rate, 'sse_prob': sse_prob, 'sse_stash': sse_stash,
    'discrete_sse': discrete_sse, 'read_tree': read_tree, 'map_attr': map_attr,
    **{name: constructor for name, (constructor, _) in RATE_CONSTRUCTORS.items()},
}


# Only registered constructors are callable; Python globals and attributes are inaccessible.
def construct(constructor_name, /, *args, **kwargs):
    name = constructor_name
    if name not in CONSTRUCTORS:
        if name.endswith('_distribution'):
            raise ValueError("Use a distribution constructor such as normal(mean, sd) directly.")
        raise ValueError(f"Unknown constructor: {name}")
    if ('n' in kwargs and name != 'iid') or 'nr' in kwargs:
        raise ValueError("n/nr no longer specify sampling batches. Use iid(n, distribution) or --runs.")
    # Script booleans historically allow quoted true/false; domain Python APIs use bool.
    if name == 'quasse':
        for key in ('cond_surv', 'cond_spn', 'cond_obs_both_sides'):
            if key in kwargs:
                kwargs[key] = flag(kwargs[key])
    return CONSTRUCTORS[name](*args, **kwargs)
