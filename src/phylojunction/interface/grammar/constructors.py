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


# Scalar arguments must not silently broadcast vectors or accept boolean values.
def scalar(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError("Expected a finite numeric scalar.")
    return value


# Retain documented string flags as well as boolean literals during script migration.
def flag(value):
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.lower() in ('true', 'false', 't', 'f'):
        return value.lower() in ('true', 't')
    raise ValueError("Expected true or false.")


# These constructors temporarily pass singleton parameters to the old Python API.
# Remove that conversion when legacy distribution consumers have been migrated.
def normal(mean, sd):
    mean, sd = scalar(mean), scalar(sd)
    if sd <= 0:
        raise ValueError("Normal sd must be positive.")
    return dn.DnNormal(1, 1, [mean], [sd], None)


# Log-space controls the interpretation of the location, retaining existing semantics.
def lognormal(meanlog, sdlog, log_space=True):
    meanlog, sdlog = scalar(meanlog), scalar(sdlog)
    log_space = flag(log_space)
    if sdlog <= 0 or (not log_space and meanlog <= 0):
        raise ValueError("Lognormal sdlog and natural-space scale must be positive.")
    return dn.DnLogNormal(1, 1, [meanlog], [sdlog], log_space)


# Keep the existing rate/scale switch; the distribution performs the conversion.
def exponential(rate, rate_parameterization=True):
    rate = scalar(rate)
    if rate <= 0:
        raise ValueError("Exponential rate or scale must be positive.")
    return dn.DnExponential(1, 1, [rate], flag(rate_parameterization))


# Gamma uses scale by default, matching the existing constructor.
def gamma(shape, scale, rate_parameterization=False):
    shape, scale = scalar(shape), scalar(scale)
    if shape <= 0 or scale <= 0:
        raise ValueError("Gamma shape and scale/rate must be positive.")
    return dn.DnGamma(1, 1, [shape], [scale], flag(rate_parameterization))


# Uniform bounds define one distribution, not a vector of parameter settings.
def unif(min, max):
    min, max = scalar(min), scalar(max)
    if min >= max:
        raise ValueError("Uniform min must be less than max.")
    return dn.DnUnif(1, 1, [min], [max])


# Category vectors are mathematical parameters, unlike the removed batch dimensions.
def categorical(cats, probs):
    if not isinstance(cats, list) or not isinstance(probs, list) or not cats or len(cats) != len(probs):
        raise ValueError("Categories and probabilities must be nonempty vectors of equal length.")
    if any(scalar(p) < 0 for p in probs) or sum(probs) <= 0:
        raise ValueError("Category weights must be nonnegative with a positive sum.")
    return dn.DnCategorical(1, 1, [cats], [probs])


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
    return DnSSE(stash, origin=flag(origin), start_states_list=[integer(start_state)],
                 stop=stop, stop_value=[scalar(stop_value)], condition_on_speciation=flag(cond_spn),
                 condition_on_survival=flag(cond_surv),
                 condition_on_obs_both_sides_root=flag(cond_obs_both_sides),
                 min_rec_taxa=integer(min_rec_taxa), max_rec_taxa=integer(max_rec_taxa),
                 abort_at_alive_count=integer(abort_at_alive_count, 1), epsilon=scalar(eps),
                 runtime_limit=scalar(runtime_limit), max_n_failed_attempts=integer(max_n_attempts, 1))


# Reading one or several trees always returns a list; [1] explicitly selects one tree.
def read_tree(file_path=None, string=None, node_name_attr='', epoch_age_ends=None):
    if (file_path is None) == (string is None):
        raise ValueError("Specify exactly one of file_path or string.")
    return CtFnTreeReader(1, 1, file_path if file_path is not None else string,
                          node_name_attr, file_path is not None,
                          [] if epoch_age_ends is None else epoch_age_ends).generate()


CONSTRUCTORS = {
    'normal': normal, 'lognormal': lognormal, 'exponential': exponential,
    'gamma': gamma, 'unif': unif, 'categorical': categorical,
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
    return CONSTRUCTORS[name](*args, **kwargs)
