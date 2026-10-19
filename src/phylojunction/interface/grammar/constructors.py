"""Construct domain objects from evaluated values, without temporary DAG nodes."""
import math
import inspect
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


# Shapes describe calling conventions; implementations retain concrete value validation.
Scalar, Vector = "Scalar", "Vector"

CONSTRUCTORS = {
    'normal': (dn.DnNormal, [{
        'mean': Scalar, 'sd': Scalar,
    }]),
    'lognormal': (dn.DnLogNormal, [{
        'meanlog': Scalar, 'sdlog': Scalar, 'log_space': Scalar,
    }]),
    'exponential': (dn.DnExponential, [{
        'rate': Scalar, 'rate_parameterization': Scalar,
    }]),
    'gamma': (dn.DnGamma, [{
        'shape': Scalar, 'scale': Scalar, 'rate_parameterization': Scalar,
    }]),
    'unif': (dn.DnUnif, [{
        'min': Scalar, 'max': Scalar,
    }]),
    'categorical': (dn.DnCategorical, [{
        'cats': Vector, 'probs': Vector,
    }]),
    'beta': (dn.DnBeta, [{
        'alpha': Scalar, 'beta': Scalar,
    }]),
    'iid': (DnIID, [{
        'n': Scalar, 'distribution': Scalar,
    }]),
    'quasse': (DnQuaSSE, [{
        'birth_rate': Scalar, 'death_rate': Scalar, 'stop': Scalar, 'stop_value': Scalar, 'start_trait': Scalar, 
        'drift': Scalar, 'diffusion': Scalar, 'method': Scalar, 'k': Scalar, 'runtime_limit': Scalar, 
        'max_steps': Scalar, 'max_alive': Scalar, 'sampling_prob': Scalar, 'cond_surv': Scalar, 
        'cond_spn': Scalar, 'cond_obs_both_sides': Scalar, 'min_rec_taxa': Scalar, 'max_rec_taxa': Scalar, 
        'max_n_attempts': Scalar, 'dt_max': Scalar, 'max_nodes': Scalar, 'block_duration': Scalar, 
        'bridge_error': Scalar, 'fossil_rate': Scalar, 'start_distribution': Scalar,
    }]),
    'sse_rate': (sse_rate, [{
        'name': Scalar, 'value': Scalar, 'event': Scalar, 'states': Vector, 'epoch': Scalar,
    }]),
    'sse_prob': (sse_prob, [{
        'name': Scalar, 'value': Scalar, 'state': Scalar, 'epoch': Scalar,
    }]),
    'sse_stash': (sse_stash, [{
        'flat_rate_mat': Vector, 'n_states': Scalar, 'n_epochs': Scalar, 'seed_age': Scalar, 
        'epoch_age_ends': Vector, 'flat_prob_mat': Vector,
    }]),
    'discrete_sse': (discrete_sse, [{
        'stash': Scalar, 'stop': Scalar, 'stop_value': Scalar, 'start_state': Scalar, 'origin': Scalar, 
        'cond_spn': Scalar, 'cond_surv': Scalar, 'cond_obs_both_sides': Scalar, 'min_rec_taxa': Scalar, 
        'max_rec_taxa': Scalar, 'abort_at_alive_count': Scalar, 'eps': Scalar, 'runtime_limit': Scalar, 
        'max_n_attempts': Scalar,
    }]),
    'read_tree': (read_tree, [{
        'file_path': Scalar, 'string': Scalar, 'node_name_attr': Scalar, 'epoch_age_ends': Vector,
    }]),
    'map_attr': (map_attr, [{
        'tree': Scalar, 'fun': Scalar, 'tip_attr_file_path': Scalar, 'maps_file_path': Scalar, 
        'attr_name': Scalar, 'n_regions': Scalar, 'geo': Scalar,
    }]),
    'quasse_constant': (RATE_CONSTRUCTORS["quasse_constant"][0], [{
        'rate': Scalar, 'cap': Scalar,
    }]),
    'quasse_logistic': (RATE_CONSTRUCTORS["quasse_logistic"][0], [{
        'y0': Scalar, 'y1': Scalar, 'midpoint': Scalar, 'slope': Scalar, 'cap': Scalar,
    }]),
    'quasse_gaussian': (RATE_CONSTRUCTORS["quasse_gaussian"][0], [{
        'baseline': Scalar, 'center_rate': Scalar, 'center': Scalar, 'width': Scalar, 'cap': Scalar,
    }]),
    'quasse_step': (RATE_CONSTRUCTORS["quasse_step"][0], [{
        'left': Scalar, 'right': Scalar, 'threshold': Scalar, 'cap': Scalar,
    }]),
    'quasse_linear': (RATE_CONSTRUCTORS["quasse_linear"][0], [{
        'intercept': Scalar, 'slope': Scalar, 'cap': Scalar,
    }]),
    'quasse_skew_gaussian': (RATE_CONSTRUCTORS["quasse_skew_gaussian"][0], [{
        'baseline': Scalar, 'amplitude': Scalar, 'location': Scalar, 'width': Scalar, 'skew': Scalar, 
        'cap': Scalar,
    }]),
    'quasse_quadratic': (RATE_CONSTRUCTORS["quasse_quadratic"][0], [{
        'baseline': Scalar, 'strength': Scalar, 'center': Scalar, 'cap': Scalar,
    }]),
    'quasse_absolute': (RATE_CONSTRUCTORS["quasse_absolute"][0], [{
        'baseline': Scalar, 'strength': Scalar, 'center': Scalar, 'cap': Scalar,
    }]),
}

# The legacy map API accepts either a tree or a collection until read_maps is separated.
CONSTRUCTORS["map_attr"][1].append({**CONSTRUCTORS["map_attr"][1][0], "tree": Vector})


# Fail at registration if metadata drifts from a callable, rather than guessing missing shapes.
def validate_registry():
    for name, (function, signatures) in CONSTRUCTORS.items():
        parameters = set(inspect.signature(function).parameters)
        if not signatures or any(set(signature) != parameters for signature in signatures):
            raise ValueError(f"Incomplete argument shapes for {name}.")
        if any(shape not in (Scalar, Vector) for signature in signatures for shape in signature.values()):
            raise ValueError(f"Unknown argument shape for {name}.")


validate_registry()


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
    return CONSTRUCTORS[name][0](*args, **kwargs)
