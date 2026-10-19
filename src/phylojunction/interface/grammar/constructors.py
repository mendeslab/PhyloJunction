"""Construct domain objects from evaluated values, without temporary DAG nodes."""
import math
from phylojunction.distribution import dn_parametric as dn
from phylojunction.distribution.dn_iid import DnIID
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


CONSTRUCTORS = {
    'normal': normal, 'lognormal': lognormal, 'exponential': exponential,
    'gamma': gamma, 'unif': unif, 'categorical': categorical,
    'beta': dn.DnBeta, 'iid': DnIID,
    **{name: constructor for name, (constructor, _) in RATE_CONSTRUCTORS.items()},
}


# Only registered constructors are callable; Python globals and attributes are inaccessible.
def construct(name, *args, **kwargs):
    if name not in CONSTRUCTORS:
        if name.endswith('_distribution'):
            raise ValueError("Use a distribution constructor such as normal(mean, sd) directly.")
        raise ValueError(f"Unknown constructor: {name}")
    if ('n' in kwargs and name != 'iid') or 'nr' in kwargs:
        raise ValueError("n/nr no longer specify sampling batches. Use iid(n, distribution) or --runs.")
    return CONSTRUCTORS[name](*args, **kwargs)
