"""Construct vectorized continuous-trait rate values for the DAG."""

from phylojunction.calculation.continuous_sse import (
    ConstantRate, LogisticRate, GaussianRate, StepRate, LinearRate,
    SkewGaussianRate, QuadraticRate, AbsoluteRate,
)
from phylojunction.pgm import pgm
from phylojunction.utility import exception_classes as ec


# Keep constructor dispatch and the grammar's allowed parameters in one place.
RATE_CONSTRUCTORS = {
    "quasse_constant": (ConstantRate, ("rate",)),
    "quasse_logistic": (LogisticRate, ("y0", "y1", "midpoint", "slope")),
    "quasse_gaussian": (GaussianRate, ("baseline", "center_rate", "center", "width")),
    "quasse_step": (StepRate, ("left", "right", "threshold")),
    "quasse_linear": (LinearRate, ("intercept", "slope")),
    "quasse_skew_gaussian": (SkewGaussianRate, ("baseline", "amplitude", "location", "width", "skew")),
    "quasse_quadratic": (QuadraticRate, ("baseline", "strength", "center")),
    "quasse_absolute": (AbsoluteRate, ("baseline", "strength", "center")),
}


# Broadcast parameters, including an optional cap, before constructing scalar rate objects.
def make_quasse_rate(det_fn_id, spec):
    constructor, names = RATE_CONSTRUCTORS[det_fn_id]
    for name in names:
        if name not in spec:
            raise ec.ParseMissingParameterError(name)
    if "cap" in spec:
        names = (*names, "cap")
    try:
        values = [list(map(float, pgm.extract_vals_as_str_from_node_dag(spec[name]))) for name in names]
        size = max(map(len, values))
        if any(len(v) not in (1, size) for v in values):
            raise ValueError("Rate parameter lengths must be one or a common length.")
        return [constructor(**{name: v[0 if len(v) == 1 else i] for name, v in zip(names, values)})
                for i in range(size)]
    except (ValueError, TypeError) as error:
        raise ec.ParseDetFnInitFailError(det_fn_id, str(error)) from error
