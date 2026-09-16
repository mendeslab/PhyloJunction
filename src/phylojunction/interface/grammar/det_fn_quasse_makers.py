"""Construct vectorized continuous-trait rate values for the DAG."""

from phylojunction.calculation.continuous_sse import LogisticRate
from phylojunction.pgm import pgm
from phylojunction.utility import exception_classes as ec


# Broadcast the four parameter vectors before constructing scalar rate objects.
def make_logistic_rate(spec):
    names = ("y0", "y1", "midpoint", "slope")
    for name in names:
        if name not in spec:
            raise ec.ParseMissingParameterError(name)
    try:
        values = [list(map(float, pgm.extract_vals_as_str_from_node_dag(spec[name]))) for name in names]
        size = max(map(len, values))
        if any(len(v) not in (1, size) for v in values):
            raise ValueError("Logistic parameter lengths must be one or a common length.")
        return [LogisticRate(*(v[0 if len(v) == 1 else i] for v in values)) for i in range(size)]
    except (ValueError, TypeError) as error:
        raise ec.ParseDetFnInitFailError("quasse_logistic", str(error)) from error
