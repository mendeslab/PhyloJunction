"""Convert QuaSSE grammar arguments to distribution parameters."""

from phylojunction.distribution.dn_quasse import DnQuaSSE
from phylojunction.pgm import pgm
from phylojunction.utility import exception_classes as ec


# Extract rate objects without string conversion; numeric model parameters retain vectors.
def make_quasse(spec):
    for name in ("birth_rate", "death_rate", "stop", "stop_value"):
        if name not in spec:
            raise ec.ParseMissingParameterError(name)
    args = {}
    vectors = {"start_trait", "drift", "diffusion", "stop_value", "sampling_prob"}
    try:
        for name, raw in spec.items():
            if name == "clamp":
                continue
            if name in ("birth_rate", "death_rate", "fossil_rate"):
                values = []
                for item in raw:
                    value = item.value if isinstance(item, pgm.NodeDAG) else item
                    values.extend(value if isinstance(value, list) else [value])
                args[name] = values
                continue
            values = pgm.extract_vals_as_str_from_node_dag(raw)
            if name in vectors:
                args[name] = [float(v) for v in values]
            else:
                if len(values) != 1:
                    raise ValueError(f"{name} requires one value.")
                if name in ("cond_surv", "cond_spn", "cond_obs_both_sides"):
                    flag = values[0].strip('"').lower()
                    if flag not in ("true", "false", "t", "f"):
                        raise ValueError(f"{name} must be true or false.")
                    args[name] = flag in ("true", "t")
                    continue
                args[name] = values[0].strip('"') if name in ("method", "stop") else float(values[0])
        return DnQuaSSE(**args)
    except (ValueError, TypeError) as error:
        raise ec.ParseDnInitFailError("quasse", str(error)) from error
