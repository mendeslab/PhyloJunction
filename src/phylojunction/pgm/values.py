"""Selection and summaries of values across model replicates."""
import numbers
import statistics
from phylojunction.data.tree import AnnotatedTree


# Preserve every nesting level, including empty lists, when assigning element paths.
def elements(value, path=()):
    if isinstance(value, list) and value:
        for i, item in enumerate(value, 1):
            yield from elements(item, (*path, i))
    else:
        yield path, value


def path_text(path):
    return ''.join(f'[{i}]' for i in path)


def shape(value):
    return tuple(shape(item) for item in value) if isinstance(value, list) else None


# Summary rows compare the same element across replicates, never elements within a run.
# Ragged vectors have no automatic correspondence; their values remain displayable/exportable.
def summaries(values):
    if len(values) < 2 or any(shape(v) != shape(values[0]) for v in values[1:]):
        return []
    rows = []
    for columns in zip(*(elements(v) for v in values)):
        path = columns[0][0]
        entries = [v for _, v in columns]
        metrics = {'value': entries}
        if all(isinstance(v, AnnotatedTree) for v in entries):
            stats = [v.get_stats_dict() for v in entries]
            metrics = {}
            for name in stats[0]:
                try:
                    metrics[name] = [float(v[name]) for v in stats]
                except (KeyError, ValueError, TypeError):
                    continue
        for metric, numbers_ in metrics.items():
            if all(isinstance(v, numbers.Real) and not isinstance(v, bool) for v in numbers_):
                rows.append({'@element': path_text(path), 'statistic': metric,
                             'mean': statistics.mean(numbers_), 'sd': statistics.stdev(numbers_)})
    return rows


# Draw one selected tree or numeric observations across replicates. Vectors must be
# selected explicitly by the caller; they are never silently pooled into a histogram.
def plot_values(axes, values, reconstructed=False):
    from phylojunction.pgm.pgm import plot_node_histogram, plot_blank
    if len(values) == 1 and isinstance(values[0], AnnotatedTree):
        values[0].plot_node(axes, draw_reconstructed=reconstructed)
    elif values and all(isinstance(v, numbers.Real) for v in values):
        plot_node_histogram(axes, values)
    else:
        plot_blank(axes)


# Comparison/coverage entries name an explicit component, shared across replicates.
def component_values(results, name):
    from phylojunction.interface.cmdbox.expressions import parse_expression, evaluate
    expression = parse_expression(name)
    return [evaluate(expression, lambda variable: results.value(variable, i), None)
            for i in range(results.n_executions)]

# Preserve tree statistic columns used by comparison plots; numeric values retain
# one row per replicate, including fixed values used in coverage experiments.
def component_frame(results, name):
    import pandas as pd
    values = component_values(results, name)
    if all(isinstance(v, AnnotatedTree) for v in values):
        rows = [{'root_age': v.root_age, 'origin_age': v.origin_age,
                 'n_total': v.n_extant_terminal_nodes + v.n_extinct_terminal_nodes,
                 'n_extant': v.n_extant_terminal_nodes, 'n_extinct': v.n_extinct_terminal_nodes,
                 'n_sa': v.n_sa_nodes} for v in values]
        frame = pd.DataFrame(rows)
    else:
        frame = pd.DataFrame({name: values})
    frame.insert(0, 'sample', range(1, len(values) + 1))
    frame.insert(1, 'replicate', 1)
    return frame

