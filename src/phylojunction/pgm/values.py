"""Selection and summaries of values across complete model executions."""
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


# Summary rows compare the same element across executions, never elements within a run.
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


# Draw one selected tree or numeric observations across executions. Vectors must be
# selected explicitly by the caller; they are never silently pooled into a histogram.
def plot_values(axes, values, reconstructed=False):
    from phylojunction.pgm.pgm import plot_node_histogram, plot_blank
    if len(values) == 1 and isinstance(values[0], AnnotatedTree):
        values[0].plot_node(axes, draw_reconstructed=reconstructed)
    elif values and all(isinstance(v, numbers.Real) for v in values):
        plot_node_histogram(axes, values)
    else:
        plot_blank(axes)
