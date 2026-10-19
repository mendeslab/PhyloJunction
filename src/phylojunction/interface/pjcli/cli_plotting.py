import typing as ty
from matplotlib.figure import Figure  # type: ignore
from matplotlib.axes import Axes  # type: ignore

# pj imports #
from phylojunction.interface.pysidegui.content_main_window \
    import ContentGUIMainWindow
from phylojunction.pgm.pgm import DirectedAcyclicGraph
from phylojunction.pgm.pgm import NodeDAG
import phylojunction.interface.cmdbox.cmd_parse as cmdp
import phylojunction.plotting.pj_organize as pjorg
import phylojunction.plotting.pj_draw as pjdraw
import phylojunction.readwrite.pj_read as pjread
import phylojunction.readwrite.pj_write as pjwrite
import phylojunction.data.tree as pjdt
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


def start_fig_and_axes(
        disabled_yticks: bool = True,
        disabled_xticks: bool = True,
        fig_width: int = 11,
        fig_height: int = 4.5) \
        -> ty.Tuple[Figure, Axes]:
    """
    """

    fig = Figure(figsize=(fig_width,fig_height))
    ax = fig.add_axes([0.075, 0.25, 0.6, 0.7])
    ax.patch.set_alpha(0.0)

    if disabled_xticks:
            ax.xaxis.set_ticks([])

    if disabled_yticks:
        ax.yaxis.set_ticks([])

    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    return fig, ax


# Export each requested execution and nested element using explicit, one-based paths.
def call_results_plot_cli(fig_dir, results, node_range_dict, fig_obj, fig_axes, prefix=''):
    from pathlib import Path
    from phylojunction.pgm.values import elements, path_text, plot_values
    for name, bounds in node_range_dict.items():
        if name not in results.model.name_node_dict:
            raise ValueError(f'Unknown model variable: {name}')
        start = bounds[0] if bounds else 1
        end = bounds[-1] if bounds else results.n_executions
        if not 1 <= start <= end <= results.n_executions:
            raise ValueError('Replicate range is out of bounds.')
        for execution in range(start, end + 1):
            for path, value in elements(results.value(name, execution - 1)):
                plot_values(fig_axes, [value])
                filename = ((prefix + '_') if prefix else '') + f'{name}_replicate{execution}{path_text(path)}.png'
                pjwrite.write_fig_to_file(str(Path(fig_dir) / filename), fig_obj)
