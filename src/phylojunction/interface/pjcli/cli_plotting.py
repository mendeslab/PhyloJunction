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


def selected_node_plot_cli(
        fig_dir: str,
        fig_obj: Figure,
        fig_axes: Axes,
        node_dag: NodeDAG,
        prefix: str = "",
        sample_idx: int = 0,
        repl_idx: int = 0,
        repl_size: int = 1) -> None:
    """
    Plot DAG node on provided Axes object (fig_axes), intended
    to be scoped to pj_cli.execute_pj_script() 
    then update canvas with new plot
    """

    outfile_path = fig_dir

    if prefix:
        outfile_path += prefix + "_"

    outfile_path += node_dag.node_name + str(sample_idx + 1) \
        + "_" + str(repl_idx + 1)
    
    print("Writing", outfile_path)
    
    # if stochastic or constant, value will be list
    if isinstance(node_dag.value, list):
        # if a tree
        if isinstance(node_dag.value[0], pjdt.AnnotatedTree):
            node_dag.plot_node(
                fig_axes,
                sample_idx=sample_idx,
                repl_idx=repl_idx,
                repl_size=repl_size)

        # when not a tree
        else:
            node_dag.plot_node(
                fig_axes,
                sample_idx=sample_idx,
                repl_size=repl_size)

        pjwrite.write_fig_to_file(outfile_path, fig_obj)

    # if deterministic, not subscriptable
    # except Exception as e:
    #     print("(cli_plotting.py) An error occurred: ", type(e).__name__, " - ", e)

    # deterministic node, value is an Object
    # else:
    #     print(type(node_dag.value))


def call_node_plot_cli(
        fig_dir: str,
        dag_obj: DirectedAcyclicGraph,
        n_samples: int,
        node_range_dict: ty.Dict[str, ty.Tuple[int]],
        fig_obj: Figure,
        fig_axes: Axes,
        prefix: str = "") -> None:
    """
    """
 
    for node_name, range_tup in node_range_dict.items():
        if node_name in dag_obj.name_node_dict:
            node_dag = dag_obj.get_node_dag_by_name(node_name)
            repl_size = node_dag.repl_size
            start_idx = 0
            end_idx = 1

            if len(range_tup) == 0:
                end_idx = n_samples

            elif len(range_tup) == 1:
                start_idx = range_tup[0] - 1 # offset
                end_idx = start_idx + 1

            elif len(range_tup) == 2:
                start_idx, end_idx = range_tup
                start_idx -= 1
                # end_idx += 1

            # debugging
            # print("start_idx = " + str(start_idx))
            # print("end_idx = " + str(end_idx))
            # print("repl_size = " + str(repl_size))

            for sample_idx in range(start_idx, end_idx):
                for repl_idx in range(repl_size):
                    try:
                        selected_node_plot_cli(
                            fig_dir,
                            fig_obj,
                            fig_axes,
                            node_dag,
                            prefix=prefix,
                            sample_idx=sample_idx,
                            repl_idx=repl_idx,
                            repl_size=repl_size
                        )

                    except IndexError:
                        raise ec.PJCLISampleOutOfRangeError(
                            str(start_idx) + "-" + str(end_idx))

        else:
            # raise ec.PJCLINodeNotInPGM(node_name)
            raise RuntimeError