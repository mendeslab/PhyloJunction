import typing as ty
import matplotlib.pyplot as plt  # type: ignore
import seaborn as sns  # type: ignore
import pandas as pd  # type: ignore
from matplotlib.figure import Figure  # type: ignore
from tabulate import tabulate  # type: ignore

# pj imports
import phylojunction.readwrite.pj_read as pjread
import phylojunction.readwrite.pj_write as pjwrite
import phylojunction.pgm.pgm as pgm
import phylojunction.interface.cmdbox.cmd_parse as cmdp
import phylojunction.plotting.pj_organize as pjorg


def plot_violins(fig: Figure,
                 ax: plt.Axes,
                 df: pd.DataFrame,
                 x: str,
                 y: str,
                 xlab: ty.Optional[str] = None,
                 ylab: ty.Optional[str] = None) -> None:
    """Draw violin plots, one variable and 2+ factors

    Args:
        fig (matplotlib.Figure): Figure object
        ax (plt.Axes): Axes object
        df (pd.DataFrame): pandas DataFrame with values for a variable
            (column 1) in two or more scenarios (column 2), i.e., two
            or more factors (e.g., different simulators)
        xlab (str): x-axis label. Defaults to None.
        ylab (str): y-axis label. Defaults to None.

    Returns:
        None
    """

    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)

    custom_palette = sns.color_palette("Spectral")
    sns.set_palette(custom_palette)
    sns.violinplot(x=x, y=y, data=df, ax=ax, linewidth=0.0)

    if xlab is not None:
        ax.set_xlabel(xlab)

    if ylab is not None:
        ax.set_ylabel(ylab)

    fig.canvas.draw()


def plot_intervals(
        fig_obj: plt.Figure,
        axes_obj: plt.Axes,
        df: pd.DataFrame,
        x: str,
        y: str,
        xlab: ty.Optional[str] = None,
        ylab: ty.Optional[str] = "Posterior mean") -> None:
    """Draw coverage plot on provided plt.Figure instance

    Args:
        fig_obj (matplotlib.Figure): Figure object
        axes_obj (plt.Axes): Axes object
        df (pd.DataFrame): pandas dataframe holding interval
            information
        x (str): Label for x-axis on pandas dataframe.
        y (str, optional): Label for y-axis on pandas dataframe.
            Defaults to "posterior_mean".
        xlab (ty.Optional[str], optional): User-provided overriding
            label for x-axis. Defaults to None.
        ylab (ty.Optional[str], optional): User-provided overriding
        label for y-axis. Defaults to None.

    Returns:
        None
    """

    axes_obj.spines['left'].set_visible(True)
    axes_obj.spines['bottom'].set_visible(True)

    if xlab is not None:
        axes_obj.set_xlabel(xlab)

    else:
        axes_obj.set_xlabel(x)

    if ylab is not None:
        axes_obj.set_ylabel(ylab)

    else:
        axes_obj.set_xlabel(y)

    # blue
    for i in range(len(df.index)):
        if df.at[i, "within_hpd"] == 1.0:
            line_color = "darkturquoise"
            axes_obj.vlines(
                float(df.at[i, x]),
                float(df.at[i, "lower_95hpd"]),
                float(df.at[i, "higher_95hpd"]),
                line_color,
                linewidth=3.0,
                alpha=0.25)
            axes_obj.plot(
                float(df.at[i, x]),
                float(df.at[i, y]),
                marker="o",
                markersize=3,
                color=line_color)

    # place red on top of blue
    for i in range(len(df.index)):
        if df.at[i, "within_hpd"] == 0.0:
            line_color = "orangered"
            axes_obj.vlines(float(df.at[i, x]), float(df.at[i, "lower_95hpd"]), float(df.at[i, "higher_95hpd"]), line_color, alpha=0.8)
            axes_obj.plot(float(df.at[i, x]), float(df.at[i, y]), marker="o", markersize=3, color=line_color)

    fig_obj.canvas.draw()


if __name__ == "__main__":
    # testing below

    # note that pjgui uses matplotlib.figure.Figure
    # (which is part of Matplotlib's OOP class library)
    # here, we instead use pyplot's figure, which is the
    # Matlab-like state-machine API
    fig = plt.figure()
    axes = fig.add_axes([0.25, 0.2, 0.5, 0.6])
    axes.patch.set_alpha(0.0)
    axes.spines['left'].set_visible(True)
    axes.spines['bottom'].set_visible(True)
    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)

    ###########
    # Violins #
    ###########
    df1 = pd.DataFrame({"Total taxon count": [1, 3, 4, 5, 5, 5, 5, 6, 7, 10],
                        "Program": ["PJ" for i in range(10)]})

    df2 = pd.DataFrame({"Total taxon count": [1, 3, 4, 5, 5, 5, 5, 6, 7, 10],
                        "Program": ["Other" for i in range(10)]})

    df = pd.concat([df1, df2], axis=0)
    print(tabulate(df, headers=df.head(), tablefmt="pretty"))

    # uncomment if testing violins
    # plot_violins(fig, axes, df, "Program", "Total taxon count",
    #              xlab="Program", ylab="Total taxon count")
    # plt.show()

    #################
    # Coverage bars #
    #################

    # initializing model
    model_fp = "examples/coverage_files/r_b_exp.pj"
    dag_obj = cmdp.script2dag(model_fp, in_pj_file=True)

    # reading true values
    scalar_constant_df = pd.DataFrame({"r_b": dag_obj.value("r_b")})
    # print(tabulate(scalar_constant_df, headers=scalar_constant_df.head(),
    #                tablefmt="pretty", showindex=False))

    # reading hpd table
    hpd_df = pjread.read_csv_into_dataframe(
        "./examples/coverage_files/r_b_exp.csv")
    # print(tabulate(hpd_df, headers=hpd_df.head(), tablefmt="pretty"))

    full_cov_df = pd.concat([scalar_constant_df, hpd_df], axis=1)
    full_cov_df = pjorg.add_within_hpd_col(full_cov_df, "r_b")
    print(
        tabulate(full_cov_df, headers=full_cov_df.head(),
                 tablefmt="pretty", showindex=False).lstrip()
    )

    # uncomment if testing coverage lines
    plot_intervals(fig, axes, full_cov_df, "r_b",
                   "posterior_mean", ylab="Posterior mean")

    plt.show()
