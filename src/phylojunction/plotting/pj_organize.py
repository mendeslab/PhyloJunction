import matplotlib.pyplot as plt  # type: ignore
import seaborn as sns  # type: ignore
import pandas as pd  # type: ignore
import typing as ty
from tabulate import tabulate  # type: ignore

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.interface.cmdbox.cmd_parse as cmdp
import phylojunction.readwrite.pj_write as pjwrite
import phylojunction.readwrite.pj_read as pjread

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


def join_dataframes(
    pj_df: pd.DataFrame,
    compare_to_df: pd.DataFrame,
    thing_to_compare: str,
    summaries_avg_over_repl: bool = False) \
        -> pd.DataFrame:
    """_summary_

    Args:
        pj_df (pd.DataFrame): _description_
        compare_to_df (pd.DataFrame): _description_
        value_to_compare (str): _description_
        summaries_avg_over_repl (bool, optional): _description_.
            Defaults to False.

    Returns:
        pd.DataFrame: _description_
    """

    # debugging
    # print("thing_to_compare = " + thing_to_compare)
    # print(tabulate(compare_to_df,
    #                compare_to_df.head(),
    #                tablefmt="plain",
    #                showindex=False).lstrip())
    # print(tabulate(pj_df, pj_df.head(),
    #                tablefmt="plain",
    #                showindex=False).lstrip())

    sub_pj_df = pd.DataFrame()
    sub_compare_to_df = pd.DataFrame()

    # TODO:
    if summaries_avg_over_repl:
        pass

    #     'pj_df' will be:
    #         program, sample, replicate, thing_to_compare
    #     (make sure value to compare is float)
    #
    #     `compare_to_df` will be:
    #         program, sample, value_to_compare
    #     (make sure value to compare is float)

    else:
        type_dict = {
            "program": str,
            "sample": int,
            "replicate": int,
            locals()["thing_to_compare"]: float
        }

        try:
            sub_pj_df = \
                pj_df.loc[:,
                          ["program",
                           "sample",
                           "replicate",
                           thing_to_compare]
                          ]
            sub_pj_df = sub_pj_df.astype(type_dict)

            # debugging
            print(
                tabulate(sub_pj_df,
                         sub_pj_df.head(),
                         tablefmt="pretty",
                         showindex=False
                         ).lstrip()
            )

        # TODO: (PEP8) This is a bare except that needs to be made more explicit
        # later. I need to see what kinds of errors I may bump into later.
        except:
            print(("Could not build joint dataframe from PJ in order to "
                   "compare " + thing_to_compare + ". For now, ignoring. "
                   "Later implement Exception."))
            return pd.DataFrame()

        try:
            sub_compare_to_df = \
                compare_to_df.loc[:,
                                  ["program",
                                   "sample",
                                   "replicate",
                                   thing_to_compare]
                                  ]
            sub_compare_to_df = sub_compare_to_df.astype(type_dict)

            # debugging
            # print(
            #     tabulate(
            #         sub_compare_to_df,
            #         sub_compare_to_df.head(),
            #         tablefmt="pretty",
            #         showindex=False
            #         ).lstrip()
            # )

        # TODO: (PEP8) This is a bare except that needs to be made more explicit
        # later. I need to see what kinds of errors I may bump into later.
        except:
            print(("Could not build sub-dataframe from other simulator data "
                   "in order to compare " + thing_to_compare + ". For now,"
                   " ignoring. Later implement Exception."))
            return pd.DataFrame()

    joint_df = pd.concat([sub_pj_df, sub_compare_to_df], axis=0)

    # debugging
    # print(
    #     tabulate(joint_df,
    #              joint_df.head(),
    #              tablefmt="pretty",
    #              showindex=False
    #              ).lstrip()
    # )

    return joint_df


def add_within_hpd_col(
        df: pd.DataFrame,
        val_col_name: str) -> pd.DataFrame:
    """Add column (1 or 0) if in or out of HPD

    Args:
        df (pd.DataFrame): DataFrame
        val_col_name (str): Name of column containing value

    Returns:
        pd.DataFrame: Input dataframe with new column
    """

    if df.empty:
        print(("Could not add a coverage column because dataframe is empty."
               " For now, ignoring. Later implement Exception"))

    else:
        for i in range(len(df.index)):
            truth = float(df.at[i, val_col_name])
            lower = float(df.at[i, "lower_95hpd"])
            upper = float(df.at[i, "higher_95hpd"])

            if (lower < truth) and (truth <= upper):
                df.at[i, "within_hpd"] = 1.0

            else:
                df.at[i, "within_hpd"] = 0.0

    return df


if __name__ == "__main__":
    # initializing figure
    #
    # note that pjgui uses matplotlib.figure.Figure
    # (which is part of Matplotlib's OOP class library)
    # here, we instead use pyplot's figure, which is the
    # Matlab-like state-machine API
    comparison_fig = plt.figure()

    comparison_fig_axes = comparison_fig.add_axes([0.25, 0.2, 0.5, 0.6])
    comparison_fig_axes.patch.set_alpha(0.0)
    comparison_fig_axes.xaxis.set_ticks([])
    # comparison_fig_axes.yaxis.set_ticks([])
    comparison_fig_axes.spines['left'].set_visible(True)
    comparison_fig_axes.spines['bottom'].set_visible(True)
    comparison_fig_axes.spines['right'].set_visible(False)
    comparison_fig_axes.spines['top'].set_visible(False)

    #################################
    # Testing parsing functions for #
    # comparison tab                #
    #################################

    # initializing model
    model_fp = "examples/multiple_scalar_tree_plated.pj"
    dag_obj = cmdp.script2dag(model_fp, in_pj_file=True)

    # getting pj dataframes
    scalar_output_stash: \
        ty.List[ty.Union[pd.DataFrame, ty.Dict[int, pd.DataFrame]]]
    tree_output_stash: ty.List[ty.Dict[str, pd.DataFrame]]
    output_df_list: ty.List[pd.DataFrame]
    scalar_output_stash, tree_output_stash = \
        pjwrite.prep_data_df(dag_obj)  # still prefixless
    _, output_df_list = \
        pjwrite.prep_data_filepaths_dfs(
            scalar_output_stash,
            tree_output_stash)

    # output_df_list[0] # df with sample, replicate, rv1, rv2
    # output_df_list[1] # df with sample, replicate, rv3, rv4
    # output_df_list[2] # df with sample, replicate, rv5
    # output_df_list[3] # df with sample, replicate, rv6
    # output_df_list[4] # df with sample, summary, rv5, rv6
    # output_df_list[5] # df with sample, replicate, trs1 complete
    # output_df_list[6] # df with sample, replicate, trs2 complete
    # output_df_list[7] # df with sample, replicate, trs1 annotated complete
    # output_df_list[8] # df with sample, replicate, trs2 annotated complete
    # output_df_list[9] # df with sample, replicate, trs1 reconstructed
    # output_df_list[10] # df with sample, replicate, trs2 reconstructed

    # df with sample, replicate, trs1 annotated reconstructed
    # output_df_list[11]

    # df with sample, replicate, trs2 annotated reconstructed
    # output_df_list[12]

    # df with sample, replicate, origin_age, root_age, n_total, n_extant,
    # n_extinct, n_sa for trs1
    # output_df_list[13]

    # df with sample, replicate, origin_age, root_age, n_total, n_extant,
    # n_extinct, n_sa for trs2
    # output_df_list[14]

    # output_df_list[15] # df with average and std. dev for each sample
    # origin_age, root_age, n_total, n_extant, n_extinct, n_sa for trs2 (?)
    # print(output_df_list[16]) # df with internal node states
    #
    # TODO: keep checking later, there are more...

    # for df in output_df_list:
    #     df["program"] = "PJ"

    output_df_list[2]["program"] = "PJ"

    # getting dataframe to compare to
    csv_fp = "examples/compare_files/" \
        + "comparison_multiple_scalar_tree_plated_rv5.csv"
    # csv_fp = "examples/compare_files/" \
    #     + "comparison_multiple_scalar_tree_plated_trs2.csv"
    compare_to_df = pjread.read_csv_into_dataframe(csv_fp)

    # try to join dataframes
    # output_df_list[2] has data for "rv5",
    # output_df_list[8] has data for "trs2"
    joint_df = join_dataframes(output_df_list[2], compare_to_df,
                               thing_to_compare="rv5",
                               summaries_avg_over_repl=False)
    # joint_df = join_dataframes(compare_to_df, output_df_list[8],
    #                            value_to_compare="root_age",
    #                            summaries_avg_over_repl=False)

    #################################
    # Testing parsing functions for #
    # validation tab                #
    #################################

    # initializing model
    model_fp = "examples/coverage_files/r_b_exp1.pj"
    dag_obj = cmdp.script2dag(model_fp, in_pj_file=True)

    # reading true values
    scalar_output_stash, tree_output_stash = pjwrite.prep_data_df(dag_obj)
    scalar_constant_df = scalar_output_stash[0]

    # print(tabulate(scalar_constant_df, headers=scalar_constant_df.head(),
    #                tablefmt="pretty", showindex=False))

    # reading hpd table
    hpd_df = pjread.read_csv_into_dataframe(
        "./examples/coverage_files/r_b_exp1.csv")
    # print(tabulate(hpd_df, headers=hpd_df.head(), tablefmt="pretty"))

    full_cov_df = pd.concat([scalar_constant_df, hpd_df], axis=1)

    full_cov_df = add_within_hpd_col(full_cov_df, "r_b")

    # print(tabulate(full_cov_df, headers=full_cov_df.head(),
    #                tablefmt="pretty", showindex=False).lstrip())

    # print(full_cov_df)
