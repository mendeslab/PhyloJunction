from phylojunction.data.trait import DiscreteTrait, ContinuousTrait
import os
import typing as ty
from phylojunction.pgm.execution import ExecutionResults
import dendropy as dp
import numpy as np
import pandas
import pandas as pd
import matplotlib.pyplot as plt
import statistics as stat
import pickle
from phylojunction.data.tree import AnnotatedTree
from phylojunction.data.attribute_transition import AttributeTransition
import phylojunction.pgm.pgm as pgm
import phylojunction.interface.cmdbox.cmd_parse as cmdp
from tabulate import tabulate

def write_str_list(outfile_handle: ty.IO, content_string_list: ty.List[str]) -> None:
    ...

def write_fig_to_file(outfile_path: str, fig_obj: plt.Figure) -> None:
    ...

def write_data_df(outfile_handle: ty.IO, data_df: pd.DataFrame, format='csv') -> None:
    ...

def initialize_scalar_dataframe(sample_size: int, n_repl: int=1, summaries_avg_over_repl: bool=False) -> pd.DataFrame:
    ...

def initialize_tree_dataframe(sample_size: int, n_repl: int=1, summaries: bool=False, summaries_avg_over_repl: bool=False) -> pd.DataFrame:
    ...

def prep_trees_rb_smap_dfs(value_records, tree_dag_node_name_list: ty.List[str], mapped_attr_name: str) -> ty.Tuple[ty.Dict[str, ty.List[pd.DataFrame]], ty.Dict[str, ty.List[str]]]:
    ...

def prep_data_df(value_records, write_nex_states: bool=False, bit_format: bool=False) -> ty.Tuple[ty.List[ty.Union[pd.DataFrame, ty.Dict[int, pd.DataFrame]]], ty.List[ty.Dict[str, pd.DataFrame]]]:
    ...

def prep_data_filepaths_dfs(scalar_output_stash: ty.List[ty.Union[pd.DataFrame, ty.Dict[int, pd.DataFrame]]], tree_output_stash: ty.List[ty.Union[ty.Dict[str, pd.DataFrame], ty.Dict[str, str]]]) -> ty.Tuple[ty.List[str], ty.List[ty.Union[pd.DataFrame, str]]]:
    ...

def dump_pgm_data(dir_string: str, dag_obj: 'ExecutionResults', prefix: str='', write_nex_states: bool=False, bit_format: bool=False) -> None:
    ...

def dump_trees_rb_smap_dfs(dir_string: str, dag_obj: 'ExecutionResults', tr_dag_node_name_list: ty.List[str], mapped_attr_name: str, prefix: str='') -> None:
    ...

def dump_serialized_pgm(file_name: str, dag_obj: 'ExecutionResults', cmd_log_list: ty.List[str], prefix: str='', to_folder: bool=False) -> None:
    ...

def get_write_inference_rev_scripts(all_sims_model_spec_list: ty.List[str], all_sims_mcmc_logging_spec_list: ty.List[str], dir_list: ty.List[str], prefix: str='', write2file: bool=False) -> ty.List[str]:
    ...
