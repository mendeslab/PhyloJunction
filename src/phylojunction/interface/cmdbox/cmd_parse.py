"""Compile PhyloJunction scripts and execute complete models.

Compilation records dependencies without drawing values. Model executions are
explicit; iid describes vectors within one execution. Interactive additions use
ExecutionResults.append and publish the returned result only after success.
"""
from pathlib import Path
from phylojunction.pgm.pgm import DirectedAcyclicGraph

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


# Compile expressions without sampling; callers choose the number of whole-model runs.
def compile_model(source):
    model = DirectedAcyclicGraph()
    for line in source.splitlines():
        model.add_definition(line)
    return model


# Read a script (or source text), compile its shared definitions, and execute it.
def script2dag(script_file_path_or_model_spec, in_pj_file=True, random_seed=None, runs=1):
    source = Path(script_file_path_or_model_spec).read_text() if in_pj_file else script_file_path_or_model_spec
    return compile_model(source).execute(runs, random_seed)
