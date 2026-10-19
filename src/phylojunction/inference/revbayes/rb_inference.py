"""RevBayes inference export for execution results."""
from phylojunction.inference.revbayes.execution import execution_to_rev


# Preserve the established entry point and directory return convention for both GUIs.
def dag_obj_to_rev_inference_spec(dag_obj, inference_root_dir, mcmc_chain_length=1000, prefix=''):
    """Convert model string to put in .Rev script.

    Args:
        dag_obj (ExecutionResults): Shared definitions and the values from
            complete executions of a script, including interactive additions.
        mcmc_chain_length (int): Number of MCMC iterations to carry
            out with RevBayes.
        inference_root_dir (str): Directory where scripts and results
            should be put, and from where RevBayes should be called.
        prefix (str): Prefix will be appended to directory and file
            names. Defaults to "".

    Returns:
        (tuple): Tuple containing three lists, with strings to be
            printed to files (2) and directory names (1).
    """
    
    return execution_to_rev(dag_obj, inference_root_dir, mcmc_chain_length, prefix)
