import sys
# sys.path.extend(["../../", "../phylojunction"]) # necessary to run it as standalone on command line (from its own path all the way up to phylojunction/)
import typing as ty
import math

# pj imports
import phylojunction.pgm.pgm as pgm

def get_rev_str_from_dn_parametric_obj(dn_obj: pgm.DistrForSampling) -> ty.Tuple[int, int, ty.List[str]]:
    
    rev_str_list: ty.List[str]
    n_sim: int = dn_obj.n_samples
    n_repl: int = dn_obj.n_repl
    
    rev_str_list = dn_obj.get_rev_inference_spec_info()
    
    return n_sim, n_repl, rev_str_list

#########################################
# Functions for converting parametric   #
# distribution PJ call into .Rev syntax #
#########################################

def get_unif_rev_inference_spec_info(n_samples: int, min_param_list: ty.List[float], max_param_list: ty.List[float], parent_node_tracker: ty.Optional[ty.Dict[str, str]]) -> ty.List[str]:
    
    rev_str_list: ty.List[str] = [] # return

    min_list: ty.List[float] = min_param_list
    max_list: ty.List[float] = max_param_list

    # so mypy won't complain
    if isinstance(parent_node_tracker, dict):
        for ith_sim in range(n_samples):
            rev_str_list.append("dnUniform(lower=" + str(min_list[ith_sim]) + ", upper=" + str(max_list[ith_sim]) + ")")

    return rev_str_list


def get_exponential_rev_inference_spec_info(n_samples: int, exp_scale_or_rate_list: ty.List[float], exp_rate_parameterization: bool, parent_node_tracker: ty.Optional[ty.Dict[str, str]]) -> ty.List[str]:
    
    rev_str_list: ty.List[str] = [] # return

    scale_or_rate_list: ty.List[float] = exp_scale_or_rate_list
    rate_parameterization: bool = exp_rate_parameterization

    # so mypy won't complain
    if isinstance(scale_or_rate_list, list) and isinstance(parent_node_tracker, dict):
        
        # if user wants scale parameterization in PJ, we need to invert it
        # for use with rev
        if not rate_parameterization:
            scale_or_rate_list = [ 1.0/float(v) for v in scale_or_rate_list ]

        for ith_sim in range(n_samples):
            ith_sim_str = "dnExponential(lambda="

            # if we can find a node that holds the value of the rate, we use it
            try:
                # key: arg in PJ syntax, value: NodeDAG name passed as arg
                ith_sim_str += parent_node_tracker["rate"]

            except:
                ith_sim_str += str(scale_or_rate_list[ith_sim])
            ith_sim_str += ")"

            rev_str_list.append(ith_sim_str)

    return rev_str_list


def get_gamma_rev_inference_spec_info(n_samples: int, gamma_shape_param_list: ty.List[float], gamma_scale_or_rate_param_list: ty.List[float], gamma_rate_parameterization: bool, parent_node_tracker: ty.Optional[ty.Dict[str, str]]) -> ty.List[str]:
    
    rev_str_list: ty.List[str] = [] # return

    shape_list: ty.List[float] = gamma_shape_param_list
    scale_or_rate_list: ty.List[float] = gamma_scale_or_rate_param_list
    rate_parameterization: bool = gamma_rate_parameterization
    
    # so mypy won't complain
    if isinstance(parent_node_tracker, dict):
        # if user wants scale parameterization in PJ, we need to invert it
        # for use with rev (rev uses rate parameterization)
        if not rate_parameterization:
            scale_or_rate_list = [ 1.0/float(v) for v in scale_or_rate_list ]

        for ith_sim in range(n_samples):
            ith_sim_str = "dnGamma("

            # if we can find a node that holds the value of the shape
            # parameter, we use it
            try:
                # returns NodeDAG, and we grab its name
                ith_sim_str += parent_node_tracker["shape"]
            
            except:
                ith_sim_str += str(shape_list[ith_sim])

            ith_sim_str += ", rate="

            # if we can find a node that holds the value of the scale parameter, we use it
            try:
                # returns NodeDAG, and we grab its name
                ith_sim_str += parent_node_tracker["scale"]

            except:
                ith_sim_str += str(scale_or_rate_list[ith_sim])
            ith_sim_str += ")"

            rev_str_list.append(ith_sim_str)

    return rev_str_list

def get_normal_rev_inference_spec_info(n_samples: int, norm_mean_param_list: ty.List[float], norm_sd_param_list: ty.List[float], parent_node_tracker: ty.Optional[ty.Dict[str, str]]) -> ty.List[str]:
    
    rev_str_list: ty.List[str] = [] # return

    real_mean_list: ty.List[float] = norm_mean_param_list
    real_sd_list: ty.List[float] = norm_sd_param_list

    # so mypy won't complain
    if isinstance(parent_node_tracker, dict):
        for ith_sim in range(n_samples):
            ith_sim_str = "dnNormal(mean="

            # if we can find a node that holds the value of the mean, we use it
            try:
                ith_sim_str += parent_node_tracker["mean"] # returns NodeDAG, and we grab its name
            except:
                ith_sim_str += str(real_mean_list[ith_sim])

            ith_sim_str += ", sd="

            # if we can find a node that holds the value of the sd, we use it
            try:
                ith_sim_str += parent_node_tracker["sd"] # returns NodeDAG, and we grab its name
            except:
                ith_sim_str += str(real_sd_list[ith_sim])
            ith_sim_str += ")"

            rev_str_list.append(ith_sim_str)

    return rev_str_list


# Translate lognormal parameters while retaining dependencies on named model nodes.
def get_ln_rev_inference_spec_info(
        n_samples: int,
        ln_mean_list: ty.List[float],
        ln_sd_list: ty.List[float],
        parent_node_tracker: ty.Optional[ty.Dict[str, str]],
        log_space: bool = True) -> ty.List[str]:
    
    rev_str_list: ty.List[str] = [] # return
    
    real_mean_list: ty.List[float] = ln_mean_list
    real_sd_list: ty.List[float] = ln_sd_list

    # so mypy won't complain
    if isinstance(parent_node_tracker, dict):
        # Rev uses the underlying normal parameters. In natural-space mode PJ's
        # location is the median exp(mu), so only that parameter needs a logarithm.
        if not log_space:
            real_mean_list = [math.log(float(value)) for value in real_mean_list]
        
        # real_mean_list and real_sd_list will have n_sim values inside, even if they are all the same
        for ith_sim in range(n_samples):
            ith_sim_str = "dnLognormal(mean="

            # if we can find a node that holds the value of the mean, we use it
            try:
                # returns NodePGM, and we grab its name
                mean_expr = parent_node_tracker["meanlog"]
                ith_sim_str += mean_expr if log_space else "ln(" + mean_expr + ")"

            except:
                ith_sim_str += str(real_mean_list[ith_sim])

            ith_sim_str += ", sd="

            # if we can find a node that holds the value of the sd, we use it
            try:
                # returns NodeDAG, and we grab its name
                ith_sim_str += parent_node_tracker["sdlog"]

            except:
                ith_sim_str += str(real_sd_list[ith_sim])
            ith_sim_str += ")"

            rev_str_list.append(ith_sim_str)

    return rev_str_list