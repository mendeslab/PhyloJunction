import sys
import time
import typing as ty
import numpy as np  # type: ignore
import pandas as pd
from collections import defaultdict

# pj imports
# import phylojunction.pgm.pgm as pgm
import phylojunction.utility.exception_classes as ec


def print_progress(idx: int, iterator_len: int) -> None:
    """Print progress bar on iteration for longer tests

    Args:
        idx (int): Index to keep track of where in iteration one is
        iterator_len (int): Size of the iterator whose progress is
            keeping track of
    """

    j = (idx + 1) / iterator_len
    sys.stdout.write('\r')
    sys.stdout.write("[%-20s] %d%%" % ('=' * int(20 * j), 100 * j))
    sys.stdout.flush()
    time.sleep(0.25)


def autovivify(levels=1, final=dict) -> ty.DefaultDict:
    return (defaultdict(final) if levels < 2 else
            defaultdict(lambda: autovivify(levels - 1, final)))


def create_str_defaultdict() -> defaultdict:
    return defaultdict(str)


def check_and_vectorize_if_must(
        param_list:
        ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]],
        dn_name: str,
        size_to_grow: int = 1) -> ty.List[ty.List[ty.Union[int, float, str]]]:
    """Check number of provided values and vectorize if necessary.

    Args:
        param_list (tuple): Parameter values (could be scalar
            integers, floats, strings, or any of those in a list).
        dn_name (str): Name of the distribution whose parameter values
            we are checking.
        size_to_grow (int): Number of values we need, which is the
            specified number of simulations ('n_samples' in
            DistributionDAG). Defaults to 1.

    Returns:
        (list): Two-dimensional list of potentially vectorized
            parameter values. First dimension are samples
            (simulations), second dimension are replicates.
    """

    # using numpy array as a hack to get nested-ness of param_list
    param_list_array: ty.Any = np.array(param_list, dtype=object)

    # each entry in shape is a dimension
    n_params = len(param_list_array.shape)

    vectorizable_param_list = param_list
    vectorized_param_list: ty.List[ty.List[ty.Union[int, float, str]]] = \
        list()

    one_scalar_provided: bool = False
    one_scalar_in_list_provided: bool = False
    one_par_right_value_count: bool = False

    ##################################
    # First we prepare the container #
    ##################################

    # scalar provided by itself, e.g., rate = 1.0
    if isinstance(vectorizable_param_list, (int, float, str)):
        vectorizable_param_list = [vectorizable_param_list]
        one_scalar_provided = True

    # dealing with weird cases -- must set some flags
    # and/or initialize variables appropriately
    elif isinstance(param_list, list):
        # single scalar provided by itself inside list,
        # e.g, param_list = [1.0]
        if len(param_list) == 1:
            one_scalar_in_list_provided = True

        # we have a single parameter inside a list of size_to_grow already,
        # but no nestedness; this happens when we do
        # rv ~ someDn(n=5, param=[1, 2, 3, 4, 5])
        #
        # here, we want to just paste those values as-is, inside
        # vectorized_param_list, but we do need to initialize this second
        # dimension already so the nestedness is correct
        elif len(param_list) == size_to_grow and n_params == 1:
            one_par_right_value_count = True
            vectorized_param_list.append([])

    ######################
    # Now we populate it #
    ######################

    # each v here will be a different parameter (e.g., mean and sd)
    # in case there is more than one parameter
    if isinstance(vectorizable_param_list, list):
        for v in vectorizable_param_list:
            # a single value was provided as scalar, either by itself, or
            # inside a list by itself
            if isinstance(v, (int, float, str)):
                # param list is a single scalar, e.g., mean=1.0, vectorize!
                if one_scalar_provided:
                    vectorized_param_list.append(
                        [v for i in range(size_to_grow)])

                # param list is a single scalar, by itself, inside a list,
                # e.g., mean=[[1.0]], vectorize!
                elif one_scalar_in_list_provided:
                    vectorized_param_list.append([v] * size_to_grow)

                # param list contains a single parameter, and we already have
                # the right number of values
                elif one_par_right_value_count:
                    vectorized_param_list[0].append(v)

            elif isinstance(v, list):
                n_val = len(v)

                # more values than specified number of samples,
                # can't vectorize!
                if n_val > size_to_grow:
                    raise ec.DimensionalityError(dn_name)

                # don't know how to multiply if more than one element,
                # can't vectorize!
                elif n_val > 1 and n_val < size_to_grow:
                    raise ec.DimensionalityError(dn_name)

                # a single value was provided as list, vectorize!
                elif n_val == 1:
                    vectorized_param_list.append(
                        [v[0] for i in range(size_to_grow)])

                # more than 1 value, but same as specified number of samples,
                # no need to vectorize!
                else:
                    vectorized_param_list.append(v)

    return vectorized_param_list


def get_ellapsed_time_in_seconds(start: float, end: float) -> int:
    """Calculate ellapsed time in seconds.

    Args:
        start (float): Start of time window
        end (float): End of time window

    Returns:
        int: Ellapsed time in seconds in time window
    """
    # returns (min, sec)
    ellapsed_minutes, ellapsed_secs = divmod(end - start, 60)

    return int(ellapsed_minutes * 60 + ellapsed_secs) # in seconds
    # return int(ellapsed_minutes)


def is_val_in_interval(
        val: ty.Union[int, float, np.float64],
        lower: ty.Union[int, float, np.float64],
        upper: ty.Union[int, float, np.float64]) -> bool:
    """Return True/False if numerical value is in (lower, upper].

    Args:
        val (ty.Union[int, float, np.float64]t): Numerical value
        lower (ty.Union[int, float, np.float64]): Lower end of interval
        upper (ty.Union[int, float, np.float64]): Upper end of interval

    Returns:
        bool: Whether val is in (lower, upper]
    """

    if lower < val and val <= upper:
        return True
    else:
        return False


def get_covg(full_cov_df: pd.DataFrame, par_name: str) -> float:
    """Return coverage from 'within_hpd' column in DataFrame.
    
    Args:
        full_cov_df (pandas.DataFrame): DataFrame object containing
            column named 'within_hpd' with 0's or 1's for parameters
            outside and inside an HPD, respectively.
        par_name (str): Name of the focal parameter for error
            printing.
            
    Returns:
        float: Bayesian coverage for the parameter
    """

    if "within_hpd" not in list(full_cov_df):
        raise ec.MissingColumnName("within_hpd",
                                   "Could not compute coverage for " + \
                                    par_name)

    int_list_in_phd = full_cov_df.loc[:, 'within_hpd']

    return \
        float(sum(i for i in int_list_in_phd)) \
            / float(int_list_in_phd.size)


def symmetric_difference(
        set1: ty.Set[ty.Any],
        set2: ty.Set[ty.Any]) -> ty.Set[ty.Any]:
    """Return symmetric difference among two sets"""
    result = set1

    for elem in set2:
        try:
            result.remove(elem)

        except KeyError:
            result.add(elem)

    return result

##############################################################################


if __name__ == "__main__":
    # can be called from utility/
    # $ python3 helper_functions.py
    #
    # can also be called from phylojunction/
    # $ python3 utility/helper_functions.py
    # or
    # $ python3 -m utility.helper_functions
    #
    # can also be called from VS Code, if open folder is phylojuction/

    d = autovivify(levels=3, final=dict)
    d[0][1][2] = 3
    print(d)

    print_progress(49, 100)

    print(get_ellapsed_time_in_seconds(0.0, 125))  # 2 minutes (and 5 seconds)
