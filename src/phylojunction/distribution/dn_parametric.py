import sys
import typing as ty
import math
import numpy as np
import random
from scipy.stats import expon, lognorm, norm, gamma, uniform  # type: ignore

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.calculation.math_utils as mu
import phylojunction.utility.exception_classes as ec
import phylojunction.utility.helper_functions as pjh
import phylojunction.inference.revbayes.rb_dn_parametric as rbpar

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class DnLogNormal(pgm.DistrForSampling):

    DN_NAME = "Log-normal"

    n_samples: int
    n_repl: int
    ln_mean_list: ty.List[float]
    ln_sd_list: ty.List[float]
    ln_log_space: bool
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.Union[bool, ty.List[ty.Union[int, float, str]]]]
    parent_node_tracker: ty.Optional[ty.Dict[str, str]]

    @staticmethod
    def draw_ln(
            n_samples: int,
            mean_param: float,
            sd_param: float,
            scale: float = 1.0,
            log_space: bool = True) -> ty.Union[np.float64, np.ndarray]:
        """Return sample from log-normal distribution.

        Args:
            n_samples (int): Number of draws (sample size).
            mean_param (float): Mean (scale) of log-normal.
            sd_param (float): Std. deviation (shape) of log-normal
                distribution.
            scale (float): Mean (scale) of log-normal distribution in
                log-space. Defaults to 1.0.
            log_space (bool, optional): If true, mean_param is the underlying normal
                mean; otherwise it is the lognormal median/scale exp(mean).
                sd_param is always the underlying normal standard deviation.

        Returns:
            (list): List of floats sampled from log-normal distribution.
        """

        if not log_space:
            # have to use scale if in natural space
            return lognorm.rvs(s=sd_param,
                               scale=mean_param,
                               loc=0.0,
                               size=n_samples)

        # s is the shape parameter (std dev)
        # scale is the mean 
        #
        # scipy assumes the scale is in natural space, so if the user passed
        # it in log-scale (default), we need to exponentiate it
        return lognorm.rvs(s=sd_param,
                           scale=math.exp(mean_param),
                           loc=0.0,
                           size=n_samples)

    # validation of pars happens in dn_grammar
    def __init__(
            self,
            n_samples: int,
            n_repl: int,
            ln_mean: ty.List[float],
            ln_sd: ty.List[float],
            ln_log_space: bool,
            parent_node_tracker: ty.Optional[ty.Dict[str, str]] = None) \
                -> None:

        self.param_dict = dict()
        self.n_samples = n_samples
        self.n_repl = n_repl
        self.ln_mean_arg = ln_mean
        self.ln_sd_arg = ln_sd
        self.ln_log_space = ln_log_space

        # one element per parameter
        check_sample_size_return = \
            self.init_check_vectorize_sample_size(
                [self.ln_mean_arg,
                 self.ln_sd_arg])

        if isinstance(check_sample_size_return, list):
            self.vectorized_params = check_sample_size_return

        # params
        self.ln_mean_list = ty.cast(ty.List[float], self.vectorized_params[0])
        self.ln_sd_list = ty.cast(ty.List[float], self.vectorized_params[1])

        # for inference, we need to keep track of parent node names
        if isinstance(parent_node_tracker, pgm.DeterministicNodeDAG):
            raise ec.ObjInitInvalidArgError(
                self.DN_NAME,
                ("One of the arguments is a deterministic node. "
                 "This is not allowed. Ignoring command..."))

        self.parent_node_tracker = parent_node_tracker

    def generate(self) -> ty.List[float]:
        sampled_values: ty.List[float] = list()

        try:
            if len(self.ln_mean_list) == 1 and len(self.ln_sd_list) == 1:
                n_draws: int = self.n_samples * self.n_repl

                # so mypy won't complain
                return ty.cast(
                    ty.List[float],
                    DnLogNormal.draw_ln(
                        n_draws,
                        self.ln_mean_list[0],
                        self.ln_sd_list[0],
                        log_space=self.ln_log_space).tolist())

            else:
                # for each simulation, draw n_repl times, using
                # one mean and sd at a time
                for i in range(self.n_samples):
                    # so mypy won't complain
                    repl = ty.cast(
                        ty.List,
                        DnLogNormal.draw_ln(
                            self.n_repl,
                            self.ln_mean_list[i],
                            self.ln_sd_list[i],
                            log_space=self.ln_log_space).tolist())
                    sampled_values += repl

                return sampled_values

        except Exception as e:
            print("An error occurred: ", type(e).__name__, " - ", e)

            raise ec.GenerateFailError(
                self.DN_NAME,
                "A \'mean_param\' and \'sd_param\' must be specified.")

    def init_check_vectorize_sample_size(
            self,
            param_list: ty.List[ty.Any] = []) \
            -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]:
        
        return pjh.check_and_vectorize_if_must(param_list,
                                               self.DN_NAME,
                                               size_to_grow=self.n_samples)

    def get_rev_inference_spec_info(self) -> ty.List[str]:

        return rbpar.get_ln_rev_inference_spec_info(
            self.n_samples,
            self.ln_mean_list,
            self.ln_sd_list,
            self.parent_node_tracker,
            log_space=self.ln_log_space)


class DnNormal(pgm.DistrForSampling):

    DN_NAME = "Normal"

    n_samples: int
    n_repl: int
    norm_mean_param_list: ty.List[float]
    norm_sd_param_list: ty.List[float]
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.Union[bool, ty.List[ty.Union[int, float, str]]]]
    parent_node_tracker: ty.Optional[ty.Dict[str, str]]

    @staticmethod
    def draw_normal(
            n_samples: int,
            mean_param: float,
            sd_param: float) -> ty.Union[np.float64, np.ndarray]:
        """Return sample from normal distribution.

        Args:
            n_samples (int): Number of draws (sample size).
            mean_param (float): Mean (location) of normal distribution.
            sd_param (float): Std. deviation (scale) of normal distribution.

        Returns:
            (list): List of values sampled from normal distribution.
        """

        # debugging
        # print("normal: mean = " + str(mean_param) + " sd = " + str(sd_param))

        return norm.rvs(scale=sd_param, loc=mean_param, size=n_samples)

    # validation of pars happens in dn_grammar
    def __init__(
            self,
            n_samples: int,
            n_repl: int,
            norm_mean_param: ty.List[float],
            norm_sd_param: ty.List[float],
            parent_node_tracker: ty.Optional[ty.Dict[str, str]]) -> None:

        self.param_dict = dict()

        # so mypy won't complain...
        # if isinstance(pars[0], int) and isinstance(pars[1], int):
        #     self.n_samples = pars[0]
        #     self.n_repl = pars[1]
        self.n_samples = n_samples
        self.n_repl = n_repl
        self.norm_mean_param_arg = norm_mean_param
        self.norm_sd_param_arg = norm_sd_param

        # makes sure these parameters are lists, or converts them into list
        # also multiplying element (when only one is provided) if necessary
        # check_sample_size_return = self.check_sample_size(pars[2:4])

        # one element per parameter
        check_sample_size_return = \
            self.init_check_vectorize_sample_size(
                [self.norm_mean_param_arg,
                 self.norm_sd_param_arg])
        if isinstance(check_sample_size_return, list):
            self.vectorized_params = check_sample_size_return

        # params
        # self.param_dict["mean_param"] = self.vectorized_params[0]
        # self.param_dict["sd_param"] = self.vectorized_params[1]
        self.norm_mean_param_list = \
            ty.cast(ty.List[float], self.vectorized_params[0])
        self.norm_sd_param_list = \
            ty.cast(ty.List[float], self.vectorized_params[1])

        # for inference, we need to keep track of parent node names
        self.parent_node_tracker = parent_node_tracker

    def generate(self) -> ty.List[float]:
        sampled_values: ty.List[float] = []

        try:
            if len(self.norm_mean_param_list) == 1 and \
                    len(self.norm_sd_param_list) == 1:
                n_draws: int = self.n_samples * self.n_repl

                # so mypy won't complain
                return ty.cast(ty.List,
                               DnNormal.draw_normal(
                                   n_draws,
                                   self.norm_mean_param_list[0],
                                   self.norm_sd_param_list[0]).tolist())

            else:
                # for each simulation, draw n_repl times, using
                # one mean and sd at a time
                for i in range(self.n_samples):
                    # so mypy won't complain
                    # for some reason, just checking type like below won't work here
                    # ... no clue why (I'm still calling isinstance() just to be safe)
                    # mypy is only satisfied if I forcefully cast into float
                    repl = ty.cast(
                        ty.List,
                        DnNormal.draw_normal(
                            self.n_repl,
                            self.norm_mean_param_list[i],
                            self.norm_sd_param_list[i]).tolist())
                    sampled_values += repl

                return sampled_values

        except Exception as e:
            print("An error occurred: ", type(e).__name__, " - ", e)

            raise ec.GenerateFailError(
                self.DN_NAME,
                ("A \'norm_mean_param\' and \'norm_sd_param\' "
                 "must be specified."))

    def init_check_vectorize_sample_size(
            self,
            param_list: ty.List[ty.Any] = []) \
            -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]:

        return pjh.check_and_vectorize_if_must(
            param_list,
            self.DN_NAME,
            size_to_grow=self.n_samples)

    def get_rev_inference_spec_info(self) -> ty.List[str]:

        return rbpar.get_normal_rev_inference_spec_info(
            self.n_samples,
            self.norm_mean_param_list,
            self.norm_sd_param_list,
            self.parent_node_tracker)


class DnExponential(pgm.DistrForSampling):

    DN_NAME = "Exponential"

    n_samples: int
    n_repl: int
    exp_scale_or_rate_list: ty.List[float]
    exp_rate_parameterization: bool
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: \
    ty.Dict[str, ty.Union[bool, ty.List[ty.Union[int, float, str]]]]
    parent_node_tracker: ty.Optional[ty.Dict[str, str]]

    @staticmethod
    def draw_exp(
            n_samples: int,
            scale_or_rate_param: float,
            rate_parameterization: bool = True) \
                -> np.ndarray:
        """Return sample from exponential distribution.

        Args:
            n_samples (int): Number of draws (sample size).
            scale_or_rate_param (float): Scale (default) or rate of
                exponential distribution.
            rate_parameterization (bool, optional): Argument of
                'scale_or_rate_param' is rate instead of scale.
                Defaults to 'True'.

        Returns:
            (list): List of floats sampled from exonential distribution.
        """

        if rate_parameterization:
            scale_or_rate_param = 1.0 / scale_or_rate_param

        return expon.rvs(loc=0.0, scale=scale_or_rate_param, size=n_samples)

    # validation of pars happens in dn_grammar
    def __init__(
            self,
            n_samples: int,
            n_repl: int,
            scale_or_rate_param: ty.List[float],
            rate_parameterization: bool,
            parent_node_tracker: ty.Optional[ty.Dict[str, str]] = None) \
                -> None:

        self.param_dict = dict()
        self.n_samples = n_samples
        self.n_repl = n_repl
        self.exp_scale_or_rate_param_arg = scale_or_rate_param
        self.exp_rate_parameterization = rate_parameterization

        check_sample_size_return = \
            self.init_check_vectorize_sample_size(
                [self.exp_scale_or_rate_param_arg])

        if isinstance(check_sample_size_return, list):
            self.vectorized_params = check_sample_size_return

        # params
        # self.param_dict["scale_or_rate_param"] = self.vectorized_params[0]
        self.exp_scale_or_rate_param_list = \
            ty.cast(ty.List[float], self.vectorized_params[0])

        # for inference, we need to keep track of parent node names
        self.parent_node_tracker = parent_node_tracker

    def generate(self) -> ty.List[float]:
        sampled_values: ty.List[float] = []

        try:
            if len(self.exp_scale_or_rate_param_list) == 1:
                n_draws: int = self.n_samples * self.n_repl

                # so mypy won't complain
                return ty.cast(
                    ty.List[float],
                    DnExponential.draw_exp(
                        n_draws,
                        self.exp_scale_or_rate_param_list[0],
                        rate_parameterization=self.exp_rate_parameterization).tolist())

            else:
                for i in range(len(self.exp_scale_or_rate_param_list)):
                    # so mypy won't complain
                    repl = ty.cast(
                        ty.List[float],
                        DnExponential.draw_exp(
                            self.n_repl,
                            self.exp_scale_or_rate_param_list[i],
                            rate_parameterization=self.exp_rate_parameterization).tolist())
                    sampled_values += repl

                return sampled_values

        except Exception as e:
            print("An error occurred: ", type(e).__name__, " - ", e)

            raise ec.GenerateFailError(
                self.DN_NAME,
                ("A \'scale_or_rate_param\' and \'rate_parameterization\' "
                 "must be specified."))

    def init_check_vectorize_sample_size(self, param_list: ty.List[ty.Any] = []) \
            -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]:

        return pjh.check_and_vectorize_if_must(
            param_list,
            self.DN_NAME,
            size_to_grow=self.n_samples)

    def get_rev_inference_spec_info(self) -> ty.List[str]:

        return rbpar.get_exponential_rev_inference_spec_info(
            self.n_samples,
            self.exp_scale_or_rate_param_list,
            self.exp_rate_parameterization,
            self.parent_node_tracker)


class DnGamma(pgm.DistrForSampling):

    DN_NAME = "Gamma"

    n_samples: int
    n_repl: int
    gamma_shape_param_list: ty.List[float]
    gamma_scale_or_rate_param_list: ty.List[float]
    gamma_rate_parameterization: bool
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.Union[bool, ty.List[ty.Union[int, float, str]]]] = dict()
    parent_node_tracker: ty.Optional[ty.Dict[str, str]]

    @staticmethod
    def draw_gamma(
        n_samples: int,
        shape_param: float,
        scale_or_rate_param: float,
        rate_parameterization: bool = False) \
            -> ty.Union[np.float64, np.ndarray]:
        """Return sample from gamma distribution.

        Args:
            n_samples (int): Number of draws (sample size).
            shape_param (float): Gamma distribution shape parameter
                (represented by alpha or kappa sometimes).
            scale_or_rate_param (float): Gamma distribution scale or
                rate parameter.
            rate_parameterization (bool, optional): Argument of
                'scale_or_rate_param' is rate instead of scale.
                Defaults to 'False'.

        Returns:
            (list): List of floats sampled from exonential distribution.
        """

        if rate_parameterization:
            scale_or_rate_param = 1.0 / scale_or_rate_param

        return gamma.rvs(
            a=shape_param,
            scale=scale_or_rate_param,
            size=n_samples)

    # validation of pars happens in dn_grammar
    def __init__(
            self,
            n_samples: int,
            n_repl: int,
            shape_param: ty.List[float],
            scale_or_rate_param: ty.List[float],
            rate_parameterization: bool,
            parent_node_tracker: ty.Optional[ty.Dict[str, str]] = None) \
                -> None:

        self.param_dict = dict()
        self.n_samples = n_samples
        self.n_repl = n_repl
        self.gamma_rate_parameterization = rate_parameterization
        self.gamma_shape_param_arg = shape_param
        self.gamma_scale_or_rate_param_arg = scale_or_rate_param

        # one element per parameter
        check_sample_size_return = \
            self.init_check_vectorize_sample_size(
                [self.gamma_shape_param_arg,
                 self.gamma_scale_or_rate_param_arg])
        if isinstance(check_sample_size_return, list):
            self.vectorized_params = check_sample_size_return

        # params
        self.gamma_shape_param_list = \
            ty.cast(ty.List[float], self.vectorized_params[0])
        self.gamma_scale_or_rate_param_list = \
            ty.cast(ty.List[float], self.vectorized_params[1])

        # for inference, we need to keep track of parent node names
        self.parent_node_tracker = parent_node_tracker

    def generate(self) -> ty.List[float]:
        sampled_values: ty.List[float] = list()

        try:
            if len(self.gamma_shape_param_list) == 1 and \
                    len(self.gamma_scale_or_rate_param_list) == 1:
                n_draws: int = self.n_samples * self.n_repl

                # so mypy won't complain
                return ty.cast(
                    ty.List[float],
                    DnGamma.draw_gamma(
                        n_draws,
                        self.gamma_shape_param_list[0],
                        self.gamma_scale_or_rate_param_list[0],
                        rate_parameterization = \
                            self.gamma_rate_parameterization).tolist())

            else:
                for i in range(self.n_samples):
                    # so mypy won't complain
                    repl = ty.cast(
                        ty.List[float],
                        DnGamma.draw_gamma(
                            self.n_repl,
                            self.gamma_shape_param_list[i],
                            self.gamma_scale_or_rate_param_list[i],
                            rate_parameterization = \
                                self.gamma_rate_parameterization).tolist())

                    sampled_values += repl

                return sampled_values

        except Exception as e:
            print("An error occurred: ", type(e).__name__, " - ", e)

            raise ec.GenerateFailError(
                self.DN_NAME,
                ("A \'shape_param\' and \'scale_or_rate_param\' "
                 "must be specified."))

    def init_check_vectorize_sample_size(
        self,
        param_list: ty.List[ty.Any] = []) \
            -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]:
        
        return pjh.check_and_vectorize_if_must(
            param_list,
            self.DN_NAME,
            size_to_grow=self.n_samples)

    def get_rev_inference_spec_info(self) -> ty.List[str]:

        return rbpar.get_gamma_rev_inference_spec_info(
            self.n_samples,
            self.gamma_shape_param_list,
            self.gamma_scale_or_rate_param_list,
            self.gamma_rate_parameterization,
            self.parent_node_tracker)


class DnUnif(pgm.DistrForSampling):

    DN_NAME = "Uniform"

    n_samples: int
    n_repl: int
    min_param_list: ty.List[float]
    max_param_list: ty.List[float]
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.Union[bool, ty.List[ty.Union[int, float, str]]]]
    parent_node_tracker: ty.Optional[ty.Dict[str, str]]

    @staticmethod
    def draw_unif(n_samples: int, min_param: float, max_param: float) \
        -> ty.Union[np.float64, np.ndarray]:
        """Return sample from uniform distribution.

        Args:
            n_samples (int): Number of draws (sample size).
            min_param (float): Minimum value.
            max_param (float): Maximum value.

        Returns:
            (list): List of floats sampled from exonential distribution.
        """

        # see definition of uniform from scipy (need to subtract min_param)
        # if we want the "usual" behavior of a uniform distribution, where
        # scale behaves as a maximum
        return uniform.rvs(loc=min_param,
                           scale=max_param - min_param,
                           size=n_samples)

    # validation of pars happens in dn_grammar
    def __init__(
            self,
            n_samples: int,
            n_repl: int,
            min_param: ty.List[float],
            max_param: ty.List[float],
            parent_node_tracker: ty.Optional[ty.Dict[str, str]] = None) \
                -> None:

        self.param_dict = dict()
        self.n_samples = n_samples
        self.n_repl = n_repl
        self.min_param_arg = min_param
        self.max_param_arg = max_param

        # one element per parameter
        check_sample_size_return = \
            self.init_check_vectorize_sample_size(
                [self.min_param_arg,
                 self.max_param_arg])

        if isinstance(check_sample_size_return, list):
            self.vectorized_params = check_sample_size_return

        # params
        self.min_param_list = \
            ty.cast(ty.List[float], self.vectorized_params[0])
        self.max_param_list = \
            ty.cast(ty.List[float], self.vectorized_params[1])

        # for inference, we need to keep track of parent node names
        self.parent_node_tracker = parent_node_tracker

    def generate(self) -> ty.List[float]:
        sampled_values: ty.List[ty.Union[float, ty.Any]] = []

        try:
            if len(self.min_param_list) == 1 and \
                    len(self.max_param_list) == 1:
                n_draws: int = self.n_samples * self.n_repl

                # so mypy won't complain
                return ty.cast(
                    ty.List[float],
                    DnUnif.draw_unif(n_draws,
                                     self.min_param_list[0],
                                     self.max_param_list[0]).tolist())

            else:
                # for i in range(len(min_params)):
                for i in range(len(self.min_param_list)):
                    # so mypy won't complain
                    repl = ty.cast(
                        ty.List[float],
                        DnUnif.draw_unif(self.n_repl,
                                         self.min_param_list[i],
                                         self.max_param_list[i]).tolist())
                    sampled_values += repl

                return sampled_values

        except Exception as e:
            print("An error occurred: ", type(e).__name__, " - ", e)

            raise ec.GenerateFailError(
                self.DN_NAME,
                ("A \'min_param\' and \'max_param\' must be specified. "
                 "Exiting..."))

    def init_check_vectorize_sample_size(
        self,
        param_list: ty.List[ty.Any] = []) \
            -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]:

        return pjh.check_and_vectorize_if_must(
            param_list,
            self.DN_NAME,
            size_to_grow=self.n_samples)

    def get_rev_inference_spec_info(self) -> ty.List[str]:

        return rbpar.get_unif_rev_inference_spec_info(
            self.n_samples,
            self.min_param_list,
            self.max_param_list,
            self.parent_node_tracker)


class DnCategorical(pgm.DistrForSampling):

    DN_NAME = "Categorical"

    n_samples: int
    n_repl: int
    cats_param_list: ty.List[int]
    probs_param_list: ty.List[float]
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.Union[bool, ty.List[ty.Union[int, float, str]]]]
    parent_node_tracker: ty.Optional[ty.Dict[str, str]]

    @staticmethod
    def draw_categorical(n_samples: int,
                         cats_param: ty.List[int],
                         probs_param: ty.List[float]) \
            -> ty.List[int]:
        """Return sample from categorical distribution.

        Args:
            n_samples (int): Number of draws (sample size).
            cats_param (int): List of categories (e.g., states) that one
                wants to sample. Must be the same size as 'probs_param'.
            probs_param (float): List of category drawing probabilities
                (i.e., weights) to use when drawing categories. Must be
                the same size as 'cats_param'.

        Returns:
            (list): List of integers provided/representing categories
                being sampled.
        """

        return random.choices(cats_param, weights=probs_param, k=n_samples)

    def __init__(
            self,
            n_samples: int,
            n_repl: int,
            cats_param: ty.List[int],
            probs_param: ty.List[float],
            parent_node_tracker: ty.Optional[ty.Dict[str, str]] = None) \
                -> None:

        self.param_dict = dict()
        self.n_samples = n_samples
        self.n_repl = n_repl
        self.cats_param_arg = cats_param
        self.probs_param_arg = probs_param

        # one element per parameter
        # check_sample_size_return = \
        #     self.init_check_vectorize_sample_size(
        #         [self.cats_param_arg,
        #          self.probs_param_arg])
        #
        # if isinstance(check_sample_size_return, list):
        #     self.vectorized_params = check_sample_size_return

        # params
        self.cats_param_list = self.cats_param_arg  # 2D list
        self.probs_param_list = self.probs_param_arg  # 2D list
        # self.cats_param_list = self.vectorized_params[0]  # 2D list
        # self.probs_param_list = self.vectorized_params[1]  # 2D list

        # for inference, we need to keep track of parent node names
        self.parent_node_tracker = parent_node_tracker

    def init_check_vectorize_sample_size(
        self,
        param_list: ty.List[ty.Any] = []) \
            -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]:

        return pjh.check_and_vectorize_if_must(
            param_list,
            self.DN_NAME,
            size_to_grow=self.n_samples)

    def generate(self) -> ty.List[float]:
        sampled_values: ty.List[int] = []

        try:
            if len(self.cats_param_list) == 1 and \
                    len(self.probs_param_list) == 1:
                n_draws: int = self.n_samples * self.n_repl

                print("self.probs_param_list", self.probs_param_list)
                # so mypy won't complain
                # ty.cast(
                sample = DnCategorical.draw_categorical(n_draws,
                                                        self.cats_param_list[0],
                                                        self.probs_param_list[0])

                print("categorical sample", sample)
                return sample

            else:
                # for i in range(len(self.cats_param_list)):
                    # so mypy won't complain
                repl = DnCategorical.draw_categorical(
                    self.n_repl,
                    self.cats_param_list[0],
                    self.probs_param_list[0])#.tolist()

                return repl
                # sampled_values += repl

                # return sampled_values

        except Exception as e:
            print("An error occurred: ", type(e).__name__, " - ", e)

            raise ec.GenerateFailError(
                self.DN_NAME,
                ("A \'cats_param\' and \'probs_param\' must be specified. "
                 "Exiting..."))

    def get_rev_inference_spec_info(self) -> ty.List[str]:
        return ""


if __name__ == "__main__":
    # can be called from distribution/
    # $ python3 dn_parametric.py
    #
    # can also be called from phylojunction/
    # $ python3 distribution/dn_parametric.py
    # or
    # $ python3 -m distribution.dn_parametric
    #
    # can also be called from VS Code, if open folder is phylojuction/

    # ln1 = DnLogNormal([1, [10,9,8], [7,6,5], True]) # should throw exception
    n_samples = 3
    n_repl = 1

    dummy_parent_node_tracker: ty.Dict[str, str] = dict()

    ln2 = DnLogNormal(n_samples, n_repl, [-2, -2.1, -2.2], [0.1, 0.2, 0.15], True, dummy_parent_node_tracker)
    ln3 = DnLogNormal(n_samples, n_repl, [-2], [0.1], True, dummy_parent_node_tracker)
    ln4 = DnLogNormal(n_samples, n_repl, [-2], [0.1], True, dummy_parent_node_tracker)

    print(ln2.generate())
    print(ln3.generate())
    print(ln4.generate())
