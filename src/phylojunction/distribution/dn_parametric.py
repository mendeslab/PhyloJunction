"""Scalar probability distributions; explicit IID and model execution own repetition."""
import math
import random
import typing as ty
import numpy as np
from scipy.stats import expon, lognorm, norm, gamma, uniform, beta
from phylojunction.pgm.pgm import DistrForSampling

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


# Distribution parameters must be finite scalars. Lists denote actual vector values,
# not a hidden batch of different parameter settings.
def finite(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite numeric scalar.")
    return float(value)


# Keep established quoted flags as well as boolean literals in the script language.
def boolean(value):
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.lower() in ('true', 'false', 't', 'f'):
        return value.lower() in ('true', 't')
    raise ValueError('Expected true or false.')


class DnNormal(DistrForSampling):
    """Normal distribution describing one scalar draw."""

    DN_NAME = 'Normal'

    # Validate and retain one parameter set, shared by all explicit IID draws.
    def __init__(self, mean, sd):
        self.mean, self.sd = finite(mean, 'mean'), finite(sd, 'sd')
        if self.sd <= 0:
            raise ValueError('Normal sd must be positive.')

    def draw(self):
        return self.draw_many(1)[0]

    def draw_many(self, count):
        self.validate_draw_count(count)
        return self.draw_normal(count, self.mean, self.sd).tolist()

    # Bulk sampler shared by scalar and IID draws; preserves the established parameterization.
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



class DnLogNormal(DistrForSampling):
    """Log-normal distribution describing one scalar draw."""

    DN_NAME = 'Log-normal'

    # Validate and retain one parameter set, shared by all explicit IID draws.
    def __init__(self, meanlog, sdlog, log_space=True):
        self.meanlog, self.sdlog = finite(meanlog, 'meanlog'), finite(sdlog, 'sdlog')
        self.log_space = boolean(log_space)
        if self.sdlog <= 0 or (not self.log_space and self.meanlog <= 0):
            raise ValueError('Lognormal sdlog and natural-space scale must be positive.')

    def draw(self):
        return self.draw_many(1)[0]

    def draw_many(self, count):
        self.validate_draw_count(count)
        return self.draw_ln(count, self.meanlog, self.sdlog, log_space=self.log_space).tolist()

    # Bulk sampler shared by scalar and IID draws; preserves the established parameterization.
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



class DnExponential(DistrForSampling):
    """Exponential distribution describing one scalar draw."""

    DN_NAME = 'Exponential'

    # Validate and retain one parameter set, shared by all explicit IID draws.
    def __init__(self, rate, rate_parameterization=True):
        self.rate = finite(rate, 'rate')
        self.rate_parameterization = boolean(rate_parameterization)
        if self.rate <= 0:
            raise ValueError('Exponential rate/scale must be positive.')

    def draw(self):
        return self.draw_many(1)[0]

    def draw_many(self, count):
        self.validate_draw_count(count)
        return self.draw_exp(count, self.rate, self.rate_parameterization).tolist()

    # Bulk sampler shared by scalar and IID draws; preserves the established parameterization.
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



class DnGamma(DistrForSampling):
    """Gamma distribution describing one scalar draw."""

    DN_NAME = 'Gamma'

    # Validate and retain one parameter set, shared by all explicit IID draws.
    def __init__(self, shape, scale, rate_parameterization=False):
        self.shape, self.scale = finite(shape, 'shape'), finite(scale, 'scale')
        self.rate_parameterization = boolean(rate_parameterization)
        if self.shape <= 0 or self.scale <= 0:
            raise ValueError('Gamma shape and scale/rate must be positive.')

    def draw(self):
        return self.draw_many(1)[0]

    def draw_many(self, count):
        self.validate_draw_count(count)
        return self.draw_gamma(count, self.shape, self.scale, self.rate_parameterization).tolist()

    # Bulk sampler shared by scalar and IID draws; preserves the established parameterization.
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



class DnUnif(DistrForSampling):
    """Uniform distribution describing one scalar draw."""

    DN_NAME = 'Uniform'

    # Validate and retain one parameter set, shared by all explicit IID draws.
    def __init__(self, min, max):
        self.min, self.max = finite(min, 'min'), finite(max, 'max')
        if self.min >= self.max:
            raise ValueError('Uniform min must be less than max.')

    def draw(self):
        return self.draw_many(1)[0]

    def draw_many(self, count):
        self.validate_draw_count(count)
        return self.draw_unif(count, self.min, self.max).tolist()

    # Bulk sampler shared by scalar and IID draws; preserves the established parameterization.
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



class DnBeta(DistrForSampling):
    """Beta distribution describing one scalar draw."""

    DN_NAME = 'Beta'

    # Validate and retain one parameter set, shared by all explicit IID draws.
    def __init__(self, alpha, beta):
        self.alpha, self.beta = finite(alpha, 'alpha'), finite(beta, 'beta')
        if self.alpha <= 0 or self.beta <= 0:
            raise ValueError('Beta shapes must be positive.')

    def draw(self):
        return self.draw_many(1)[0]

    def draw_many(self, count):
        self.validate_draw_count(count)
        return beta.rvs(self.alpha, self.beta, size=count).tolist()


class DnCategorical(DistrForSampling):
    """Categorical distribution describing one scalar draw."""

    DN_NAME = 'Categorical'

    # Validate and retain one parameter set, shared by all explicit IID draws.
    def __init__(self, cats, probs):
        if not isinstance(cats, list) or not isinstance(probs, list) or not cats or len(cats) != len(probs):
            raise ValueError('Categories and probabilities must be nonempty vectors of equal length.')
        self.cats = list(cats)
        self.probs = [finite(p, 'probability') for p in probs]
        if any(p < 0 for p in self.probs) or sum(self.probs) <= 0:
            raise ValueError('Category weights must be nonnegative with a positive sum.')

    def draw(self):
        return self.draw_many(1)[0]

    def draw_many(self, count):
        self.validate_draw_count(count)
        return self.draw_categorical(count, self.cats, self.probs)

    # Bulk sampler shared by scalar and IID draws; preserves the established parameterization.
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



if __name__ == "__main__":
    # python -m phylojunction.distribution.dn_parametric
    # Different parameter pairs describe different distributions; draw_many repeats
    # one distribution. Neither operation introduces a model-execution dimension.
    distributions = [DnLogNormal(mean, sd) for mean, sd in
                     zip([-2, -2.1, -2.2], [0.1, 0.2, 0.15])]
    print([distribution.draw() for distribution in distributions])
    distribution = DnLogNormal(-2, 0.1)
    print(distribution.draw_many(3))
    print(distribution.draw_many(3))
