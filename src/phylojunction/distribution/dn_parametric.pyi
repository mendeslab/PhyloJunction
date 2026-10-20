import math
import random
import typing as ty
import numpy as np
from scipy.stats import expon, lognorm, norm, gamma, uniform, beta
from phylojunction.pgm.pgm import DistrForSampling

def finite(value, name):
    ...

def boolean(value):
    ...

class DnNormal(DistrForSampling):
    DN_NAME = 'Normal'

    def __init__(self, mean, sd):
        ...

    def draw(self):
        ...

    def draw_many(self, count):
        ...

    @staticmethod
    def draw_normal(n_samples: int, mean_param: float, sd_param: float) -> ty.Union[np.float64, np.ndarray]:
        ...

class DnLogNormal(DistrForSampling):
    DN_NAME = 'Log-normal'

    def __init__(self, meanlog, sdlog, log_space=True):
        ...

    def draw(self):
        ...

    def draw_many(self, count):
        ...

    @staticmethod
    def draw_ln(n_samples: int, mean_param: float, sd_param: float, scale: float=1.0, log_space: bool=True) -> ty.Union[np.float64, np.ndarray]:
        ...

class DnExponential(DistrForSampling):
    DN_NAME = 'Exponential'

    def __init__(self, rate, rate_parameterization=True):
        ...

    def draw(self):
        ...

    def draw_many(self, count):
        ...

    @staticmethod
    def draw_exp(n_samples: int, scale_or_rate_param: float, rate_parameterization: bool=True) -> np.ndarray:
        ...

class DnGamma(DistrForSampling):
    DN_NAME = 'Gamma'

    def __init__(self, shape, scale, rate_parameterization=False):
        ...

    def draw(self):
        ...

    def draw_many(self, count):
        ...

    @staticmethod
    def draw_gamma(n_samples: int, shape_param: float, scale_or_rate_param: float, rate_parameterization: bool=False) -> ty.Union[np.float64, np.ndarray]:
        ...

class DnUnif(DistrForSampling):
    DN_NAME = 'Uniform'

    def __init__(self, min, max):
        ...

    def draw(self):
        ...

    def draw_many(self, count):
        ...

    @staticmethod
    def draw_unif(n_samples: int, min_param: float, max_param: float) -> ty.Union[np.float64, np.ndarray]:
        ...

class DnBeta(DistrForSampling):
    DN_NAME = 'Beta'

    def __init__(self, alpha, beta):
        ...

    def draw(self):
        ...

    def draw_many(self, count):
        ...

class DnCauchy(DistrForSampling):
    DN_NAME = 'Cauchy'
    def __init__(self, location, scale): ...
    def draw(self): ...
    def draw_many(self, count): ...

class DnCategorical(DistrForSampling):
    DN_NAME = 'Categorical'

    def __init__(self, cats, probs):
        ...

    def draw(self):
        ...

    def draw_many(self, count):
        ...

    @staticmethod
    def draw_categorical(n_samples: int, cats_param: ty.List[int], probs_param: ty.List[float]) -> ty.List[int]:
        ...
