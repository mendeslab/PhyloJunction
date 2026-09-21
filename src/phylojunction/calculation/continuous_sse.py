"""Continuous-trait rate functions shared by QuaSSE simulators."""

import numpy as np
from scipy.special import expit, ndtr


class RateFunction:
    """A nonnegative trait-dependent rate, optionally capped explicitly."""

    # Validate scalar parameters once, sharing cap and sign rules across rate families.
    def _initialize(self, parameters, nonnegative=(), positive=(), cap=None):
        for name, value in parameters.items():
            value = float(value)
            if not np.isfinite(value):
                raise ValueError(f"{name} must be finite.")
            if (name in nonnegative and value < 0) or (name in positive and value <= 0):
                raise ValueError(f"Invalid {name}: {value}.")
            setattr(self, name, value)
        self.cap = None if cap is None else float(cap)
        if self.cap is not None and (not np.isfinite(self.cap) or self.cap < 0):
            raise ValueError("cap must be finite and nonnegative.")
        self._parameter_names = tuple(parameters)

    # Apply the explicit cap after evaluation; cap=0 avoids needless overflow entirely.
    def __call__(self, x):
        x = np.asarray(x, dtype=float)
        if self.cap == 0:
            return np.zeros_like(x)
        with np.errstate(over="ignore", under="ignore"):
            result = self._evaluate(x)
        return result if self.cap is None else np.minimum(self.cap, result)

    @property
    def is_identically_zero(self):
        return self.upper_bound == 0

    # A cap is itself a bound; None means no finite representable bound is available.
    @property
    def upper_bound(self):
        bound = self._bound()
        if self.cap is not None:
            return min(bound, self.cap)
        return bound if np.isfinite(bound) else None

    # Include the cap only when supplied, keeping the original logistic description.
    def __str__(self):
        parameters = [f"{name}={getattr(self, name)}" for name in self._parameter_names]
        if self.cap is not None:
            parameters.append(f"cap={self.cap}")
        return f"{self._label}({', '.join(parameters)})"


class LogisticRate(RateFunction):
    """The four-parameter logistic link used by the Java SSE package."""

    # Store validated scalar parameters; each DAG sample gets its own rate object.
    _label = "logistic"

    def __init__(self, y0, y1, midpoint, slope, cap=None):
        self._initialize(dict(y0=y0, y1=y1, midpoint=midpoint, slope=slope), ("y0", "y1"), cap=cap)

    def _bound(self):
        return max(self.y0, self.y1)

    # The logistic weight interpolates the plateaus, including decreasing curves.
    # expit handles either tail without overflowing an exponential.
    def _evaluate(self, x):
        x = np.asarray(x, dtype=float)
        if self.slope == 0 or self.y0 == self.y1:
            return np.full_like(x, self.y0 if self.y0 == self.y1 else self.y0 / 2 + self.y1 / 2)
        with np.errstate(over="ignore"):
            z = self.slope * (x - self.midpoint)
        # expit(-z) equals 1 - expit(z), but retains the small positive tail when
        # expit(z) rounds to one. Evaluate both weights directly to avoid cancellation.
        return self.y0 * expit(-z) + self.y1 * expit(z)


class ConstantRate(RateFunction):
    """A rate constant in both trait and time."""
    _label = "constant"

    def __init__(self, rate, cap=None):
        self._initialize(dict(rate=rate), ("rate",), cap=cap)

    def _bound(self):
        return self.rate

    def _evaluate(self, x):
        return np.full_like(x, self.rate)


# Rescale before subtracting only when subtraction overflowed, preserving wide profiles.
def _standardized(x, center, width):
    delta = x - center
    # np.where evaluates both arms, including unused inf-inf expressions.
    with np.errstate(invalid="ignore"):
        return np.where(np.isinf(delta), x / width - center / width, delta / width)


class GaussianRate(RateFunction):
    """A Gaussian peak or trough between nonnegative center and tail rates."""
    _label = "gaussian"

    def __init__(self, baseline, center_rate, center, width, cap=None):
        self._initialize(dict(baseline=baseline, center_rate=center_rate, center=center, width=width),
                         ("baseline", "center_rate"), ("width",), cap)

    def _bound(self):
        return max(self.baseline, self.center_rate)

    # Convex weights keep peaks and troughs nonnegative, including zero center rates.
    def _evaluate(self, x):
        z = _standardized(x, self.center, self.width)
        exponent = -0.5 * z * z
        return self.baseline * -np.expm1(exponent) + self.center_rate * np.exp(exponent)


class StepRate(RateFunction):
    """Two constant rates; the threshold belongs to the right interval."""
    _label = "step"

    def __init__(self, left, right, threshold, cap=None):
        self._initialize(dict(left=left, right=right, threshold=threshold), ("left", "right"), cap=cap)

    def _bound(self):
        return max(self.left, self.right)

    def _evaluate(self, x):
        return np.where(x < self.threshold, self.left, self.right)


class LinearRate(RateFunction):
    """A linear rate with negative values replaced by zero, not an implicit upper cap."""
    _label = "linear"

    def __init__(self, intercept, slope, cap=None):
        self._initialize(dict(intercept=intercept, slope=slope), cap=cap)

    def _bound(self):
        return float("inf") if self.slope else max(0., self.intercept)

    def _evaluate(self, x):
        return np.maximum(0., self.intercept + self.slope * x)


class SkewGaussianRate(RateFunction):
    """An asymmetric peak; location and amplitude need not be its maximum coordinates."""
    _label = "skew_gaussian"

    def __init__(self, baseline, amplitude, location, width, skew, cap=None):
        self._initialize(dict(baseline=baseline, amplitude=amplitude, location=location, width=width, skew=skew),
                         ("baseline", "amplitude"), ("width",), cap)

    def _bound(self):
        return self.baseline + 2 * self.amplitude

    # Twice the Gaussian times a normal CDF gives the symmetric peak when skew=0.
    def _evaluate(self, x):
        z = _standardized(x, self.location, self.width)
        weight = np.exp(-0.5 * z * z)
        if self.skew:
            weight *= 2 * ndtr(self.skew * z)
        return self.baseline + self.amplitude * weight


# Scale distance before squaring, avoiding overflow of (x-center)^2 for tiny strengths.
# The fallback subtraction handles opposite-sign traits near the floating-point limit.
def _distance_growth(x, center, strength, power):
    if strength == 0:
        return np.zeros_like(x)
    scale = np.sqrt(strength) if power == 2 else strength
    delta = x - center
    with np.errstate(invalid="ignore"):
        scaled = np.where(np.isinf(delta), x * scale - center * scale, delta * scale)
    return scaled * scaled if power == 2 else np.abs(scaled)


class QuadraticRate(RateFunction):
    """A nonnegative parabola, unbounded in trait unless explicitly capped."""
    _label = "quadratic"

    def __init__(self, baseline, strength, center, cap=None):
        self._initialize(dict(baseline=baseline, strength=strength, center=center),
                         ("baseline", "strength"), cap=cap)

    def _bound(self):
        return float("inf") if self.strength else self.baseline

    def _evaluate(self, x):
        return self.baseline + _distance_growth(x, self.center, self.strength, 2)


class AbsoluteRate(RateFunction):
    """A V-shaped rate, unbounded in trait unless explicitly capped."""
    _label = "absolute"

    def __init__(self, baseline, strength, center, cap=None):
        self._initialize(dict(baseline=baseline, strength=strength, center=center),
                         ("baseline", "strength"), cap=cap)

    def _bound(self):
        return float("inf") if self.strength else self.baseline

    def _evaluate(self, x):
        return self.baseline + _distance_growth(x, self.center, self.strength, 1)
