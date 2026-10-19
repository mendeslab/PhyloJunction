"""Registered continuous-trait rate constructors."""

from phylojunction.calculation.continuous_sse import (
    ConstantRate, LogisticRate, GaussianRate, StepRate, LinearRate,
    SkewGaussianRate, QuadraticRate, AbsoluteRate,
)


# Keep constructor dispatch and the grammar's allowed parameters in one place.
RATE_CONSTRUCTORS = {
    "quasse_constant": (ConstantRate, ("rate",)),
    "quasse_logistic": (LogisticRate, ("y0", "y1", "midpoint", "slope")),
    "quasse_gaussian": (GaussianRate, ("baseline", "center_rate", "center", "width")),
    "quasse_step": (StepRate, ("left", "right", "threshold")),
    "quasse_linear": (LinearRate, ("intercept", "slope")),
    "quasse_skew_gaussian": (SkewGaussianRate, ("baseline", "amplitude", "location", "width", "skew")),
    "quasse_quadratic": (QuadraticRate, ("baseline", "strength", "center")),
    "quasse_absolute": (AbsoluteRate, ("baseline", "strength", "center")),
}
