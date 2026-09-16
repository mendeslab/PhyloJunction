"""Continuous-trait rate functions shared by QuaSSE simulators."""

import numpy as np
from scipy.special import expit


class LogisticRate:
    """The four-parameter logistic link used by the Java SSE package."""

    # Store validated scalar parameters; each DAG sample gets its own rate object.
    def __init__(self, y0, y1, midpoint, slope):
        self.y0, self.y1, self.midpoint, self.slope = map(float, (y0, y1, midpoint, slope))
        if not np.all(np.isfinite([self.y0, self.y1, self.midpoint, self.slope])):
            raise ValueError("Logistic parameters must be finite.")
        if min(self.y0, self.y1) < 0:
            raise ValueError("Logistic rate plateaus must be nonnegative.")

    # The logistic weight interpolates the plateaus, including decreasing curves.
    # expit handles either tail without overflowing an exponential.
    def __call__(self, x):
        x = np.asarray(x, dtype=float)
        if self.slope == 0 or self.y0 == self.y1:
            return np.full_like(x, self.y0 if self.y0 == self.y1 else self.y0 / 2 + self.y1 / 2)
        with np.errstate(over="ignore"):
            z = self.slope * (x - self.midpoint)
        # expit(-z) equals 1 - expit(z), but retains the small positive tail when
        # expit(z) rounds to one. Evaluate both weights directly to avoid cancellation.
        return self.y0 * expit(-z) + self.y1 * expit(z)

    def __str__(self):
        return f"logistic(y0={self.y0}, y1={self.y1}, midpoint={self.midpoint}, slope={self.slope})"
