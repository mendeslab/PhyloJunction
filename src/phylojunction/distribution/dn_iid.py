"""Explicit IID vectors, including vectors of distributions with vector-valued draws."""

from phylojunction.pgm.pgm import DistrForSampling


class DnIID(DistrForSampling):
    """One draw is a vector of n independent draws from distribution."""

    DN_NAME = "IID"

    # Vector length belongs to this distribution; model repetition does not.
    def __init__(self, n, distribution):
        self.validate_draw_count(n)
        if not isinstance(distribution, DistrForSampling):
            raise TypeError("iid requires a distribution object.")
        self.n = n
        self.distribution = distribution

    def draw(self):
        return self.distribution.draw_many(self.n)
