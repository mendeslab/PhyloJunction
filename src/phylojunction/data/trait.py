"""Describe the primary trait stored on an annotated tree's nodes."""

from dataclasses import dataclass


@dataclass(frozen=True)
class DiscreteTrait:
    """A finite state space stored in the named node attribute."""

    states: int
    # Name of the node attribute holding the trait value.
    name: str = "state"


@dataclass(frozen=True)
class ContinuousTrait:
    """A real-valued trait stored in the named node attribute."""

    # Name of the node attribute holding the trait value.
    name: str = "trait"
