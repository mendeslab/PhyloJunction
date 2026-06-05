import typing as ty
from abc import ABC, abstractmethod

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class RegionStatus():
    _comm_class_id: int

    def __init__(self) -> None:
        pass

    @property
    def comm_class_id(self) -> int:
        return self._comm_class_id

    @comm_class_id.setter
    def comm_class_id(self, id: int) -> None:
        self._comm_class_id = id



class EvolRelevantEvent(ABC):
    """Evolution-relevant event.

    Examples of evolution-relevant events are those represented by
    stochastic maps (e.g., range contraction, range expansion,
    extinction, dispersal) or the (paleogeographic) appearance and
    disappearance of barriers.

    Parameters:
        n_chars (int): Number of characters involved in event.
        age (float): Age of event (age of present moment is 0.0).
        time (float, optional): Time of event (0.0 at origin or
            root). Defaults to None.
        char_status_dict (dict): Dictionary with integer as keys
            representing the position of a character in a bit
            pattern (e.g., if there are three regions A, B and C,
            0 = A, 1 = B, 2 = C), and a RegionStatus as values. This
            RegionStatus represents some type of status the region is
            in (e.g., if the character represents a region being
            part of a biogeographic range, then 'DISCONNECTED' could
            mean that that region is disconnected from others relative
            to a range-splitting event upon speciation).
    """

    _n_chars: int
    _age: ty.Optional[float]
    _time: ty.Optional[float]
    _char_status_dict: ty.Dict[int, bool]
    _short_str: str

    @abstractmethod
    def __init__(self,
                 n_chars: int,
                 age: ty.Optional[float],
                 time: ty.Optional[float] = None) -> None:

        self._n_chars = n_chars
        self._age = None
        if age is not None:
            self._age = float(age)

        self._time = None
        if time is not None:
            self._time = float(time)

        self._char_status_dict = dict()
        for idx in range(n_chars):
            self._char_status_dict[idx] = False

    @abstractmethod
    def short_str(self) -> str:
        pass

    @property
    def age(self) -> float:
        return self._age

    @property
    def time(self) -> float:
        return self._time

    @property
    def char_status_dict(self):
        return self._char_status_dict

    @char_status_dict.setter
    def char_status_dict(self,
                         char_pos_in_bit_patt: int,
                         char_status: RegionStatus) -> \
            ty.Dict[int, RegionStatus]:

        self._char_status_dict[char_pos_in_bit_patt] = char_status

    def char_status(self, char_pos_in_bit_patt) -> bool:
        return self._char_status_dict[char_pos_in_bit_patt]

    def __lt__(self, other) -> bool:
        return self._age > other._age
