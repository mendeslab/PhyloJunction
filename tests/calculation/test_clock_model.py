import unittest

# pj imports
import phylojunction.calculation.clock_model as pjclock

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class TestClockModel(unittest.TestCase):
    """Test the 'ClockModel' abstract base class' own contract.

    Concrete clock models (e.g., a Gamma-spike clock model) are
    tested on their own, once implemented. This test file only
    guards the shape of the abstraction itself.
    """

    def test_cannot_instantiate_abstract_class(self) -> None:
        """'ClockModel' must remain abstract (no concrete model yet)."""

        with self.assertRaises(TypeError):
            pjclock.ClockModel()  # type: ignore[abstract]


if __name__ == "__main__":
    unittest.main()
