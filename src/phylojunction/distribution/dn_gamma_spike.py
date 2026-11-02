"""A draw attaches gamma-spike rates to a copy of one tree's reconstructed edges."""
import copy
from phylojunction.pgm.pgm import DistrForSampling
from phylojunction.data.tree import AnnotatedTree
from phylojunction.distribution.dn_parametric import finite, boolean
from phylojunction.calculation.gamma_spike_clock_model import populate_gamma_spike_mol_rates_on_rec_tr_branches


class DnGammaSpikeClock(DistrForSampling):
    """Relative rates are supplied values; each draw samples the punctuated component."""
    DN_NAME = 'GammaSpikeClock'

    # Validate without RNG calls or changes to the caller's tree. A private reconstruction
    # supplies the vector length and is copied for every draw, including draws through IID.
    def __init__(self, tree, attr_name, clock_rate, spike_mean, spike_shape,
                 br_specific_rel_rate, strict_clock_when_no_stubs=False):
        if not isinstance(tree, AnnotatedTree):
            raise TypeError('gamma_spike_clock requires one tree.')
        if not isinstance(attr_name, str) or not attr_name.isidentifier():
            raise ValueError('attr_name must be a nonempty Python identifier.')
        self.clock_rate = finite(clock_rate, 'clock_rate')
        self.spike_mean = finite(spike_mean, 'spike_mean')
        self.spike_shape = finite(spike_shape, 'spike_shape')
        if self.clock_rate <= 0 or self.spike_mean < 0 or self.spike_shape <= 0:
            raise ValueError('clock_rate and spike_shape must be positive; spike_mean must be nonnegative.')
        if not isinstance(br_specific_rel_rate, list):
            raise TypeError('br_specific_rel_rate must be a vector.')
        self.rates = [finite(rate, 'relative rate') for rate in br_specific_rel_rate]
        if any(rate < 0 for rate in self.rates):
            raise ValueError('Relative rates must be nonnegative.')
        self.strict = boolean(strict_clock_when_no_stubs)
        self.attr_name = attr_name
        self.tree = copy.deepcopy(tree)
        if len(self.rates) != self.tree.n_rec_branches():
            raise ValueError('Provide one relative rate per reconstructed edge (n_rec_branches).')

    # The stored tree is a template, never the returned value. Preserve times, traits and
    # maps; only the copy's reconstructed edges receive rates and serializable annotations.
    def draw(self):
        result = copy.deepcopy(self.tree)
        populate_gamma_spike_mol_rates_on_rec_tr_branches(
            result, self.clock_rate, self.spike_mean, self.spike_shape, self.rates, self.strict)
        for node in result.tree_reconstructed:
            setattr(node, self.attr_name, node.mol_rate)
            node.annotations.add_bound_attribute(self.attr_name)
        result.clock_model_attr_name = self.attr_name
        return result
