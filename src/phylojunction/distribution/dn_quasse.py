"""QuaSSE simulation with interchangeable numerical methods.

The initial method follows diversitree/R/simulate-quasse.R. Events precede trait
updates. The final age-limited interval and its event probability are shortened;
size stopping omits the terminating birth and its subsequent trait update.
"""
from phylojunction.data.trait import ContinuousTrait

import time
import numpy as np
import dendropy as dp

from phylojunction.calculation.continuous_sse import RateFunction
from phylojunction.data.tree import AnnotatedTree
from phylojunction.pgm.pgm import DistrForSampling
from phylojunction.utility import exception_classes as ec


# Reject fractional counts rather than silently truncating Python or grammar inputs.
def _integer(value, name, minimum=1):
    number = float(value)
    if not np.isfinite(number) or not number.is_integer() or number < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}.")
    return int(number)


class DnQuaSSE(DistrForSampling):
    DN_NAME = "DnQuaSSE"

    # Keep model vectors separate from scalar execution settings, then validate both.
    def __init__(self, birth_rate, death_rate, stop, stop_value, start_trait=0.0,
                 drift=0.0, diffusion=0.0, n=1, nr=1, method="diversitree", k=500,
                 runtime_limit=300, max_steps=1000000, max_alive=100000, rng_seed=None,
                 sampling_prob=1.0, cond_surv=True, cond_spn=False, cond_obs_both_sides=False,
                 min_rec_taxa=None, max_rec_taxa=None, max_n_attempts=200, dt_max=None, max_nodes=None):
        self.n_sim = _integer(n, "n")
        self.n_repl = _integer(nr, "nr")
        self.method = method
        self.stop = stop
        self.k = _integer(k, "k")
        self.max_steps = _integer(max_steps, "max_steps")
        self.max_alive = _integer(max_alive, "max_alive")
        self.max_nodes = None if max_nodes is None else _integer(max_nodes, "max_nodes")
        self.dt_max = None if dt_max is None else float(dt_max)
        if self.dt_max is not None and (not np.isfinite(self.dt_max) or self.dt_max <= 0):
            raise ValueError("dt_max must be positive and finite.")
        self.runtime_limit = float(runtime_limit)
        if not np.isfinite(self.runtime_limit) or self.runtime_limit <= 0:
            raise ValueError("runtime_limit must be positive and finite.")
        if method != "diversitree":
            raise ValueError(f"Unknown QuaSSE method {method!r}; supported: diversitree.")
        if stop not in ("age", "size"):
            raise ValueError('stop must be "age" or "size".')
        self.birth_rate, self.death_rate = birth_rate, death_rate
        self.start_trait, self.drift, self.diffusion = start_trait, drift, diffusion
        self.stop_value = stop_value
        self.sampling_prob = sampling_prob
        self.cond_surv, self.cond_spn = cond_surv, cond_spn
        self.cond_obs_both_sides = cond_obs_both_sides
        if not all(isinstance(v, bool) for v in (cond_surv, cond_spn, cond_obs_both_sides)):
            raise ValueError("Condition flags must be booleans.")
        if stop == "size" and (min_rec_taxa is not None or max_rec_taxa is not None):
            raise ValueError("Sampled-tip count bounds apply only to age stopping.")
        self.min_rec_taxa = 0 if min_rec_taxa is None else _integer(min_rec_taxa, "min_rec_taxa", 0)
        self.max_rec_taxa = (float("inf") if max_rec_taxa is None
                             else _integer(max_rec_taxa, "max_rec_taxa", 0))
        if self.min_rec_taxa > self.max_rec_taxa:
            raise ValueError("min_rec_taxa exceeds max_rec_taxa.")
        self.max_n_attempts = _integer(max_n_attempts, "max_n_attempts")
        self.init_check_vectorize_sample_size()
        if rng_seed is not None:
            np.random.seed(rng_seed)

    # A sample has one parameter set, shared by its independent replicates and retries.
    def init_check_vectorize_sample_size(self, param_list=None):
        for name in ("birth_rate", "death_rate", "start_trait", "drift", "diffusion", "stop_value", "sampling_prob"):
            raw = getattr(self, name)
            values = list(raw) if isinstance(raw, (list, tuple, np.ndarray)) else [raw]
            if len(values) not in (1, self.n_sim):
                raise ValueError(f"{name} must have length 1 or n={self.n_sim}.")
            if name in ("birth_rate", "death_rate"):
                if not all(isinstance(v, RateFunction) for v in values):
                    raise ValueError(f"{name} requires QuaSSE rate objects.")
            else:
                values = [float(v) for v in values]
                if not np.all(np.isfinite(values)):
                    raise ValueError(f"{name} must be finite.")
                if name in ("diffusion", "stop_value") and min(values) < 0:
                    raise ValueError(f"{name} must be nonnegative.")
                if name == "sampling_prob":
                    if min(values) < 0 or max(values) > 1:
                        raise ValueError("sampling_prob must be between zero and one.")
                    if self.stop == "size" and any(v != 1 for v in values):
                        raise ValueError("Size stopping requires complete sampling.")
                if name == "stop_value" and self.stop == "size":
                    values = [_integer(v, name) for v in values]
                    if max(values) > self.max_alive:
                        raise ValueError("Size target exceeds max_alive.")
            setattr(self, name, values * self.n_sim if len(values) == 1 else values)

    # Nodes are the only persistent simulation representation; internal traits freeze at birth.
    def _new_node(self, tree, label, trait, parent=None):
        node = dp.Node(label=label, edge_length=0.0,
                       taxon=tree.taxon_namespace.require_taxon(label=label))
        node.trait = float(trait)
        node.alive = node.sampled = True
        node.is_sa = node.is_sa_dummy_parent = node.is_sa_lineage = False
        node.annotations.add_bound_attribute("trait")
        if parent is not None:
            parent.add_child(node)
        return node

    # Enforce a common deadline inside attempts, so a single large tree cannot bypass it.
    def _check_runtime(self, deadline):
        if time.monotonic() >= deadline:
            raise ec.RunTimeLimit(self.runtime_limit)

    # Each step freezes rates at the current traits, selects at most one event, then evolves
    # the resulting living nodes. A full step has expected event count R/(kR) = 1/k.
    # Truncating a final step also scales this probability. Equal horizon/step times finish
    # the simulation; size stopping happens at the selected birth, before any trait update.
    def _simulate_diversitree(self, sample_idx, deadline):
        birth, death = self.birth_rate[sample_idx], self.death_rate[sample_idx]
        drift, variance = self.drift[sample_idx], self.diffusion[sample_idx]
        target = self.stop_value[sample_idx]
        if self.max_nodes is not None and self.max_nodes < 2:
            raise ec.GenerateFailError(self.DN_NAME, "Resource limit max_nodes cannot hold origin and initial lineage.")
        tree = dp.Tree(is_rooted=True)
        origin = self._new_node(tree, "origin", self.start_trait[sample_idx])
        origin.alive = origin.sampled = False
        tree.seed_node = origin
        living = [self._new_node(tree, "brosc", origin.trait, origin)]
        elapsed = 0.0
        node_count = 0
        steps = 0
        zero_rates = birth.is_identically_zero and death.is_identically_zero
        if self.stop == "size" and zero_rates:
            raise ec.GenerateFailError(self.DN_NAME, "Zero rates cannot reach the size-stopping birth.")
        while living and (self.stop == "size" or elapsed < target):
            self._check_runtime(deadline)
            if steps >= self.max_steps:
                raise ec.GenerateFailError(self.DN_NAME, "Exceeded max_steps.")
            steps += 1
            traits = np.array([nd.trait for nd in living])
            lx, mx = birth(traits), death(traits)
            with np.errstate(over="ignore"):
                birth_sum, death_sum = float(np.sum(lx)), float(np.sum(mx))
                total = birth_sum + death_sum
            if not np.isfinite(total):
                raise ec.GenerateFailError(self.DN_NAME, "Nonfinite total event rate.")
            if total == 0 and not zero_rates and self.dt_max is None:
                raise ec.GenerateFailError(self.DN_NAME,
                                          "Zero total rate from nonzero functions: supply dt_max to advance traits.")
            remaining = target - elapsed if self.stop == "age" else float("inf")
            full_step = (1.0 / self.k) / total if total else float("inf")
            # Locally zero rates need periodic trait updates; identically zero rates can
            # use a single exact Brownian transition. Shortened steps have event chance R*dt.
            dt_limit = self.dt_max if self.dt_max is not None and not zero_rates else float("inf")
            dt = min(full_step, dt_limit, remaining)
            if not np.isfinite(dt) or dt <= 0 or elapsed + dt == elapsed:
                raise ec.GenerateFailError(self.DN_NAME, "Cannot advance simulation time.")
            probability = 1.0 / self.k if dt == full_step else total * dt
            final_step = self.stop == "age" and dt == remaining
            if total and np.random.random() < probability:
                is_birth = np.random.random() < birth_sum / total
                weights = lx if is_birth else mx
                chosen = int(np.random.choice(len(living), p=weights / np.sum(weights)))
                if is_birth and self.stop == "size" and len(living) == target:
                    break
                if is_birth and len(living) >= self.max_alive:
                    raise ec.GenerateFailError(self.DN_NAME, "Exceeded max_alive.")
                # node_count counts all daughters ever allocated, including extinct lineages;
                # origin and brosc add two more. Check before allocating either new daughter.
                if is_birth and self.max_nodes is not None and node_count + 4 > self.max_nodes:
                    raise ec.GenerateFailError(self.DN_NAME, "Exceeded resource limit max_nodes.")
                node = living.pop(chosen)
                node.alive = node.sampled = False
                if is_birth:
                    if node.label == "brosc":
                        node.label = node.taxon.label = "root"
                    for _ in range(2):
                        node_count += 1
                        living.append(self._new_node(tree, f"nd{node_count}", node.trait, node))
            if not living:
                break
            # Brownian variance is variance-per-time multiplied by elapsed duration.
            # No draw is needed for deterministic evolution, including zero diffusion.
            with np.errstate(over="ignore", invalid="ignore"):
                increments = (np.random.normal(drift * dt, np.sqrt(variance * dt), len(living))
                              if variance else np.full(len(living), drift * dt))
                new_traits = np.array([nd.trait for nd in living]) + increments
            if not np.all(np.isfinite(new_traits)):
                raise ec.GenerateFailError(self.DN_NAME, "Nonfinite Brownian trait update.")
            for node, trait in zip(living, new_traits):
                node.edge_length += dt
                node.trait = float(trait)
            elapsed = target if final_step else elapsed + dt
        horizon = float(target if self.stop == "age" else elapsed)
        return tree, horizon

    # Dispatch without changing the tree contract when additional methods are introduced.
    def simulate(self, sample_idx=0, deadline=None):
        if deadline is None:
            deadline = time.monotonic() + self.runtime_limit
        tree, horizon = self._simulate_diversitree(sample_idx, deadline)
        probability = self.sampling_prob[sample_idx]
        for node in tree.leaf_node_iter():
            node.sampled = node.alive and (probability == 1 or
                                          (probability > 0 and np.random.random() < probability))
        return AnnotatedTree(tree, ContinuousTrait(), start_at_origin=True, max_age=horizon,
                             condition_on_obs_both_sides_root=self.cond_obs_both_sides,
                             tree_died=not any(nd.alive for nd in tree.leaf_node_iter()),
                             tree_invalid=False)

    # Check observations on the complete tree: pruning can replace the original root.
    def _is_tree_accepted(self, tree):
        if self.cond_surv and tree.tree_died:
            return False
        if self.cond_spn and tree.root_node is None:
            return False
        if self.cond_obs_both_sides:
            if tree.root_node is None or not all(
                    any(nd.alive and nd.sampled for nd in child.leaf_iter())
                    for child in tree.root_node.child_node_iter()):
                return False
        return (self.stop != "age" or
                self.min_rec_taxa <= tree.n_extant_sampled_terminal_nodes <= self.max_rec_taxa)

    # Rejection resamples only the tree, holding each sample's model parameters fixed.
    # Numerical and resource errors propagate; treating them as rejection would add conditioning.
    def generate(self):
        deadline = time.monotonic() + self.runtime_limit
        output = []
        rejected = 0
        for sample in range(self.n_sim):
            for _ in range(self.n_repl):
                while True:
                    self._check_runtime(deadline)
                    tree = self.simulate(sample, deadline)
                    self._check_runtime(deadline)
                    if self._is_tree_accepted(tree):
                        output.append(tree)
                        break
                    rejected += 1
                    if rejected >= self.max_n_attempts:
                        raise ec.MaxNFailedAttemptsLimit(self.max_n_attempts)
        return output

    def get_rev_inference_spec_info(self):
        raise NotImplementedError("QuaSSE RevBayes inference export is not supported.")
