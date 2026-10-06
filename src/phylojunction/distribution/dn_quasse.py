"""QuaSSE simulation with interchangeable numerical methods.

The initial method follows diversitree/R/simulate-quasse.R. Events precede trait
updates. The final age-limited interval and its event probability are shortened;
size stopping omits the terminating birth and its subsequent trait update.
"""
from phylojunction.data.trait import ContinuousTrait

import time
import math
import heapq
from dataclasses import dataclass
from itertools import count
import numpy as np
import dendropy as dp

from phylojunction.calculation.continuous_sse import RateFunction, ConstantRate
from phylojunction.data.tree import AnnotatedTree
from phylojunction.pgm.pgm import DistrForSampling
from phylojunction.utility import exception_classes as ec


# Reject fractional counts rather than silently truncating Python or grammar inputs.
def _integer(value, name, minimum=1):
    number = float(value)
    if not np.isfinite(number) or not number.is_integer() or number < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}.")
    return int(number)


@dataclass
class _LineageInterval:
    """One unfinished Brownian continuation; current trait/edge remain on the node."""
    node: dp.Node
    time: float
    end_time: float
    end_trait: float
    rate_bound: float


class _LocalBoundExceeded(Exception):
    """Only a finite rate exceeding a finite proposal bound permits numerical retry."""


class DnQuaSSE(DistrForSampling):
    DN_NAME = "DnQuaSSE"

    # Keep model vectors separate from scalar execution settings, then validate both.
    def __init__(self, birth_rate, death_rate, stop, stop_value, start_trait=0.0,
                 drift=0.0, diffusion=0.0, n=1, nr=1, method="diversitree", k=500,
                 runtime_limit=300, max_steps=1000000, max_alive=100000, rng_seed=None,
                 sampling_prob=1.0, cond_surv=True, cond_spn=False, cond_obs_both_sides=False,
                 min_rec_taxa=None, max_rec_taxa=None, max_n_attempts=200, dt_max=None, max_nodes=None,
                 block_duration=None, bridge_error=None, fossil_rate=None):
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
        if method not in ("diversitree", "local_thinning"):
            raise ValueError(f"Unknown QuaSSE method {method!r}; supported: diversitree, local_thinning.")
        self.block_duration = None if block_duration is None else float(block_duration)
        self.bridge_error = None if bridge_error is None else float(bridge_error)
        if method == "local_thinning":
            if self.block_duration is None or not math.isfinite(self.block_duration) or self.block_duration <= 0:
                raise ValueError("local_thinning requires a positive finite block_duration.")
            if self.bridge_error is None:
                self.bridge_error = 1e-8
            if not math.isfinite(self.bridge_error) or not 0 < self.bridge_error < 1:
                raise ValueError("bridge_error must be finite and strictly between zero and one.")
            if self.dt_max is not None or self.k != 500:
                raise ValueError("local_thinning uses block_duration, not dt_max or k (default k=500 is ignored).")
        elif self.block_duration is not None or self.bridge_error is not None:
            raise ValueError("block_duration and bridge_error apply only to local_thinning.")
        if stop not in ("age", "size"):
            raise ValueError('stop must be "age" or "size".')
        self.birth_rate, self.death_rate = birth_rate, death_rate
        self.fossil_rate = ConstantRate(0) if fossil_rate is None else fossil_rate
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
        for name in ("birth_rate", "death_rate", "fossil_rate", "start_trait", "drift", "diffusion",
                     "stop_value", "sampling_prob"):
            raw = getattr(self, name)
            values = list(raw) if isinstance(raw, (list, tuple, np.ndarray)) else [raw]
            if len(values) not in (1, self.n_sim):
                raise ValueError(f"{name} must have length 1 or n={self.n_sim}.")
            if name in ("birth_rate", "death_rate", "fossil_rate"):
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
                    if self.stop == "size" and any(v == 0 for v in values):
                        raise ValueError("Size stopping requires positive sampling_prob.")
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

    # Preserve the continuing lineage's identity/mark; the second child records the observation.
    # Child order distinguishes continuation from fossil even after pruning zero-length edges.
    def _record_fossil(self, tree, node, fossil_index):
        label, sampled = node.label, node.sampled
        node.label = f"dummy{fossil_index}"
        node.taxon = tree.taxon_namespace.require_taxon(label=node.label)
        node.alive = node.sampled = node.is_sa = node.is_sa_lineage = False
        node.is_sa_dummy_parent = True
        continuation = self._new_node(tree, label, node.trait, node)
        continuation.sampled, continuation.is_sa_lineage = sampled, True
        fossil = self._new_node(tree, f"sa{fossil_index}", node.trait, node)
        fossil.alive, fossil.is_sa = False, True
        return continuation

    # One daughter inherits and the other gets a fresh Bernoulli mark, so K increases by at most one.
    # Randomize orientation only for unequal marks; complete sampling must preserve old RNG order.
    def _birth_sampling_flags(self, parent_sampled, rho):
        if rho == 1:
            return True, True
        fresh = bool(np.random.random() < rho)
        if fresh == parent_sampled or np.random.random() < .5:
            return parent_sampled, fresh
        return fresh, parent_sampled

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
        fossil = self.fossil_rate[sample_idx]
        fossils_enabled = not fossil.is_identically_zero
        rho = self.sampling_prob[sample_idx]
        if self.max_nodes is not None and self.max_nodes < 2:
            raise ec.GenerateFailError(self.DN_NAME, "Resource limit max_nodes cannot hold origin and initial lineage.")
        tree = dp.Tree(is_rooted=True)
        origin = self._new_node(tree, "origin", self.start_trait[sample_idx])
        origin.alive = origin.sampled = False
        tree.seed_node = origin
        living = [self._new_node(tree, "brosc", origin.trait, origin)]
        if self.stop == "size":
            living[0].sampled = rho == 1 or bool(np.random.random() < rho)
        sampled_count = int(living[0].sampled)
        elapsed = 0.0
        node_count = fossil_count = 0
        steps = 0
        zero_rates = birth.is_identically_zero and death.is_identically_zero
        if self.stop == "size" and zero_rates:
            raise ec.GenerateFailError(self.DN_NAME, "Zero rates cannot reach the size-stopping birth.")
        zero_rates = zero_rates and not fossils_enabled
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
                if fossils_enabled:
                    fx = fossil(traits)
                    total += float(np.sum(fx))
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
                if fossils_enabled:
                    mark = np.random.random()
                    is_birth = mark < birth_sum / total
                    is_fossil = mark >= (birth_sum + death_sum) / total
                    weights = fx if is_fossil else (lx if is_birth else mx)
                else:
                    is_birth = np.random.random() < birth_sum / total
                    is_fossil = False
                    weights = lx if is_birth else mx
                chosen = int(np.random.choice(len(living), p=weights / np.sum(weights)))
                node = living[chosen]
                flags = (True, True)
                if is_birth and self.stop == "size":
                    flags = self._birth_sampling_flags(node.sampled, rho)
                    next_count = sampled_count - int(node.sampled) + sum(flags)
                    if next_count > target:
                        break
                    sampled_count = next_count
                if is_birth and len(living) >= self.max_alive:
                    raise ec.GenerateFailError(self.DN_NAME, "Exceeded max_alive.")
                # Each birth/fossil adds two nodes, in addition to origin and initial lineage.
                if ((is_birth or is_fossil) and self.max_nodes is not None
                    and node_count + 2 * fossil_count + 4 > self.max_nodes):
                    raise ec.GenerateFailError(self.DN_NAME, "Exceeded resource limit max_nodes.")
                if is_fossil:
                    fossil_count += 1
                    living[chosen] = self._record_fossil(tree, node, fossil_count)
                else:
                    living.pop(chosen)
                    if not is_birth:
                        sampled_count -= int(node.sampled)
                    node.alive = node.sampled = False
                    if is_birth:
                        if node.label == "brosc":
                            node.label = node.taxon.label = "root"
                        for flag in flags:
                            node_count += 1
                            daughter = self._new_node(tree, f"nd{node_count}", node.trait, node)
                            daughter.sampled = flag
                            living.append(daughter)
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

    # Reveal the existing bridge at a later time, retaining its endpoint and all past draws.
    # Given z at t and y at e, the mean is (1-f)z+fy and variance v*u*(1-f), f=u/(e-t).
    # Drift is already in y. Equal times do nothing; the endpoint is assigned without a draw.
    def _advance_bridge(self, interval, new_time, variance):
        if not math.isfinite(new_time) or not interval.time <= new_time <= interval.end_time:
            raise ec.GenerateFailError(self.DN_NAME, "Invalid Brownian bridge time.")
        if new_time == interval.time:
            return
        duration = new_time - interval.time
        if new_time == interval.end_time:
            trait = interval.end_trait
        else:
            fraction = duration / (interval.end_time - interval.time)
            trait = (1 - fraction) * interval.node.trait + fraction * interval.end_trait
            if variance:
                trait += math.sqrt(variance * duration * (1 - fraction)) * np.random.normal()
        if not math.isfinite(trait):
            raise ec.GenerateFailError(self.DN_NAME, "Nonfinite Brownian bridge trait.")
        interval.node.trait = float(trait)
        interval.node.edge_length += duration
        interval.time = new_time

    # Allocate epsilon_j=delta/[j(j+1)] before sampling an unconditional endpoint.
    # A bridge's two-sided deviation exceeds r with probability <=2*exp(-2*r*r/(v*h)).
    # Enlarging the endpoint range by r therefore gives a likely region, not a hard barrier.
    def _start_bridge_interval(self, node, time, end_time, birth, death, drift, variance, block_index, fossil_rate=None):
        if not math.isfinite(end_time) or end_time <= time:
            raise ec.GenerateFailError(self.DN_NAME, "Cannot advance Brownian block time.")
        duration = end_time - time
        scale = math.sqrt(variance * duration)
        endpoint = node.trait + drift * duration
        if variance:
            endpoint += scale * np.random.normal()
        # Summable allocation survives arbitrary numbers of blocks and retries. Computing
        # its log directly avoids underflow of epsilon_j; the index belongs to the full call.
        log_inverse = math.log(2) - math.log(self.bridge_error) + math.log(block_index) + math.log(block_index + 1)
        radius = scale * math.sqrt(log_inverse / 2)
        left, right = min(node.trait, endpoint) - radius, max(node.trait, endpoint) + radius
        if not all(math.isfinite(x) for x in (endpoint, left, right)):
            raise ec.GenerateFailError(self.DN_NAME, "Nonfinite Brownian endpoint or region.")
        bound = birth.bound_on_interval(left, right) + death.bound_on_interval(left, right)
        if fossil_rate is not None:
            bound += fossil_rate.bound_on_interval(left, right)
        if not math.isfinite(bound):
            raise ec.GenerateFailError(self.DN_NAME, "Nonfinite total proposal bound.")
        return _LineageInterval(node, time, end_time, float(endpoint), bound)

    # Each living node owns one bridge and one queued candidate/endpoint. Pop chronologically,
    # reveal only that trait, then thin at its actual rate. Serial numbers break ties; a
    # candidate at the endpoint is an endpoint action. Birth/death removes the popped record,
    # so there are no stale entries. Conditional on staying in its region, Poisson thinning
    # is exact; the call-wide sum of excursion probabilities bounds first divergence.
    def _simulate_local_thinning(self, sample_idx, deadline, block_indices, work_indices):
        birth, death = self.birth_rate[sample_idx], self.death_rate[sample_idx]
        drift, variance = self.drift[sample_idx], self.diffusion[sample_idx]
        target = self.stop_value[sample_idx]
        fossil = self.fossil_rate[sample_idx]
        fossils_enabled = not fossil.is_identically_zero
        rho = self.sampling_prob[sample_idx]
        self._check_runtime(deadline)
        if self.max_nodes is not None and self.max_nodes < 2:
            raise ec.GenerateFailError(self.DN_NAME, "Resource limit max_nodes cannot hold origin and initial lineage.")
        tree = dp.Tree(is_rooted=True)
        origin = self._new_node(tree, "origin", self.start_trait[sample_idx])
        origin.alive = origin.sampled = False
        tree.seed_node = origin
        initial = self._new_node(tree, "brosc", origin.trait, origin)
        if self.stop == "size":
            initial.sampled = rho == 1 or bool(np.random.random() < rho)
        sampled_count = int(initial.sampled)
        if self.stop == "age" and target == 0:
            return tree, 0.0
        zero_biological = birth.is_identically_zero and death.is_identically_zero
        if zero_biological and self.stop == "size":
            raise ec.GenerateFailError(self.DN_NAME, "Zero rates cannot reach the size-stopping birth.")
        if zero_biological and not fossils_enabled:
            trait = initial.trait + drift * target
            if variance:
                trait += math.sqrt(variance * target) * np.random.normal()
            if not math.isfinite(trait):
                raise ec.GenerateFailError(self.DN_NAME, "Nonfinite Brownian trait update.")
            initial.trait, initial.edge_length = float(trait), float(target)
            return tree, float(target)

        living = {}
        queue = []
        serials = count()

        # Queue an exponential proposal or the fixed endpoint. No accepted future event is
        # sampled here: size stopping must still be able to reveal every pending bridge.
        def schedule(interval):
            candidate_time = interval.end_time
            is_candidate = False
            if interval.rate_bound:
                wait = float(np.random.exponential(1.0 / interval.rate_bound))
                if math.isnan(wait) or wait <= 0:
                    raise ec.GenerateFailError(self.DN_NAME, "Cannot advance candidate time.")
                if wait < interval.end_time - interval.time:
                    candidate_time = min(interval.time + wait, interval.end_time)
                    if candidate_time <= interval.time:
                        raise ec.GenerateFailError(self.DN_NAME, "Cannot advance candidate time.")
                    is_candidate = candidate_time < interval.end_time
            heapq.heappush(queue, (candidate_time, next(serials), is_candidate, interval))

        # The same block creation is used for initialization, endpoint renewal, and daughters.
        # Consuming the index before drawing prevents endpoint-dependent budget allocation.
        def start(node, time):
            end_time = time + self.block_duration
            if self.stop == "age":
                end_time = min(end_time, target)
            fossil_args = {"fossil_rate": fossil} if fossils_enabled else {}
            interval = self._start_bridge_interval(node, time, end_time, birth, death,
                                                   drift, variance, next(block_indices), **fossil_args)
            living[node] = interval
            schedule(interval)

        start(initial, 0.0)
        node_count = fossil_count = 0
        elapsed = 0.0
        while queue:
            self._check_runtime(deadline)
            if next(work_indices) >= self.max_steps:
                raise ec.GenerateFailError(self.DN_NAME, "Exceeded max_steps (local queue actions).")
            elapsed, _, is_candidate, interval = heapq.heappop(queue)
            self._advance_bridge(interval, elapsed, variance)
            node = interval.node
            if not is_candidate:
                if self.stop != "age" or elapsed < target:
                    start(node, elapsed)
                continue
            lx, mx = float(birth(node.trait)), float(death(node.trait))
            fx = float(fossil(node.trait)) if fossils_enabled else 0.0
            total = lx + mx + fx if fossils_enabled else lx + mx
            if not all(math.isfinite(x) and x >= 0 for x in (lx, mx, fx, total)):
                raise ec.GenerateFailError(self.DN_NAME, "Nonfinite or negative event rate.")
            if total > interval.rate_bound:
                raise _LocalBoundExceeded(f"At time {elapsed}, trait {node.trait}: rate {total} "
                                          f"exceeds bound {interval.rate_bound}.")
            mark = np.random.uniform(0, interval.rate_bound)
            if mark >= total:
                schedule(interval)
                continue
            is_birth = mark < lx
            flags = (True, True)
            next_count = sampled_count
            if is_birth and self.stop == "size":
                flags = self._birth_sampling_flags(node.sampled, rho)
                next_count = sampled_count - int(node.sampled) + sum(flags)
            if is_birth and self.stop == "size" and next_count > target:
                # Other nodes have revealed no events after this time. Conditional on their
                # endpoints/past traits, advance each existing bridge to the stopping birth.
                for pending in living.values():
                    self._check_runtime(deadline)
                    self._advance_bridge(pending, elapsed, variance)
                return tree, elapsed
            if is_birth and len(living) >= self.max_alive:
                raise ec.GenerateFailError(self.DN_NAME, "Exceeded max_alive.")
            is_fossil = mark >= lx + mx
            if ((is_birth or is_fossil) and self.max_nodes is not None
                    and node_count + 2 * fossil_count + 4 > self.max_nodes):
                raise ec.GenerateFailError(self.DN_NAME, "Exceeded resource limit max_nodes.")
            del living[node]
            if is_fossil:
                fossil_count += 1
                continuation = self._record_fossil(tree, node, fossil_count)
                # Sampling reveals the existing bridge; its endpoint, bound and budget remain valid.
                interval.node = continuation
                living[continuation] = interval
                schedule(interval)
                continue
            sampled_count = next_count if is_birth else sampled_count - int(node.sampled)
            node.alive = node.sampled = False
            if is_birth:
                if node.label == "brosc":
                    node.label = node.taxon.label = "root"
                for flag in flags:
                    node_count += 1
                    daughter = self._new_node(tree, f"nd{node_count}", node.trait, node)
                    daughter.sampled = flag
                    start(daughter, elapsed)
        return tree, float(target if self.stop == "age" else elapsed)

    # Dispatch without changing downstream tree/observation contracts or default RNG order.
    # Numerical retries discard only the current tree; budget, work, RNG and deadline persist.
    def simulate(self, sample_idx=0, deadline=None, *, _block_indices=None):
        if deadline is None:
            deadline = time.monotonic() + self.runtime_limit
        if self.method == "diversitree":
            tree, horizon = self._simulate_diversitree(sample_idx, deadline)
        else:
            block_indices = count(1) if _block_indices is None else _block_indices
            work_indices = count()
            while True:
                try:
                    tree, horizon = self._simulate_local_thinning(sample_idx, deadline, block_indices, work_indices)
                    break
                except _LocalBoundExceeded:
                    # A detected violation already implies an excursion charged to this call's
                    # budget. Restarting is approximate recovery, not exact rejection sampling.
                    continue
        probability = self.sampling_prob[sample_idx]
        has_fossils = False
        for node in tree.leaf_node_iter():
            if node.is_sa:
                has_fossils = True
            elif self.stop == "age":
                node.sampled = node.alive and (probability == 1 or
                                              (probability > 0 and np.random.random() < probability))
        if has_fossils:
            for node in tree:
                node.annotations.add_bound_attribute("is_sa")
                node.annotations.add_bound_attribute("is_sa_dummy_parent")
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
    # Other numerical/resource errors propagate. Bound-exceedance retries happen inside simulate.
    def generate(self):
        deadline = time.monotonic() + self.runtime_limit
        output = []
        rejected = 0
        block_indices = count(1)
        for sample in range(self.n_sim):
            for _ in range(self.n_repl):
                while True:
                    self._check_runtime(deadline)
                    tree = (self.simulate(sample, deadline, _block_indices=block_indices)
                            if self.method == "local_thinning" else self.simulate(sample, deadline))
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
