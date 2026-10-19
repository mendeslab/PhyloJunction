"""Values from complete executions of one shared dependency graph."""
from phylojunction.utility.exception_classes import ScriptEvaluationError
import copy
import random
import numpy as np
from phylojunction.pgm.pgm import DistrForSampling
from phylojunction.interface.cmdbox.expressions import evaluate, value_type
from phylojunction.interface.grammar.constructors import construct


class ExecutionResults:
    """One graph, with values and sampling distributions indexed by execution and definition.

    Python execution indices are zero-based. Script vector indices are one-based.
    Each result owns its graph, preventing later interactive additions from changing
    another result. Definitions are shared by all executions within the result.
    """

    # Execute every definition in order for each model run. Failed execution restores
    # both RNGs; successful execution advances them, including when no seed is specified.
    def __init__(self, model, runs=1, random_seed=None):
        if type(runs) is not int or runs < 1:
            raise ValueError("The number of model executions must be a positive integer.")
        self.model = copy.deepcopy(model)
        self.random_seed = random_seed
        self._values = [[] for _ in range(runs)]
        self._sampled_from = [[] for _ in range(runs)]
        state = random.getstate(), np.random.get_state()
        try:
            if random_seed is not None:
                random.seed(random_seed)
                np.random.seed(random_seed)
            self._execute_from(0)
        except Exception:
            random.setstate(state[0])
            np.random.set_state(state[1])
            raise

    @property
    def n_executions(self):
        return len(self._values)

    # Bind lookup to definition IDs, so a reassigned name cannot change prior expressions.
    def _execute_from(self, start):
        for execution, (row, sampled) in enumerate(zip(self._values, self._sampled_from), 1):
            for node in self.model.definitions[start:]:
                try:
                    value = evaluate(node.expression, lambda name: row[node.bindings[name]], construct)
                    distribution = None
                    if node.operator == '~':
                        validate_distributions(value)
                        distribution = value
                        value = draw_distributions(distribution)
                        value_type(value)
                    row.append(value)
                    sampled.append(distribution)
                except ScriptEvaluationError:
                    raise
                except Exception as error:
                    raise ScriptEvaluationError(
                        str(error), source=node.source, source_name=node.source_name,
                        line_number=node.line_number, variable=node.node_name,
                        execution=execution) from error

    # Reject Python's implicit negative indexing for execution selection as well.
    def value(self, name, execution_index=0):
        if type(execution_index) is not int or not 0 <= execution_index < self.n_executions:
            raise IndexError("Execution index out of range.")
        return self._values[execution_index][self.model.name_node_dict[name].definition_id]

    def values(self, name):
        return [self.value(name, i) for i in range(self.n_executions)]

    # Interactive additions are transactional snapshots, not mutations of a shared graph.
    # Only new definitions run: earlier descendants retain their previously sampled values.
    def append(self, source, source_name="<interactive>"):
        candidate = copy.deepcopy(self)
        start = len(candidate.model.definitions)
        state = random.getstate(), np.random.get_state()
        try:
            for number, line in enumerate(source.splitlines(), 1):
                candidate.model.add_definition(line, source_name, number)
            candidate._execute_from(start)
        except Exception:
            random.setstate(state[0])
            np.random.set_state(state[1])
            raise
        return candidate


# Check all leaves before sampling, so an invalid later element never starts partial draws.
def validate_distributions(value):
    if isinstance(value, list):
        for item in value:
            validate_distributions(item)
    elif not isinstance(value, DistrForSampling):
        raise TypeError("~ requires a distribution or nested lists of distributions, without ordinary values.")


# Recursion applies to specifications only: a distribution's draw is already a value.
def draw_distributions(value):
    if isinstance(value, list):
        return [draw_distributions(item) for item in value]
    return value.draw()
