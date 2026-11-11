"""Values from model replicates of one shared dependency graph."""
from phylojunction.utility.exception_classes import ScriptEvaluationError
import copy
import random
import numpy as np
from phylojunction.pgm.pgm import DistrForSampling, _Condition
from phylojunction.interface.cmdbox.expressions import evaluate, value_type
from phylojunction.interface.grammar.constructors import construct


class ExecutionResults:
    """One graph, with values and sampling distributions indexed by replicate and definition.

    Python replicate indices are zero-based. Script vector indices are one-based.
    Each result owns its graph, preventing later interactive additions from changing
    another result. Definitions are shared by all replicates within the result.
    """

    # Execute every definition in order for each model run. Failed execution restores
    # both RNGs; successful execution advances them, including when no seed is specified.
    def __init__(self, model, runs=1, random_seed=None, on_replicate_complete=None, max_attempts=1000):
        if type(runs) is not int or runs < 1:
            raise ValueError("The number of replicates must be a positive integer.")
        if max_attempts is not None and (type(max_attempts) is not int or max_attempts < 1):
            raise ValueError("max_attempts must be a positive integer or None.")
        self.model = copy.deepcopy(model)
        self.random_seed = random_seed
        self._values = [[] for _ in range(runs)]
        self._sampled_from = [[] for _ in range(runs)]
        state = random.getstate(), np.random.get_state()
        try:
            if random_seed is not None:
                random.seed(random_seed)
                np.random.seed(random_seed)
            self._execute_from(0, on_replicate_complete, max_attempts)
        except Exception:
            random.setstate(state[0])
            np.random.set_state(state[1])
            raise

    @property
    def n_executions(self):
        return len(self._values)

    # Each rejected attempt discards all its values but keeps advancing both RNGs.
    # Append executes only new statements, which cannot contain acceptance conditions.
    def _execute_from(self, start, on_replicate_complete=None, max_attempts=1000):
        statements = self.model.statements[start:]
        for execution, (row, sampled) in enumerate(zip(self._values, self._sampled_from), 1):
            attempt = 0
            while True:
                attempt += 1
                for node in statements:
                    try:
                        value = evaluate(node.expression, lambda name: row[node.bindings[name]], construct)
                        if isinstance(node, _Condition):
                            if type(value) is not bool:
                                raise TypeError('condition(...) requires a Boolean scalar.')
                            if value:
                                continue
                            if max_attempts is not None and attempt >= max_attempts:
                                raise ValueError(f'No accepted replicate after {attempt} attempts; '
                                                 'increase max_attempts or revise the condition.')
                            row.clear()
                            sampled.clear()
                            break
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
                            line_number=node.line_number, variable=getattr(node, 'node_name', None),
                            execution=execution) from error
                else:
                    break
            if on_replicate_complete is not None:
                on_replicate_complete(execution)

    # Reject Python's implicit negative indexing for replicate selection as well.
    def value(self, name, execution_index=0):
        if type(execution_index) is not int or not 0 <= execution_index < self.n_executions:
            raise IndexError("Replicate index out of range.")
        return self._values[execution_index][self.model.name_node_dict[name].definition_id]

    def values(self, name):
        return [self.value(name, i) for i in range(self.n_executions)]

    # Interactive additions are transactional snapshots, not mutations of a shared graph.
    # Only new definitions run: earlier descendants retain their previously sampled values.
    def append(self, source, source_name="<interactive>"):
        candidate = copy.deepcopy(self)
        start = len(candidate.model.statements)
        state = random.getstate(), np.random.get_state()
        try:
            for number, line in enumerate(source.splitlines(), 1):
                candidate.model.add_definition(line, source_name, number)
            for statement in candidate.model.statements[start:]:
                if isinstance(statement, _Condition):
                    raise ScriptEvaluationError(
                        'Cannot append conditions to executed results; rerun the complete model.',
                        source=statement.source, source_name=statement.source_name,
                        line_number=statement.line_number)
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
