"""Structural expressions shared by script parsing and model evaluation.

Only the explicitly supported nodes are interpreted; Python code is never evaluated.
Variable names in arguments are dependencies, whereas constructor names are dispatch keys.
"""
import ast
import re


# Preserve the language's assignments while delegating nested expression syntax to Python's AST.
def parse_statement(source):
    source = source.strip()
    if not source or source.startswith('#'):
        return None
    match = re.fullmatch(r'([A-Za-z][A-Za-z_0-9]*)\s*(<-|:=|~)\s*(.+)', source)
    if match is None:
        raise ValueError('Expected name <- value, name := expression, or name ~ distribution.')
    name, operator, text = match.groups()
    expression = parse_expression(text)
    return name, operator, expression


# Reject executable Python constructs before a model is compiled or any values are sampled.
def parse_expression(source):
    expression = ast.parse(source, mode='eval').body
    validate(expression)
    return expression


# Binding patterns are syntax, not list expressions. Flat tuple unpacking supports zip
# without adding arbitrary Python tuples, assignments, or iterator execution to PJ.
def binding_names(target):
    targets = target.elts if isinstance(target, ast.Tuple) else [target]
    if not targets or any(not isinstance(item, ast.Name) for item in targets):
        raise ValueError("A comprehension target must be a name or flat tuple of names.")
    names = [item.id for item in targets]
    if len(set(names)) != len(names) or any(name in ('true', 'false') for name in names):
        raise ValueError("Comprehension binding names must be distinct and cannot be true or false.")
    return names


# Only the explicitly supported subset of comprehensions enters evaluation or export.
def comprehension(expression):
    if len(expression.generators) != 1:
        raise ValueError("Use nested comprehensions instead of multiple for clauses.")
    generator = expression.generators[0]
    if generator.is_async or generator.ifs:
        raise ValueError("Comprehension filters and asynchronous iteration are not supported.")
    binding_names(generator.target)
    return generator


# Validate structure recursively. Domain-specific argument checks belong to constructors.
def validate(expression):
    if isinstance(expression, ast.Constant) and type(expression.value) in (str, int, float, bool):
        return
    if isinstance(expression, ast.Name):
        return
    if isinstance(expression, ast.List):
        for item in expression.elts:
            validate(item)
        return
    if isinstance(expression, ast.ListComp):
        generator = comprehension(expression)
        validate(generator.iter)
        validate(expression.elt)
        return
    if isinstance(expression, ast.UnaryOp) and isinstance(expression.op, (ast.UAdd, ast.USub)):
        if isinstance(expression.operand, ast.Constant) and type(expression.operand.value) in (int, float):
            return
    if isinstance(expression, ast.Subscript):
        validate(expression.value)
        validate(expression.slice)
        return
    if isinstance(expression, ast.Call) and isinstance(expression.func, ast.Name):
        seen = set()
        for keyword in expression.keywords:
            if keyword.arg is None or keyword.arg in seen:
                raise ValueError('Duplicate keyword arguments and ** expansion are not supported.')
            seen.add(keyword.arg)
            validate(keyword.value)
        for argument in expression.args:
            validate(argument)
        return
    raise ValueError('Expected a literal, variable, vector, indexing expression, or named constructor.')


# Dependency discovery excludes call targets and includes every collection/index expression.
def variable_names(expression, local=frozenset()):
    if isinstance(expression, ast.Name):
        if expression.id not in ('true', 'false') and expression.id not in local:
            yield expression.id
    elif isinstance(expression, ast.ListComp):
        generator = comprehension(expression)
        yield from variable_names(generator.iter, local)
        yield from variable_names(expression.elt, local | set(binding_names(generator.target)))
    elif isinstance(expression, ast.Call):
        for argument in expression.args:
            yield from variable_names(argument, local)
        for keyword in expression.keywords:
            yield from variable_names(keyword.value, local)
    else:
        for child in ast.iter_child_nodes(expression):
            yield from variable_names(child, local)


# Evaluate the validated AST against one execution's values. Subscripts are always
# one-based language indices; Python negative indexing and slicing are not inherited.
def evaluate(expression, lookup, construct):
    if isinstance(expression, ast.Constant):
        return expression.value
    if isinstance(expression, ast.Name):
        if expression.id in ("true", "false"):
            return expression.id == "true"
        return lookup(expression.id)
    if isinstance(expression, ast.List):
        value = [evaluate(item, lookup, construct) for item in expression.elts]
        value_type(value)
        return value
    if isinstance(expression, ast.ListComp):
        generator = comprehension(expression)
        iterable = evaluate(generator.iter, lookup, construct)
        if not isinstance(iterable, list):
            raise TypeError("Comprehensions require a list iterable.")
        names = binding_names(generator.target)
        result = []
        for item in iterable:
            if isinstance(generator.target, ast.Tuple):
                if not isinstance(item, (list, tuple)) or len(item) != len(names):
                    raise ValueError("Comprehension unpacking requires exactly one entry per binding.")
                bindings = dict(zip(names, item))
            else:
                bindings = {names[0]: item}
            result.append(evaluate(expression.elt,
                                   lambda name: bindings[name] if name in bindings else lookup(name), construct))
        value_type(result)
        return result
    if isinstance(expression, ast.UnaryOp):
        value = expression.operand.value
        return -value if isinstance(expression.op, ast.USub) else value
    if isinstance(expression, ast.Subscript):
        vector = evaluate(expression.value, lookup, construct)
        index = evaluate(expression.slice, lookup, construct)
        if not isinstance(vector, (list, tuple)):
            raise TypeError("Only vectors can be indexed.")
        if type(index) is not int:
            raise TypeError("A vector index must be an integer.")
        if not 1 <= index <= len(vector):
            raise IndexError(f"Index {index} is outside the range 1..{len(vector)}.")
        return vector[index - 1]
    if isinstance(expression, ast.Call):
        args = [evaluate(arg, lookup, construct) for arg in expression.args]
        kwargs = {item.arg: evaluate(item.value, lookup, construct) for item in expression.keywords}
        if expression.func.id == 'zip_strict':
            if kwargs or len(args) < 2 or any(not isinstance(arg, list) for arg in args):
                raise TypeError("zip_strict requires at least two positional list arguments.")
            if len({len(arg) for arg in args}) != 1:
                raise ValueError("zip_strict requires equal lengths.")
            return list(zip(*args))
        result = construct(expression.func.id, *args, **kwargs)
        value_type(result)
        return result
    raise ValueError("Unsupported expression.")


# Empty lists carry an unknown element type; nonempty siblings constrain it. Numeric
# values share a type, but Boolean values do not. Distribution result types are deferred:
# distribution objects are compatible here without making claims about their draws.
def value_type(value):
    from phylojunction.pgm.pgm import DistrForSampling
    from phylojunction.calculation.continuous_sse import RateFunction
    if isinstance(value, list):
        element = None
        for item in value:
            element = merge_types(element, value_type(item))
        return ('vector', element)
    if isinstance(value, tuple):
        return ('record', tuple(value_type(item) for item in value))
    if isinstance(value, DistrForSampling):
        return ('distribution',)
    if isinstance(value, RateFunction):
        return ('rate',)
    if type(value) in (int, float):
        return ('number',)
    return (type(value),)


# Joining recursive list types checks nesting as well as leaf compatibility, not lengths.
# Tuples are zip records, whose fields may differ but must match corresponding fields.
def merge_types(left, right):
    if left is None:
        return right
    if right is None or left == right:
        return left
    if left[0] == right[0] == 'vector':
        return ('vector', merge_types(left[1], right[1]))
    if left[0] == right[0] == 'record' and len(left[1]) == len(right[1]):
        return ('record', tuple(merge_types(a, b) for a, b in zip(left[1], right[1])))
    raise TypeError("List elements must have compatible types and nesting; lengths may differ.")
