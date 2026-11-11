"""Structural expressions shared by script parsing and model evaluation.

Only the explicitly supported nodes are interpreted; Python code is never evaluated.
Variable names in arguments are dependencies, whereas constructor names are dispatch keys.
"""
import ast
import io
import operator
import re
import tokenize


_BINARY_OPERATORS = {ast.Add: operator.add, ast.Sub: operator.sub,
                     ast.Mult: operator.mul, ast.Div: operator.truediv}
_UNARY_OPERATORS = {ast.UAdd: operator.pos, ast.USub: operator.neg}
_COMPARISONS = {ast.Lt: operator.lt, ast.LtE: operator.le, ast.Gt: operator.gt,
                ast.GtE: operator.ge, ast.Eq: operator.eq, ast.NotEq: operator.ne}


# Translate only tokens, leaving strings intact. Wrap each ! operand through comparisons
# but stop at logical operators/delimiters, matching R precedence even in x == !y.
def _python_expression(source):
    raw = list(tokenize.generate_tokens(io.StringIO(source).readline))
    tokens = []
    index = 0
    while index < len(raw):
        token = raw[index]
        index += 1
        if token.type in (tokenize.ENDMARKER, tokenize.NEWLINE, tokenize.NL, tokenize.COMMENT):
            continue
        if token.type == tokenize.ERRORTOKEN and token.string.isspace():
            continue
        if token.type == tokenize.NAME and token.string in ('and', 'or', 'not'):
            raise ValueError('Use &&, || and ! for Boolean operators.')
        if token.type == tokenize.OP and token.string in ('&', '|'):
            if index >= len(raw) or raw[index].string != token.string or raw[index].start != token.end:
                raise ValueError('Use scalar && and ||; vectorized & and | are not supported.')
            index += 1
            tokens.append('and' if token.string == '&' else 'or')
        else:
            tokens.append(token.string)
    closing = {}
    for index, token in enumerate(tokens):
        if token != '!':
            continue
        depth = 0
        end = index + 1
        while end < len(tokens):
            following = tokens[end]
            if depth == 0 and following in ('and', 'or', ',', ')', ']', 'for', 'if'):
                break
            if following in ('(', '['):
                depth += 1
            elif following in (')', ']'):
                depth -= 1
            end += 1
        closing[end] = closing.get(end, 0) + 1
    pieces = []
    for index, token in enumerate(tokens):
        pieces.extend(')' for _ in range(closing.get(index, 0)))
        pieces.append('( not' if token == '!' else token)
    pieces.extend(')' for _ in range(closing.get(len(tokens), 0)))
    return ' '.join(pieces)


# Apply only scalar numeric operations, sharing the same rules with symbolic export.
# Exact type checks exclude booleans and Python sequence concatenation/repetition.
def arithmetic(operation, *operands):
    if any(type(value) not in (int, float) for value in operands):
        raise TypeError("Arithmetic requires numeric scalar operands (not booleans).")
    operators = _UNARY_OPERATORS if len(operands) == 1 else _BINARY_OPERATORS
    return operators[type(operation)](*operands)


# Preserve the language's assignments while delegating nested expression syntax to Python's AST.
def parse_statement(source):
    source = source.strip()
    if not source or source.startswith('#'):
        return None
    match = re.fullmatch(r'([A-Za-z][A-Za-z_0-9]*)\s*(<-|:=|~)\s*(.+)', source)
    if match is None:
        expression = ast.parse(_python_expression(source), mode='eval').body
        if (not isinstance(expression, ast.Call) or not isinstance(expression.func, ast.Name)
                or expression.func.id != 'condition' or len(expression.args) != 1 or expression.keywords):
            raise ValueError('Expected an assignment or condition(expression).')
        validate(expression.args[0])
        return None, 'condition', expression.args[0]
    name, operator, text = match.groups()
    if name == 'condition':
        raise ValueError('condition is reserved for standalone condition(expression) statements.')
    expression = parse_expression(text)
    return name, operator, expression


# Reject executable Python constructs before a model is compiled or any values are sampled.
def parse_expression(source):
    expression = ast.parse(_python_expression(source), mode='eval').body
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
    if isinstance(expression, ast.UnaryOp) and (type(expression.op) in _UNARY_OPERATORS
                                              or isinstance(expression.op, ast.Not)):
        validate(expression.operand)
        return
    if isinstance(expression, ast.Compare):
        if len(expression.ops) != 1 or type(expression.ops[0]) not in _COMPARISONS:
            raise ValueError('Use separate comparisons joined by && or ||, not chained comparisons.')
        validate(expression.left)
        validate(expression.comparators[0])
        return
    if isinstance(expression, ast.BoolOp) and isinstance(expression.op, (ast.And, ast.Or)):
        for value in expression.values:
            validate(value)
        return
    if isinstance(expression, ast.BinOp) and type(expression.op) in _BINARY_OPERATORS:
        validate(expression.left)
        validate(expression.right)
        return
    if isinstance(expression, ast.Subscript):
        validate(expression.value)
        validate(expression.slice)
        return
    if isinstance(expression, ast.Call) and isinstance(expression.func, ast.Name):
        if expression.func.id == 'condition':
            raise ValueError('condition(...) is allowed only as a standalone statement.')
        seen = set()
        for keyword in expression.keywords:
            if keyword.arg is None or keyword.arg in seen:
                raise ValueError('Duplicate keyword arguments and ** expansion are not supported.')
            seen.add(keyword.arg)
            validate(keyword.value)
        for argument in expression.args:
            validate(argument)
        return
    raise ValueError('Expected a literal, variable, vector, scalar arithmetic/comparison/Boolean expression, '
                     'indexing expression, or named constructor.')


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


# Evaluate the validated AST against one replicate's values. Subscripts are always
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
        if isinstance(expression.op, ast.Not):
            value = evaluate(expression.operand, lookup, construct)
            if type(value) is not bool:
                raise TypeError('! requires a Boolean scalar operand.')
            return not value
        return arithmetic(expression.op, evaluate(expression.operand, lookup, construct))
    if isinstance(expression, ast.BoolOp):
        for operand in expression.values:
            value = evaluate(operand, lookup, construct)
            if type(value) is not bool:
                raise TypeError('&& and || require Boolean scalar operands.')
            if isinstance(expression.op, ast.And) and not value:
                return False
            if isinstance(expression.op, ast.Or) and value:
                return True
        return value
    if isinstance(expression, ast.Compare):
        left = evaluate(expression.left, lookup, construct)
        right = evaluate(expression.comparators[0], lookup, construct)
        numeric = type(left) in (int, float) and type(right) in (int, float)
        equality = isinstance(expression.ops[0], (ast.Eq, ast.NotEq))
        if not numeric and not (equality and type(left) is type(right) and type(left) in (bool, str)):
            raise TypeError('Comparisons require numeric scalars; == and != also accept matching Boolean or string scalars.')
        return _COMPARISONS[type(expression.ops[0])](left, right)
    if isinstance(expression, ast.BinOp):
        left = evaluate(expression.left, lookup, construct)
        right = evaluate(expression.right, lookup, construct)
        return arithmetic(expression.op, left, right)
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


if __name__ == "__main__":
    # python -m phylojunction.interface.cmdbox.expressions
    # Inspect the parsed call and its dependencies without constructing or sampling it.
    source = 'a <- sse_rate(name="lambda0",value=rate,event="w_speciation",states=[0,0,0])'
    name, operator, expression = parse_statement(source)
    print(name, operator, ast.dump(expression, indent=2))
    print('Dependencies:', list(variable_names(expression)))
