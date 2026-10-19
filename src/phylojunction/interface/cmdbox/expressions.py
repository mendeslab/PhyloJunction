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
def variable_names(expression):
    if isinstance(expression, ast.Name):
        if expression.id not in ('true', 'false'):
            yield expression.id
    elif isinstance(expression, ast.Call):
        for argument in expression.args:
            yield from variable_names(argument)
        for keyword in expression.keywords:
            yield from variable_names(keyword.value)
    else:
        for child in ast.iter_child_nodes(expression):
            yield from variable_names(child)
