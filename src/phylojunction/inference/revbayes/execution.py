"""RevBayes scripts from shared definitions and replicate-specific observations."""
import ast
import json
import inspect
from dataclasses import dataclass, replace
from pathlib import Path
from phylojunction.pgm.pgm import DistrForSampling
from phylojunction.distribution.dn_iid import DnIID
from phylojunction.distribution import dn_parametric as dn
from phylojunction.interface.grammar.constructors import CONSTRUCTORS, bind_call, construct, flag
from phylojunction.interface.cmdbox.expressions import binding_names, comprehension


# Empty literals have no inferable element type in Rev. Typed empty draws are handled
# separately through dnIID, rather than making up an element type for arbitrary [].
def literal(value):
    if isinstance(value, (list, tuple)):
        if not value:
            raise ValueError('RevBayes export cannot translate an untyped [] literal.')
        return '[' + ', '.join(literal(item) for item in value) + ']'
    if isinstance(value, (str, int, float, bool)):
        return json.dumps(value, allow_nan=False)
    raise ValueError(f'RevBayes export does not support {type(value).__name__} values.')


@dataclass(frozen=True)
class _RevExpression:
    """A symbolic expression and its fixed collection structure, local to the exporter.

    value is an replicate exemplar, never a substitute for a random parameter expression.
    fixed means the value is independent of stochastic nodes. items records list structure;
    inner retains an IID component expression because Rev distribution objects are opaque.
    """

    code: str
    value: object
    fixed: bool
    items: tuple | None = None
    inner: object = None

    # Naming a list preserves its structure but references its entries through the alias.
    def named(self, code):
        return replace(self, code=code, items=None if self.items is None else
                       tuple(item.named(f'{code}[{i}]') for i, item in enumerate(self.items, 1)))


# Emit separate scripts for replicates. Definitions retain symbolic dependencies; only
# observations and proven fixed dimensions come from a replicate's realized values.
def execution_to_rev(results, root, mcmc_chain_length=1000, prefix=''):
    if type(mcmc_chain_length) is not int or mcmc_chain_length < 1:
        raise ValueError('MCMC chain length must be a positive integer.')
    definitions = results.model.definitions
    names = {}
    occupied = set(results.model.name_node_dict)
    for node in definitions:
        name = node.node_name
        if results.model.name_node_dict[name] is not node:
            name += '__definition' + str(node.definition_id + 1)
            while name in occupied:
                name += '_'
        occupied.add(name)
        names[node.definition_id] = name
        occupied.update(item.id for item in ast.walk(node.expression) if isinstance(item, ast.Name))
    helper = 'pjInternal'
    while any(name.startswith(helper) for name in occupied):
        helper += 'X'
    scripts, logging = [], []
    stem = prefix + '_' if prefix else ''
    directories = [str(Path(root)), str(Path(root) / (stem + 'rev_scripts')),
                   str(Path(root) / (stem + 'inference_results'))]
    for execution, row in enumerate(results._values):
        lines, moves, stochastic = [], [], []
        translated = {}
        counter = 0

        # All generated names share one deterministic allocator, independent of Rev scope.
        def fresh():
            nonlocal counter
            counter += 1
            return f'{helper}{counter}'

        moves_name, model_name, monitors_name, mcmc_name = (fresh() for _ in range(4))

        # Collections have fixed structure even when their elements depend on random values.
        # Constructing a literal collection does not evaluate its elements a second time.
        def collection(items):
            if not items:
                raise ValueError('Rev export cannot infer a type for an empty collection; typed IID draws are supported.')
            return _RevExpression('[' + ', '.join(item.code for item in items) + ']',
                                  [item.value for item in items], all(item.fixed for item in items), tuple(items))

        # Convert parameterizations without freezing stochastic switches or dimensions.
        def distribution(name, arguments):
            def text(key):
                return arguments[key].code

            def boolean(key):
                if not arguments[key].fixed:
                    raise ValueError('Rev export requires fixed distribution parameterization switches.')
                return flag(arguments[key].value)

            values = {key: item.value for key, item in arguments.items()}
            value = construct(name, **values)
            fixed = all(item.fixed for item in arguments.values())
            inner = None
            if name == 'normal':
                code = f'dnNormal(mean={text("mean")}, sd={text("sd")})'
            elif name == 'cauchy':
                code = f'dnCauchy(location={text("location")}, scale={text("scale")})'
            elif name == 'beta':
                code = f'dnBeta(alpha={text("alpha")}, beta={text("beta")})'
            elif name == 'unif':
                code = f'dnUniform(lower={text("min")}, upper={text("max")})'
            elif name == 'lognormal':
                mean = text('meanlog') if boolean('log_space') else f'ln({text("meanlog")})'
                code = f'dnLognormal(mean={mean}, sd={text("sdlog")})'
            elif name in ('exponential', 'gamma'):
                parameter = 'rate' if name == 'exponential' else 'scale'
                rate = text(parameter) if boolean('rate_parameterization') else f'1.0/({text(parameter)})'
                code = (f'dnExponential(lambda={rate})' if name == 'exponential' else
                        f'dnGamma(shape={text("shape")}, rate={rate})')
            elif name == 'iid':
                if not arguments['n'].fixed:
                    raise ValueError('Rev export does not support stochastic IID dimensions.')
                inner = arguments['distribution']
                code = f'dnIID({text("n")},{inner.code})'
            else:
                raise ValueError(f'RevBayes export does not support {name}.')
            return _RevExpression(code, value, fixed, inner=inner)

        # Translate expressions using lexical bindings, not textual replacement. Fixed-size
        # comprehensions construct lists of symbolic entries; random values stay references.
        # This specializes known structure without specializing any random parameter value.
        def render(expression, owner, local=None):
            local = {} if local is None else local
            if isinstance(expression, ast.Constant):
                return _RevExpression(literal(expression.value), expression.value, True)
            if isinstance(expression, ast.Name):
                if expression.id in ('true', 'false'):
                    return _RevExpression(expression.id, expression.id == 'true', True)
                if expression.id in local:
                    return local[expression.id]
                return translated[owner.bindings[expression.id]]
            if isinstance(expression, ast.List):
                if not expression.elts:
                    literal([])
                return collection([render(item, owner, local) for item in expression.elts])
            if isinstance(expression, ast.UnaryOp):
                item = render(expression.operand, owner, local)
                negative = isinstance(expression.op, ast.USub)
                return _RevExpression(('-' if negative else '+') + item.code,
                                      -item.value if negative else item.value, item.fixed)
            if isinstance(expression, ast.Subscript):
                vector = render(expression.value, owner, local)
                index = render(expression.slice, owner, local)
                if vector.items is None:
                    raise ValueError('Cannot index an unknown collection structure in Rev export.')
                item = vector.items[index.value - 1]
                if not vector.code:
                    # Zip records have symbolic fields but no Rev mixed-type container.
                    # A fixed selection references that field directly, just as unpacking does.
                    if not index.fixed:
                        raise ValueError('Rev export requires fixed indices into zip records.')
                    return item
                if not index.fixed:
                    # Numeric indexing can remain symbolic. Choosing a distribution or a
                    # subcollection would require proving more than the realized selection.
                    if isinstance(item.value, DistrForSampling) or item.items is not None:
                        raise ValueError('Rev export requires fixed indices for distributions and subcollections.')
                    item = replace(item, fixed=False)
                return item.named(f'{vector.code}[{index.code}]')
            if isinstance(expression, ast.ListComp):
                generator = comprehension(expression)
                iterable = render(generator.iter, owner, local)
                if iterable.items is None:
                    raise ValueError('Rev export requires a fixed comprehension structure.')
                names_local = binding_names(generator.target)
                entries = []
                for item in iterable.items:
                    if isinstance(generator.target, ast.Tuple):
                        if item.items is None or len(item.items) != len(names_local):
                            raise ValueError('Comprehension unpacking has the wrong number of entries.')
                        bindings = dict(zip(names_local, item.items))
                    else:
                        bindings = {names_local[0]: item}
                    entries.append(render(expression.elt, owner, {**local, **bindings}))
                return collection(entries)
            if not isinstance(expression, ast.Call):
                raise ValueError('Unsupported Rev expression.')
            name = expression.func.id
            args = [render(arg, owner, local) for arg in expression.args]
            kwargs = {kw.arg: render(kw.value, owner, local) for kw in expression.keywords}
            if name == 'zip_strict':
                if kwargs or len(args) < 2 or any(arg.items is None for arg in args):
                    raise ValueError('zip_strict requires at least two list arguments.')
                if len({len(arg.items) for arg in args}) != 1:
                    raise ValueError('zip_strict requires equal lengths.')
                # Zip records are compiler bindings, not Rev mixed-type vectors. Their
                # fields are rendered independently when unpacked by a comprehension.
                records = tuple(_RevExpression('', tuple(item.value for item in entries),
                                                all(item.fixed for item in entries), tuple(entries))
                                for entries in zip(*(arg.items for arg in args)))
                return _RevExpression('', [item.value for item in records],
                                      all(arg.fixed for arg in args), records)
            if name not in ('normal', 'lognormal', 'unif', 'gamma', 'exponential', 'beta', 'cauchy', 'iid'):
                raise ValueError(f'RevBayes export does not support {name}.')
            bound, traversed = bind_call(name, [item.value for item in args],
                                         {key: item.value for key, item in kwargs.items()})
            # Reuse Python binding for symbolic arguments; defaults have no dependencies.
            symbolic = inspect.signature(CONSTRUCTORS[name][0]).bind(*args, **kwargs)
            arguments = {key: _RevExpression(literal(value), value, True)
                         for key, value in bound.arguments.items() if key not in symbolic.arguments}
            arguments.update(symbolic.arguments)
            if not traversed:
                return distribution(name, arguments)
            return collection([distribution(name, {key: item.items[i] if key in traversed else item
                                                    for key, item in arguments.items()})
                               for i in range(len(bound.arguments[traversed[0]]))])

        # Output structure comes from the distribution specification, never from assuming
        # a realized variable-length draw has a fixed shape. Supported scalar families and
        # fixed-count IID are sufficient to prove the shapes accepted by this exporter.
        def sampled_value(spec, observed):
            if spec.items is not None:
                children = tuple(sampled_value(item, value) for item, value in zip(spec.items, observed))
            elif isinstance(spec.value, DnIID):
                children = tuple(sampled_value(spec.inner, value) for value in observed)
            else:
                children = None
            return _RevExpression('', observed, False, children)

        # Register moves and observations using fixed output paths. Keeping this separate
        # from sampling loops permits different scalar families in one distribution list.
        def register(target, spec, observed, clamped):
            if spec.items is not None:
                for i, (item, value) in enumerate(zip(spec.items, observed), 1):
                    register(f'{target}[{i}]', item, value, clamped)
                return
            if isinstance(spec.value, DnIID) and spec.value.n:
                for i, value in enumerate(observed, 1):
                    register(f'{target}[{i}]', spec.inner, value, clamped)
                return
            if isinstance(spec.value, DnIID):
                # A zero-count IID node has a typed empty value and no movable elements.
                stochastic.append(target)
                if clamped:
                    # A zero-count draw is deterministically empty; its own typed value
                    # supplies the observation without an untyped [] literal.
                    lines.append(f'{target}.clamp({target})')
                return
            if not isinstance(spec.value, (dn.DnNormal, dn.DnLogNormal, dn.DnUnif,
                                           dn.DnExponential, dn.DnGamma, dn.DnBeta, dn.DnCauchy)):
                raise ValueError(f'Rev export does not support {type(spec.value).__name__}.')
            stochastic.append(target)
            if clamped:
                lines.append(f'{target}.clamp({literal(observed)})')
            else:
                move = 'mvScale' if isinstance(spec.value, (dn.DnLogNormal, dn.DnExponential, dn.DnGamma)) else 'mvSlide'
                tuning = 'lambda=1.0' if move == 'mvScale' else 'delta=1.0'
                moves.append(f'{moves_name}[{len(moves)+1}] = {move}({target}, {tuning}, weight=1.0)')

        # A uniform sampling recipe allows one loop despite unequal list lengths. IID
        # parameters are opaque inside a Rev object, so mixed recipes use explicit paths.
        def recipe(spec):
            if spec.items is not None:
                recipes = [recipe(item) for item in spec.items]
                return ('list', recipes[0]) if recipes and all(r == recipes[0] for r in recipes) else None
            if isinstance(spec.value, DnIID):
                return ('iid', spec.value.n, spec.inner.code, recipe(spec.inner))
            return ('scalar',)

        # Literal list syntax is traversed directly. Other collections use loops when their
        # recipes agree; heterogeneous IID recipes retain explicit element paths instead.
        def sample(target, spec, syntax=None, indent=''):
            if spec.items is not None:
                if isinstance(syntax, ast.List) or recipe(spec) is None:
                    for i, item in enumerate(spec.items, 1):
                        child = syntax.elts[i-1] if isinstance(syntax, ast.List) else None
                        sample(f'{target}[{i}]', item, child, indent)
                    return
                if not spec.items:
                    raise ValueError('Rev export cannot construct an untyped empty sampled list.')
                # Rev cannot call .size() directly on a list literal. Materialize a
                # generated collection once, including inside a source literal list.
                if spec.code.startswith('['):
                    temporary = fresh()
                    lines.append(f'{indent}{temporary} = {spec.code}')
                    spec = spec.named(temporary)
                index = fresh()
                lines.append(f'{indent}if ({spec.code}.size() > 0) {{')
                lines.append(f'{indent}    for ({index} in 1:{spec.code}.size()) {{')
                sample(f'{target}[{index}]', spec.items[0].named(f'{spec.code}[{index}]'), indent=indent+'        ')
                lines.extend([indent+'    }', indent+'}'])
            elif isinstance(spec.value, DnIID) and spec.value.n:
                # IID preserves one component distribution; no need to recover its source AST.
                index = fresh()
                lines.append(f'{indent}for ({index} in 1:{spec.value.n}) {{')
                sample(f'{target}[{index}]', spec.inner, indent=indent+'    ')
                lines.append(indent+'}')
            else:
                lines.append(f'{indent}{target} ~ {spec.code}')

        for node in definitions:
            value = row[node.definition_id]
            name = names[node.definition_id]
            spec = render(node.expression, node)
            if node.is_sampled:
                # Nonliteral collections are evaluated once before traversing. Literal
                # constructors are pure and can be sampled directly without a temporary.
                if spec.items is not None and not isinstance(node.expression, ast.List):
                    temporary = fresh()
                    lines.append(f'{temporary} = {spec.code}')
                    spec = spec.named(temporary)
                sample(name, spec, node.expression)
                register(name, spec, value, node.is_clamped)
                translated[node.definition_id] = sampled_value(spec, value).named(name)
            else:
                if not spec.code:
                    # A zip alias is a structural record binding, expanded when consumed.
                    translated[node.definition_id] = spec
                    continue
                contains_distribution = isinstance(spec.value, DistrForSampling)
                pending = list(spec.items or ())
                while pending:
                    item = pending.pop()
                    contains_distribution |= isinstance(item.value, DistrForSampling)
                    pending.extend(item.items or ())
                operator = '=' if contains_distribution else (':=' if node.operator == ':=' else '<-')
                lines.append(f'{name} {operator} {spec.code}')
                translated[node.definition_id] = spec.named(name)
        if not stochastic or not moves:
            raise ValueError('Rev inference export needs at least one unclamped numeric stochastic value.')
        scripts.append('\n'.join(lines))
        log_path = str(Path(directories[2]) / f'{stem}run{execution+1}.log')
        logging.append('\n'.join(moves + [
            f'{model_name} = model(' + ', '.join(stochastic) + ')',
            f'{monitors_name}[1] = mnModel(filename={literal(log_path)}, printgen=1)',
            f'{mcmc_name} = mcmc({model_name}, {moves_name}, {monitors_name})',
            f'{mcmc_name}.run({mcmc_chain_length})', 'quit()']))
    return scripts, logging, directories
