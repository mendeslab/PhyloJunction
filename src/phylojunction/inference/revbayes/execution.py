"""RevBayes scripts from shared definitions and execution-specific observations."""
import ast
import inspect
import json
from pathlib import Path
from phylojunction.pgm.pgm import DistrForSampling
from phylojunction.distribution.dn_iid import DnIID
from phylojunction.distribution import dn_parametric as dn
from phylojunction.interface.grammar.constructors import CONSTRUCTORS, flag
from phylojunction.interface.cmdbox.expressions import evaluate


# Empty literals have no inferable element type in Rev; reject rather than emit invalid [].
def literal(value):
    if isinstance(value, list):
        if not value:
            raise ValueError('RevBayes export requires nonempty vectors; CSV/MCON supports empty vectors.')
        return '[' + ', '.join(literal(item) for item in value) + ']'
    if isinstance(value, (str, int, float, bool)):
        return json.dumps(value, allow_nan=False)
    raise ValueError(f'RevBayes export does not support {type(value).__name__} values.')


# Emit separate scripts for executions. Definitions retain symbolic dependencies; only
# observations and fixed loop dimensions come from an execution's realized values.
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
    helper = 'pj'
    while any(name.startswith(helper) for name in occupied):
        helper += 'X'
    scripts, logging = [], []
    stem = prefix + '_' if prefix else ''
    directories = [str(Path(root)), str(Path(root) / (stem + 'rev_scripts')),
                   str(Path(root) / (stem + 'inference_results'))]
    for execution, row in enumerate(results._values):
        lines, moves, stochastic = [], [], []

        # Resolve aliases to distribution expressions while retaining their defining scope.
        def resolve(expression, owner):
            while isinstance(expression, ast.Name):
                target = definitions[owner.bindings[expression.id]]
                if not isinstance(row[target.definition_id], DistrForSampling):
                    break
                expression, owner = target.expression, target
            return expression, owner

        # Render the restricted expression language as Rev, converting parameterizations.
        def render(expression, owner):
            if isinstance(expression, ast.Constant):
                return literal(expression.value)
            if isinstance(expression, ast.Name):
                if expression.id in ('true', 'false'):
                    return expression.id
                return names[owner.bindings[expression.id]]
            if isinstance(expression, ast.List):
                if not expression.elts:
                    return literal([])
                return '[' + ', '.join(render(item, owner) for item in expression.elts) + ']'
            if isinstance(expression, ast.UnaryOp):
                return ('-' if isinstance(expression.op, ast.USub) else '+') + render(expression.operand, owner)
            if isinstance(expression, ast.Subscript):
                return render(expression.value, owner) + '[' + render(expression.slice, owner) + ']'
            if not isinstance(expression, ast.Call):
                raise ValueError('Unsupported Rev expression.')
            constructor = expression.func.id
            supported = {'normal', 'lognormal', 'unif', 'gamma', 'exponential', 'beta'}
            if constructor not in supported:
                raise ValueError(f'RevBayes export does not support {constructor}.')
            bound = inspect.signature(CONSTRUCTORS[constructor]).bind(
                *expression.args, **{kw.arg: kw.value for kw in expression.keywords})
            bound.apply_defaults()
            args = bound.arguments

            def text(name):
                value = args[name]
                return render(value, owner) if isinstance(value, ast.AST) else literal(value)

            # Parameterization switches are configuration values evaluated in this execution.
            def boolean(name):
                value = args[name]
                if isinstance(value, ast.AST):
                    value = evaluate(value, lambda n: row[owner.bindings[n]], None)
                return flag(value)

            if constructor == 'normal':
                return f'dnNormal(mean={text("mean")}, sd={text("sd")})'
            if constructor == 'beta':
                return f'dnBeta(alpha={text("alpha")}, beta={text("beta")})'
            if constructor == 'unif':
                return f'dnUniform(lower={text("min")}, upper={text("max")})'
            if constructor == 'lognormal':
                mean = text('meanlog') if boolean('log_space') else f'ln({text("meanlog")})'
                return f'dnLognormal(mean={mean}, sd={text("sdlog")})'
            parameter = 'rate' if constructor == 'exponential' else 'scale'
            rate = text(parameter) if boolean('rate_parameterization') else f'1.0/({text(parameter)})'
            return f'dnExponential(lambda={rate})' if constructor == 'exponential' else \
                f'dnGamma(shape={text("shape")}, rate={rate})'

        # Descend explicit IID dimensions; leaf statements remain independent stochastic
        # nodes in Rev, while PhyloJunction retains one vector-valued DAG definition.
        def sample(target, expression, owner, distribution, observed, clamped):
            expression, owner = resolve(expression, owner)
            if isinstance(distribution, DnIID):
                if not isinstance(expression, ast.Call) or expression.func.id != 'iid':
                    raise ValueError('Cannot resolve IID distribution expression for Rev export.')
                bound = inspect.signature(DnIID).bind(*expression.args,
                    **{kw.arg: kw.value for kw in expression.keywords})
                count_expression = bound.arguments['n']
                from phylojunction.interface.cmdbox.expressions import variable_names
                pending = [owner.bindings[n] for n in variable_names(count_expression)]
                while pending:
                    dependency = definitions[pending.pop()]
                    if dependency.is_sampled:
                        raise ValueError('Rev export does not support stochastic IID dimensions.')
                    pending.extend(dependency.bindings.values())
                if distribution.n == 0:
                    literal([])
                for i, value in enumerate(observed, 1):
                    sample(f'{target}[{i}]', bound.arguments['distribution'], owner,
                           distribution.distribution, value, clamped)
                return
            if not isinstance(distribution, (dn.DnNormal, dn.DnLogNormal, dn.DnUnif,
                                             dn.DnExponential, dn.DnGamma, dn.DnBeta)):
                raise ValueError(f'Rev export does not support {type(distribution).__name__}.')
            lines.append(f'{target} ~ {render(expression, owner)}')
            stochastic.append(target)
            if clamped:
                lines.append(f'{target}.clamp({literal(observed)})')
            else:
                move = 'mvScale' if isinstance(distribution, (dn.DnLogNormal, dn.DnExponential, dn.DnGamma)) else 'mvSlide'
                tuning = 'lambda=1.0' if move == 'mvScale' else 'delta=1.0'
                moves.append(f'{helper}Moves[{len(moves)+1}] = {move}({target}, {tuning}, weight=1.0)')

        for node in definitions:
            value = row[node.definition_id]
            name = names[node.definition_id]
            if node.is_sampled:
                distribution = results._sampled_from[execution][node.definition_id]
                sample(name, node.expression, node, distribution, value, node.is_clamped)
            elif isinstance(value, DistrForSampling):
                # Resolve distribution aliases at sampling sites rather than constructing
                # unsupported Rev vectors of distributions or a fictional dnIID primitive.
                continue
            else:
                lines.append(f'{name} {node.operator if node.operator == ":=" else "<-"} {render(node.expression, node)}')
        if not stochastic or not moves:
            raise ValueError('Rev inference export needs at least one unclamped numeric stochastic value.')
        scripts.append('\n'.join(lines))
        log_path = str(Path(directories[2]) / f'{stem}run{execution+1}.log')
        logging.append('\n'.join(moves + [
            f'{helper}Model = model(' + ', '.join(stochastic) + ')',
            f'{helper}Monitors[1] = mnModel(filename={literal(log_path)}, printgen=1)',
            f'{helper}Mcmc = mcmc({helper}Model, {helper}Moves, {helper}Monitors)',
            f'{helper}Mcmc.run({mcmc_chain_length})', 'quit()']))
    return scripts, logging, directories
