"""Export execution values as separate CSV files, combined CSV, or JSON/MCON."""
import io
import json
import re
from pathlib import Path
import pandas as pd
from phylojunction.data.tree import AnnotatedTree
from phylojunction.pgm.values import elements, path_text, shape, summaries

FORMATS = {'separate-csv': 'Separate files (CSV)', 'combined-csv': 'Combined (CSV)', 'json': 'JSON (MCON)'}


# Domain construction objects are model definitions, not observations. Trees use annotated
# Newick in tabular/JSON fields; dedicated files retain tree statistics and node observations.
def json_value(value):
    if isinstance(value, list):
        return [json_value(item) for item in value]
    if isinstance(value, AnnotatedTree):
        return value.tree.as_string(schema='newick', suppress_annotations=False).strip()
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise TypeError(f'{type(value).__name__} is a model construction object, not an exportable value.')


# Collect before writing so unsupported shapes cannot leave a misleading partial CSV.
def export_values(results):
    columns = {}
    for name in results.model.name_node_dict:
        values = results.values(name)
        try:
            encoded = [json_value(value) for value in values]
        except TypeError:
            continue
        columns[name] = encoded
    return columns


# JSON is used only for list-valued CSV fields, including empty vectors; csv handles quoting.
def csv_value(value):
    return json.dumps(value, allow_nan=False) if isinstance(value, list) else value


# Reuse the established tree serializers, remapping their row ordinals to execution/path
# coordinates. No temporary DAG nodes or second per-execution graph are created.
def tree_files(results, write_nex_states=False, bit_format=False):
    from phylojunction.readwrite.pj_write import prep_data_df, prep_data_filepaths_dfs
    outputs = {}
    for name in results.model.name_node_dict:
        entries = [(i, path_text(path), value) for i, values in enumerate(results.values(name), 1)
                   for path, value in elements(values) if isinstance(value, AnnotatedTree)]
        if not entries:
            continue
        records = [(name, [v for _, _, v in entries], 1, True, True, len(entries))]
        scalar, trees = prep_data_df(records, write_nex_states, bit_format)
        filenames, contents = prep_data_filepaths_dfs(scalar, trees)
        for filename, content in zip(filenames, contents):
            if filename.endswith(('_traits.tsv', '_fossils.tsv')):
                content = pd.read_csv(io.StringIO(content), sep='\t')
            if isinstance(content, pd.DataFrame) and 'sample' in content:
                indices = [int(i) - 1 for i in content['sample']]
                content = content.drop(columns=['sample', 'replicate'], errors='ignore')
                content.insert(0, '@element', [entries[i][1] for i in indices])
                content.insert(0, '@execution', [entries[i][0] for i in indices])
            match = re.search(r'_sample(\d+)_repl1', filename)
            if match:
                execution, path, _ = entries[int(match.group(1)) - 1]
                filename = filename[:match.start()] + f'_execution{execution}{path}' + filename[match.end():]
            outputs[filename] = content
    return outputs


# Stage every output in memory before opening files. Combined CSV requires identical
# complete shapes; MCON records preserve lists and omit synthetic execution variables.
def write_results(directory, results, prefix='', output_format='separate-csv',
                  write_nex_states=False, bit_format=False):
    if output_format not in FORMATS:
        raise ValueError(f'Unknown output format: {output_format}')
    columns = export_values(results)
    outputs = {}
    if output_format == 'combined-csv':
        rows = []
        for name, values in columns.items():
            if any(shape(v) != shape(values[0]) for v in values[1:]):
                raise ValueError(f'{name} changes vector shape between executions; use JSON (MCON).')
        for i in range(results.n_executions):
            row = {'@execution': i + 1}
            for name, values in columns.items():
                row.update({name + path_text(path): csv_value(v) for path, v in elements(values[i])})
            rows.append(row)
        outputs['combined.csv'] = pd.DataFrame(rows)
    elif output_format == 'json':
        header = {'format': 'MCON', 'version': '0.2', 'nested': False, 'atomic': False}
        records = [header] + [{name: values[i] for name, values in columns.items()}
                              for i in range(results.n_executions)]
        outputs['samples.mcon'] = ''.join(json.dumps(row, allow_nan=False) + '\n' for row in records)
    else:
        for name, values in columns.items():
            rows = [{'@execution': i, '@element': path_text(path), name: csv_value(value)}
                    for i, vector in enumerate(values, 1) for path, value in elements(vector)]
            outputs[name + '.csv'] = pd.DataFrame(rows)
    outputs.update(tree_files(results, write_nex_states, bit_format))
    for name in columns:
        stats = summaries(results.values(name))
        if stats:
            outputs[name + '_stats_summary.csv'] = pd.DataFrame(stats)
    destination = Path(directory)
    destination.mkdir(parents=True, exist_ok=True)
    for filename, content in outputs.items():
        path = destination / (((prefix + '_') if prefix else '') + filename)
        if isinstance(content, pd.DataFrame):
            content.to_csv(path, index=False, sep='\t' if path.suffix == '.tsv' else ',')
        else:
            path.write_text(content)
    return [destination / (((prefix + '_') if prefix else '') + name) for name in outputs]


# Preserve the existing stochastic-map serializer and identify each selected tree by path.
def write_maps(directory, results, names, attribute, prefix=''):
    from phylojunction.readwrite.pj_write import prep_trees_rb_smap_dfs
    records = []
    for name in names:
        for i, vector in enumerate(results.values(name), 1):
            for path, value in elements(vector):
                if not isinstance(value, AnnotatedTree):
                    raise TypeError(f'{name}{path_text(path)} is not a tree.')
                key = f'{name}_execution{i}{path_text(path)}'
                records.append((key, [value], 1, 1))
    frames, _ = prep_trees_rb_smap_dfs(records, [r[0] for r in records], attribute)
    destination = Path(directory)
    destination.mkdir(parents=True, exist_ok=True)
    for key, tables in frames.items():
        filename = ((prefix + '_') if prefix else '') + key + '_smaps.tsv'
        tables[0].to_csv(destination / filename, sep='\t', index=False)
