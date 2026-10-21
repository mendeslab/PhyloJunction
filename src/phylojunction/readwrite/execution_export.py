"""Export replicate values as separate CSV files, combined CSV, or JSON/MCON."""
import io
import json
import re
from pathlib import Path
import pandas as pd
from phylojunction.data.tree import AnnotatedTree
from phylojunction.pgm.values import elements, path_text, shape, summaries

FORMATS = {'separate-csv': 'Separate files (CSV)', 'combined-csv': 'Combined (CSV)', 'json': 'JSON (MCON)'}


# Reject competing producers before a dictionary assignment can discard an earlier file.
# Existing files from previous exports are separate: normal overwrite behavior is unchanged.
def add_output(outputs, filename, content):
    if filename in outputs:
        raise ValueError(f"Conflicting export filename {filename!r}; rename the competing model variables.")
    outputs[filename] = content


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


# Reuse the established tree serializers, remapping their row ordinals to replicate/path
# coordinates. No temporary DAG nodes or second per-replicate graph are created.
def tree_files(results, write_nex_states=False, bit_format=False, write_stub_counts=False):
    from phylojunction.readwrite.pj_write import prep_data_df, prep_data_filepaths_dfs
    outputs = {}
    for name in results.model.name_node_dict:
        entries = [(i, path_text(path), value) for i, values in enumerate(results.values(name), 1)
                   for path, value in elements(values) if isinstance(value, AnnotatedTree)]
        if not entries:
            continue
        if write_stub_counts:
            # Node labels identify edges within each tree; the union of columns preserves
            # different topologies, with absent labels left empty rather than assigned zero.
            rows = []
            for replicate, path, tree in entries:
                tree.populate_n_stubs_on_branches_of_rec_tr()
                counts = {nd.label: nd.n_stub for nd in tree.tree_reconstructed if nd.label is not None}
                if '@replicate' in counts or '@element' in counts:
                    raise ValueError("Tree labels conflict with stub-table identifier columns.")
                rows.append({'@replicate': replicate, '@element': path, **counts})
            table = pd.DataFrame(rows)
            for column in table.columns[2:]:
                table[column] = table[column].astype('Int64')
            add_output(outputs, name + '_n_stub.tsv', table)
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
                content.insert(0, '@replicate', [entries[i][0] for i in indices])
            match = re.search(r'_sample(\d+)_repl1', filename)
            if match:
                execution, path, _ = entries[int(match.group(1)) - 1]
                filename = filename[:match.start()] + f'_replicate{execution}{path}' + filename[match.end():]
            add_output(outputs, filename, content)
    return outputs


# Stage every output in memory before opening files. Combined CSV requires identical
# complete shapes; MCON records preserve lists and omit synthetic replicate variables.
def write_results(directory, results, prefix='', output_format='separate-csv',
                  write_nex_states=False, bit_format=False, write_stub_counts=False):
    if output_format not in FORMATS:
        raise ValueError(f'Unknown output format: {output_format}')
    columns = export_values(results)
    outputs = {}
    if output_format == 'combined-csv':
        rows = []
        for name, values in columns.items():
            if any(shape(v) != shape(values[0]) for v in values[1:]):
                raise ValueError(f'{name} changes vector shape between replicates; use JSON (MCON).')
        for i in range(results.n_executions):
            row = {'@replicate': i + 1}
            for name, values in columns.items():
                row.update({name + path_text(path): csv_value(v) for path, v in elements(values[i])})
            rows.append(row)
        add_output(outputs, 'combined.csv', pd.DataFrame(rows))
    elif output_format == 'json':
        header = {'format': 'MCON', 'version': '0.2', 'nested': False, 'atomic': False}
        records = [header] + [{name: values[i] for name, values in columns.items()}
                              for i in range(results.n_executions)]
        add_output(outputs, 'samples.mcon', ''.join(json.dumps(row, allow_nan=False) + '\n' for row in records))
    else:
        for name, values in columns.items():
            rows = [{'@replicate': i, '@element': path_text(path), name: csv_value(value)}
                    for i, vector in enumerate(values, 1) for path, value in elements(vector)]
            add_output(outputs, name + '.csv', pd.DataFrame(rows))
    for filename, content in tree_files(results, write_nex_states, bit_format, write_stub_counts).items():
        add_output(outputs, filename, content)
    for name in columns:
        stats = summaries(results.values(name))
        if stats:
            add_output(outputs, name + '_stats_summary.csv', pd.DataFrame(stats))
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
    if isinstance(names, str):
        raise TypeError("Map export requires a list of variable names, not a string.")
    names = list(names)
    if not names or any(not isinstance(name, str) or not name for name in names):
        raise ValueError("Select at least one tree variable; empty names are not allowed.")
    if len(set(names)) != len(names):
        raise ValueError("Map export variable names must be distinct.")
    records = []
    record_names = set()
    for name in names:
        if name not in results.model.name_node_dict:
            raise ValueError(f"Unknown map export variable {name!r}.")
        for i, vector in enumerate(results.values(name), 1):
            for path, value in elements(vector):
                if not isinstance(value, AnnotatedTree):
                    raise TypeError(f'{name}{path_text(path)}, replicate {i}, is not a tree.')
                key = f'{name}_replicate{i}{path_text(path)}'
                if key in record_names:
                    raise ValueError(f'Conflicting map export identifier {key!r}.')
                record_names.add(key)
                records.append((key, [value], 1, 1))
    frames, _ = prep_trees_rb_smap_dfs(records, [r[0] for r in records], attribute)
    outputs = {}
    for key, tables in frames.items():
        filename = ((prefix + '_') if prefix else '') + key + '_smaps.tsv'
        add_output(outputs, filename, tables[0])
    destination = Path(directory)
    destination.mkdir(parents=True, exist_ok=True)
    for filename, table in outputs.items():
        table.to_csv(destination / filename, sep='\t', index=False)
