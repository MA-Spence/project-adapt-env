"""Validate ANA-072 deliverables; --check never changes a recorded artifact."""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import hashlib
import json
import re
from urllib.parse import unquote, urlsplit

import nbformat
import numpy as np
import pandas as pd

from build_notebooks import source_cells

HERE = Path(__file__).resolve().parent
ANA = HERE.parent
ROOT = ANA.parents[1]
TABLES = ANA / 'tables'
OUTPUT = TABLES / 'artifact_validation.json'
STEMS = ['validation', 'analysis', 'architecture', 'difficulty']
FIGURES = ['01_average_fitness_distribution', '02_average_local_behaviour',
           '03_local_genetic_architecture', '04_latent_readout_observation',
           '05_static_learning_and_navigation_proxies',
           '06_resolution_and_neutral_connectivity', '07_joint_search']
MODES = {'none', 'independent_difference', 'shared_stop_sensitivity'}
LAYERS = {('none', 'intrinsic'), ('independent_difference', 'intrinsic'),
          ('independent_difference', 'observed'),
          ('shared_stop_sensitivity', 'intrinsic'),
          ('shared_stop_sensitivity', 'observed')}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalize(source):
    # labproj's paired export adds a header and may add boundary blank lines.
    # Cell content, code indentation, and all nontrailing whitespace remain checked.
    return '\n'.join(line.rstrip() for line in source.replace('\r\n', '\n').splitlines()).strip('\n')


def main(check_only=False):
    errors, checks, checked_files = [], {}, {}

    def file_record(path):
        assert path.is_file() and path.stat().st_size > 0, f'Missing/empty file: {path}'
        name = str(path.relative_to(ROOT))
        checked_files[name] = dict(bytes=path.stat().st_size, sha256=sha(path))
        return path

    def category(name, function):
        try:
            checks[name] = function()
        except Exception as exc:
            errors.append(dict(check=name, error=f'{type(exc).__name__}: {exc}'))
            checks[name] = dict(passed=False)

    def notebook_checks():
        output = []
        required_images = {'validation': 1, 'analysis': 2, 'architecture': 2, 'difficulty': 2}
        for stem in STEMS:
            source_path = file_record(HERE / f'{stem}.py')
            notebook_path = file_record(HERE / f'{stem}.ipynb')
            nb = nbformat.read(notebook_path, as_version=4)
            nbformat.validate(nb)
            paired = source_cells(source_path)
            assert len(paired) == len(nb.cells), f'{stem}: cell count differs'
            for i, ((kind, text), cell) in enumerate(zip(paired, nb.cells)):
                assert kind == cell.cell_type, f'{stem}: cell {i} type differs'
                assert normalize(text) == normalize(cell.source), f'{stem}: cell {i} source differs'
            code = [cell for cell in nb.cells if cell.cell_type == 'code']
            assert code, f'{stem}: no code cells'
            counts = [cell.execution_count for cell in code]
            assert counts == list(range(1, len(code) + 1)), f'{stem}: incomplete/nonsequential execution {counts}'
            outputs = [entry for cell in code for entry in cell.outputs]
            assert not any(entry.output_type == 'error' for entry in outputs), f'{stem}: error output'
            images = sum('image/png' in entry.get('data', {}) for entry in outputs)
            assert images >= required_images[stem], f'{stem}: expected embedded figures'
            output.append(dict(notebook=f'notebooks/{stem}.ipynb', cells=len(nb.cells),
                code_cells=len(code), embedded_pngs=images, source_parity=True,
                sequential_execution=True, error_outputs=0))
        return output

    def figure_checks():
        rows = []
        for stem in FIGURES:
            png, pdf = (file_record(ANA / 'figures' / f'{stem}.{ext}') for ext in ('png', 'pdf'))
            assert png.read_bytes().startswith(b'\x89PNG\r\n\x1a\n'), f'{png.name}: not PNG'
            assert pdf.read_bytes().startswith(b'%PDF-'), f'{pdf.name}: not PDF'
            rows.append(dict(stem=stem, png_bytes=png.stat().st_size, pdf_bytes=pdf.stat().st_size))
        return rows

    def load_table(name):
        return pd.read_csv(file_record(TABLES / name))

    def cohort_check(frame, rows_per_seed):
        synthetic = frame[frame.seed.ne(-1)]
        counts = synthetic.groupby(['seed', 'cohort']).size()
        assert len(counts) == 24 and counts.eq(rows_per_seed).all(), counts.to_dict()
        cohort_seeds = synthetic.groupby('cohort').seed.nunique().to_dict()
        assert cohort_seeds == {'fresh': 8, 'reused': 16}, cohort_seeds
        return cohort_seeds

    def table_checks():
        result = {}
        trials = load_table('validation_search_summary.csv')
        assert len(trials) == 3 and set(trials.observation) == MODES
        assert trials.complete.eq(11712).all() and trials.failed.eq(0).all()
        assert trials.loss_improvement_over_warm_start.eq(0).all()
        assert trials.off_grid_trials.eq(11706).all()
        result['search'] = dict(studies=3, completed_trials=int(trials.complete.sum()), failed_trials=0)
        latent = load_table('validation_latent_replay.csv')
        assert len(latent) == latent.seed.nunique() == 24
        assert latent.maximum_a_difference.eq(0).all() and latent.maximum_g_difference.max() < 1e-13
        average = load_table('average_full_library_by_seed.csv')
        assert len(average) == 120
        assert set(zip(average.observation, average.layer)) == LAYERS
        assert not average.duplicated(['seed', 'observation', 'layer']).any()
        result['average'] = dict(rows=120, cohorts=cohort_check(average, 5))
        replay = load_table('average_reconstruction_checks.csv')
        assert len(replay) == 120 and replay.nodes.eq(160000).all()
        assert int(replay.saved_cube_values.sum()) == 15360
        assert replay.maximum_absolute_difference.max() < 1e-12
        for name in ['average_distribution_summary.csv', 'average_local_summary.csv']:
            table = load_table(name)
            assert set(table.cohort) == {'fresh', 'reused', 'all'}
            assert table.loc[table.cohort.eq('fresh'), 'n'].eq(8).all()
            assert table.loc[table.cohort.eq('reused'), 'n'].eq(16).all()
            assert table.loc[table.cohort.eq('all'), 'n'].eq(24).all()
        architecture = load_table('architecture_seed_metrics.csv')
        observed = architecture[architecture.layer.ne('latent')]
        assert len(observed) == 363 and len(architecture) == 435
        assert set(observed.relative_resolution) == {0., .001, .01}
        assert len(observed[observed.seed.eq(-1)]) == 3
        result['architecture'] = dict(rows=len(architecture), archived_layer_rows=363,
            cohorts=cohort_check(observed, 15))
        edges = load_table('difficulty_edge_seed_metrics.csv')
        assert len(edges) == 121 and edges.nodes.eq(159129).all() and edges.edges.eq(6016433).all()
        assert edges.neighbor_correlation.between(-1, 1).all()
        result['edges'] = dict(rows=121, cohorts=cohort_check(edges, 5), nodes=159129, undirected_edges=6016433)
        graphs = load_table('difficulty_graph_seed_metrics.csv')
        assert len(graphs) == 363 and set(graphs.relative_resolution) == {0., .001, .01}
        result['graphs'] = dict(rows=363, cohorts=cohort_check(graphs, 15))
        for name in ['architecture_cohort_summary.csv', 'difficulty_edge_summary.csv', 'difficulty_graph_summary.csv']:
            table = load_table(name)
            assert len(table) > 0, f'Empty summary: {name}'
        return result

    def scientific_checks():
        names = ['validation_checks.json', 'average_reconstruction.json',
                 'architecture_validation.json', 'difficulty_verification.json',
                 'analysis_environment.json']
        docs = {name: json.loads(file_record(TABLES / name).read_text()) for name in names}
        assert docs['validation_checks.json']['archive']['all_hashes_match'] is True
        assert docs['validation_checks.json']['model_replay']['selected_parameters_equal_RUN180'] is True
        assert docs['average_reconstruction.json']['reconstructed_maps'] == 120
        assert set(docs['architecture_validation.json']['self_checks'].values()) == {'pass'}
        assert docs['difficulty_verification.json']['maps_summarized'] == 121
        assert docs['analysis_environment.json']['python'] == '3.11.15'
        for path in [ROOT / 'data/processed/aggregate-thermodynamic-gnk-calibration/RUN-186.dvc',
                     HERE / 'requirements.txt', HERE / 'environment.lock.txt',
                     HERE / 'build_notebooks.py', HERE / 'execute.py', HERE / 'check_artifacts.py',
                     HERE / 'average_tools.py', HERE / 'architecture_tools.py', HERE / 'difficulty_tools.py']:
            file_record(path)
        return dict(records=names, exact_recorded_python='3.11.15', recorded_checks_present=True)

    def link_checks():
        results = []
        for filename in ['analysis.md', 'README.md']:
            path = file_record(ANA / filename)
            text = path.read_text()
            assert len(text) > 1000, f'{filename}: report still appears to be a template'
            links = []
            for match in re.finditer(r'!?\[[^\]]*\]\(([^)]+)\)', text):
                target = match.group(1).strip()
                target = target[1:-1] if target.startswith('<') and target.endswith('>') else target
                parts = urlsplit(target)
                if parts.scheme or parts.netloc or not parts.path:
                    continue
                relative = unquote(parts.path)
                resolved = (ANA / relative).resolve()
                pending_this_output = resolved == OUTPUT.resolve() and not check_only
                assert resolved.exists() or pending_this_output, f'{filename}: broken link {target}'
                links.append(relative)
            assert links, f'{filename}: no checked local links'
            if filename == 'analysis.md':
                for stem in FIGURES:
                    assert f'figures/{stem}.png' in links, f'Report does not link figure {stem}'
            results.append(dict(document=filename, local_links_checked=len(links),
                                links=links, all_targets_exist=True))
        return results

    for name, function in [('notebooks', notebook_checks), ('figures', figure_checks),
                           ('table_coverage', table_checks), ('scientific_records', scientific_checks),
                           ('report_links', link_checks)]:
        category(name, function)
    if check_only:
        try:
            previous = json.loads(OUTPUT.read_text())
            assert previous['passed'] is True, 'Recorded artifact validation did not pass'
            assert previous['checked_files'] == checked_files, 'Artifact bytes changed since recorded validation; rerun without --check before archiving'
        except Exception as exc:
            errors.append(dict(check='recorded_validation', error=f'{type(exc).__name__}: {exc}'))
    result = dict(analysis='ANA-072', checked_at=datetime.now(timezone.utc).isoformat(),
                  passed=not errors, scope='Artifact/provenance checks, not biological adequacy or task performance.',
                  checks=checks, checked_files=checked_files, errors=errors)
    if not check_only:
        TABLES.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(dict(passed=result['passed'], mode='check-only' if check_only else 'record',
        files=len(checked_files), notebooks=len(checks.get('notebooks', [])),
        figure_pairs=len(checks.get('figures', [])), errors=errors), indent=2))
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Validate without changing artifact_validation.json; require its recorded hashes to match.')
    raise SystemExit(main(parser.parse_args().check))
