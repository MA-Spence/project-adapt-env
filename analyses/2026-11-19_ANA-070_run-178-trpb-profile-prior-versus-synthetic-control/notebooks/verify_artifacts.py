"""Check executed notebooks, embedded figures, and source pairing; write the archive manifest."""
from pathlib import Path
import base64
import hashlib
import json
import re
import nbformat
from PIL import Image

ANA = Path(__file__).resolve().parents[1]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def source_cells(path):
    result = []
    for block in re.split(r'^# %%', path.read_text(), flags=re.MULTILINE)[1:]:
        first, text = block.split('\n', 1)
        text = text.strip('\n')
        kind = 'markdown' if '[markdown]' in first else 'code'
        if kind == 'markdown':
            text = '\n'.join(line[2:] if line.startswith('# ') else line[1:] if line == '#' else line for line in text.splitlines())
        result.append((kind, text))
    return result

checks = {}
embedded = []
for stem, expected in [('analysis', 8), ('validation', 3)]:
    path = ANA / 'notebooks' / f'{stem}.ipynb'
    nb = nbformat.read(path, as_version=4)
    nbformat.validate(nb)
    assert [(cell.cell_type, cell.source) for cell in nb.cells] == source_cells(path.with_suffix('.py'))
    code = [cell for cell in nb.cells if cell.cell_type == 'code']
    assert len(code) == expected
    assert [cell.execution_count for cell in code] == list(range(1, expected + 1))
    outputs = [output for cell in code for output in cell.outputs]
    assert not any(output.output_type == 'error' for output in outputs)
    images = [hashlib.sha256(base64.b64decode(output.data['image/png'])).hexdigest()
              for output in outputs if 'image/png' in output.get('data', {})]
    embedded.extend(images)
    checks[stem] = {'code_cells': len(code), 'executed': True, 'errors': 0,
                    'embedded_pngs': len(images), 'paired_source_matches': True,
                    'sha256': sha(path)}
assert len(embedded) == 5
figures = {}
for path in sorted((ANA / 'figures').glob('*.png')):
    with Image.open(path) as picture:
        picture.verify()
    with Image.open(path) as picture:
        size = picture.size
        assert min(size) > 1000
    pdf = path.with_suffix('.pdf')
    assert pdf.read_bytes().startswith(b'%PDF') and pdf.stat().st_size > 10000
    assert sha(path) in embedded
    figures[path.stem] = {'pixels': size, 'png_sha256': sha(path), 'pdf_sha256': sha(pdf)}
assert len(figures) == 5
# Validate local links in the report and README. The storage sidecar is plain text,
# not a link, and is created after archiving to avoid circular identities.
for name in ['analysis.md', 'README.md']:
    for target in re.findall(r'\]\(([^)]+)\)', (ANA / name).read_text()):
        assert (ANA / target).exists(), (name, target)
validation = {'notebooks': checks, 'figures': figures,
              'local_report_links_exist': True,
              'visual_review': 'All five PNG figures inspected for legibility, labels, target placement and clipping.',
              'scientific_adequacy_gate': False}
(ANA / 'tables/artifact_validation.json').write_text(json.dumps(validation, indent=2) + '\n')
paths = [ANA / 'README.md', ANA / 'analysis.md']
for directory in ['notebooks', 'figures', 'tables']:
    paths.extend(p for p in (ANA / directory).rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix != '.pyc')
manifest = {str(p.relative_to(ANA)): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted(paths)}
(ANA / 'artifact_manifest.json').write_text(json.dumps({'analysis': 'ANA-070', 'files': manifest,
    'exclusions': ['metadata.yaml', 'storage_verification.json', 'artifact_manifest.json', '__pycache__'],
    'source_runs': ['RUN-178', 'RUN-175']}, indent=2) + '\n')
print(json.dumps({'notebooks': checks, 'figures': len(figures), 'manifest_files': len(manifest)}, indent=2))
