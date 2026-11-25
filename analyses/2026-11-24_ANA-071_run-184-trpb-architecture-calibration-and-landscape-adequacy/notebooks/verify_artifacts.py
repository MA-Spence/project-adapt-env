"""Check executed notebooks, paired sources, embedded figures and report links."""
from pathlib import Path
import base64,hashlib,json,re
import nbformat
from PIL import Image
from build_notebooks import source_cells
ANA=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
checks={};embedded=[]
for stem,n in [('validation',5),('analysis',4),('architecture',5)]:
    path=ANA/'notebooks'/f'{stem}.ipynb'
    nb=nbformat.read(path,as_version=4);nbformat.validate(nb)
    assert [(c.cell_type,c.source) for c in nb.cells]==source_cells(path.with_suffix('.py'))
    code=[c for c in nb.cells if c.cell_type=='code']
    assert len(code)==n and [c.execution_count for c in code]==list(range(1,n+1))
    outputs=[o for c in code for o in c.outputs]
    assert not any(o.output_type=='error' for o in outputs)
    images=[hashlib.sha256(base64.b64decode(o.data['image/png'])).hexdigest()
        for o in outputs if 'image/png' in o.get('data',{})]
    embedded.extend(images)
    checks[stem]={'executed_code_cells':n,'errors':0,'paired_source_matches':True,
        'embedded_pngs':len(images),'sha256':sha(path)}
assert len(embedded)==5
figures={}
for path in sorted((ANA/'figures').glob('*.png')):
    with Image.open(path) as im:im.verify()
    with Image.open(path) as im:size=im.size
    assert min(size)>1000
    pdf=path.with_suffix('.pdf')
    assert pdf.read_bytes().startswith(b'%PDF') and pdf.stat().st_size>10000
    assert sha(path) in embedded
    figures[path.stem]={'pixels':size,'png_sha256':sha(path),'pdf_sha256':sha(pdf)}
assert len(figures)==5
for name in ['analysis.md','README.md']:
    for target in re.findall(r'\]\(([^)]+)\)',(ANA/name).read_text()):
        if '://' not in target:assert (ANA/target).exists(),(name,target)
assert json.loads((ANA/'tables/validation_checks.json').read_text())['status']=='passed'
validation={'notebooks':checks,'figures':figures,'local_report_links_exist':True,
    'visual_review':'All five rendered PNGs inspected for labels, readability, empirical references, all seeds and clipping.',
    'adequacy_gate':'No registered numerical pass threshold; report gives a scoped qualitative assessment.'}
(ANA/'tables/artifact_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
paths=[ANA/'README.md',ANA/'analysis.md',ANA/'metadata.yaml']
for directory in ['notebooks','figures','tables']:
    paths.extend(p for p in (ANA/directory).rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc')
manifest={str(p.relative_to(ANA)):{'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(paths)}
(ANA/'artifact_manifest.json').write_text(json.dumps({'analysis':'ANA-071','files':manifest},indent=2)+'\n')
print('Verified 3 executed notebooks, 14 code cells, 5 embedded PNG/PDF pairs and report links.')
