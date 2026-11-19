"""Execute ANA-068 with this Python interpreter and a temporary kernel spec."""
from pathlib import Path
import json
import os
import sys
import tempfile
import nbformat
from nbclient import NotebookClient

notebooks = Path(__file__).resolve().parent
root = notebooks.parents[2]
path = notebooks / 'analysis.ipynb'
original = os.environ.get('JUPYTER_PATH')
with tempfile.TemporaryDirectory(prefix='ana068-kernel-') as directory:
    spec = Path(directory) / 'kernels' / 'python3'
    spec.mkdir(parents=True)
    (spec / 'kernel.json').write_text(json.dumps({
        'argv': [sys.executable, '-m', 'ipykernel_launcher', '-f', '{connection_file}'],
        'display_name': 'Python 3 (ANA-068)', 'language': 'python'}))
    os.environ['JUPYTER_PATH'] = directory + (os.pathsep + original if original else '')
    try:
        document = nbformat.read(path, as_version=4)
        client = NotebookClient(document, timeout=900, kernel_name='python3', resources={'metadata': {'path': str(root)}})
        client.execute()
        temporary = path.with_suffix('.executed.ipynb')
        nbformat.write(document, temporary)
        os.replace(temporary, path)
    finally:
        if original is None:
            os.environ.pop('JUPYTER_PATH', None)
        else:
            os.environ['JUPYTER_PATH'] = original
print('Executed ANA-068:', sum(c.cell_type == 'code' for c in document.cells), 'code cells')
