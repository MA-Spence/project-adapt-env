"""Execute ANA-074 notebooks in fresh kernels using this exact interpreter."""
from pathlib import Path
import json,os,sys,tempfile
import nbformat
from nbclient import NotebookClient
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
original=os.environ.get('JUPYTER_PATH')
with tempfile.TemporaryDirectory(prefix='ana074-kernel-') as directory:
    spec=Path(directory)/'kernels/python3';spec.mkdir(parents=True)
    (spec/'kernel.json').write_text(json.dumps({'argv':[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'],'display_name':'Python 3 (ANA-074)','language':'python'}))
    os.environ['JUPYTER_PATH']=directory+(os.pathsep+original if original else '')
    try:
        for stem in (sys.argv[1:] or ['analysis']):
            path=HERE/f'{stem}.ipynb';nb=nbformat.read(path,as_version=4)
            NotebookClient(nb,timeout=900,kernel_name='python3',resources={'metadata':{'path':str(ROOT)}}).execute()
            temporary=path.with_suffix('.executed.ipynb');nbformat.write(nb,temporary);os.replace(temporary,path)
            print('Executed',stem,sum(c.cell_type=='code' for c in nb.cells),'code cells',flush=True)
    finally:
        if original is None:os.environ.pop('JUPYTER_PATH',None)
        else:os.environ['JUPYTER_PATH']=original
