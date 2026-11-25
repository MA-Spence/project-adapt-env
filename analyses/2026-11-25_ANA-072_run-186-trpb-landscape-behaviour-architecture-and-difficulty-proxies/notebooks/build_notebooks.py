"""Build notebooks from reviewable percent-format sources, without executing."""
from pathlib import Path
import hashlib
import nbformat
HERE = Path(__file__).resolve().parent
def source_cells(path):
    cells=[];mode=None;lines=[]
    def flush():
        if mode is None:return
        text='\n'.join(lines).strip('\n')
        if mode=='markdown':
            text='\n'.join(line[2:] if line.startswith('# ') else line[1:] if line=='#' else line for line in text.splitlines())
        cells.append((mode,text))
    for line in path.read_text().splitlines():
        if line.startswith('# %%'):
            flush();lines=[];mode='markdown' if '[markdown]' in line else 'code'
        else:lines.append(line)
    flush()
    return cells
def build(stem):
    cells=[]
    for i,(kind,source) in enumerate(source_cells(HERE/f'{stem}.py')):
        make=nbformat.v4.new_markdown_cell if kind=='markdown' else nbformat.v4.new_code_cell
        cell=make(source);cell.id=hashlib.sha256(f'{stem}:{i}:{source}'.encode()).hexdigest()[:12];cells.append(cell)
    nb=nbformat.v4.new_notebook(cells=cells,metadata={'kernelspec':{'name':'python3','display_name':'Python 3 (ANA-072)','language':'python'},'language_info':{'name':'python','version':'3.11.15'}})
    nbformat.write(nb,HERE/f'{stem}.ipynb');print('Built',stem,len(cells),'cells')
if __name__=='__main__':
    import sys
    for stem in (sys.argv[1:] or ['validation','analysis','architecture','difficulty']):build(stem)
