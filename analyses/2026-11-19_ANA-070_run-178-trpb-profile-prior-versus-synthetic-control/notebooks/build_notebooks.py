"""Build both notebooks from their reviewable percent-format Python sources."""
from pathlib import Path
import nbformat
HERE = Path(__file__).resolve().parent
for stem in ['analysis', 'validation']:
    cells=[]; mode=None; lines=[]
    def flush():
        if mode is None: return
        text='\n'.join(lines).strip('\n')
        if mode=='markdown':
            text='\n'.join(line[2:] if line.startswith('# ') else line[1:] if line=='#' else line for line in text.splitlines())
            cells.append(nbformat.v4.new_markdown_cell(text))
        else: cells.append(nbformat.v4.new_code_cell(text))
    for line in (HERE/f'{stem}.py').read_text().splitlines():
        if line.startswith('# %%'):
            flush(); lines=[]; mode='markdown' if '[markdown]' in line else 'code'
        else: lines.append(line)
    flush()
    nb=nbformat.v4.new_notebook(cells=cells, metadata={'kernelspec':{'name':'python3','display_name':'Python 3 (ANA-070)','language':'python'},'language_info':{'name':'python','version':'3.11.15'}})
    nbformat.write(nb,HERE/f'{stem}.ipynb')
    print('Built',stem,len(cells),'cells')
