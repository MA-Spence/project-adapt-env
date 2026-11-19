"""Execute the canonical notebook and retain cell outputs; use its Python environment."""
from pathlib import Path
import os
import tempfile
import nbformat
from nbclient import NotebookClient

folder = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "adaptenv-ana064-mpl"))
os.environ.setdefault("JUPYTER_RUNTIME_DIR", str(Path(tempfile.gettempdir()) / "adaptenv-ana064-jupyter"))
os.environ.setdefault("IPYTHONDIR", str(Path(tempfile.gettempdir()) / "adaptenv-ana064-ipython"))
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
path = folder / "analysis.ipynb"
notebook = nbformat.read(path, as_version=4)
client = NotebookClient(notebook, timeout=240, kernel_name="python3", resources={"metadata": {"path": str(folder)}})
try:
    client.execute()
finally:
    nbformat.write(notebook, path)
assert not [o for c in notebook.cells for o in c.get("outputs", []) if o.output_type == "error"]
print(f"Executed {sum(c.cell_type == 'code' for c in notebook.cells)} code cells without errors: {path}")
