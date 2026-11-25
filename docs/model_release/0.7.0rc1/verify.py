"""Reproducible source/wheel engineering verification (not a biological RUN).

Use Python -S to exclude unrelated editable-install import hooks; explicitly
supply the isolated environment and base scientific dependency directories.
"""
from pathlib import Path
import argparse
import sys

parser = argparse.ArgumentParser(__doc__)
parser.add_argument("--package-path", type=Path, required=True)
parser.add_argument("--site", type=Path, action="append", required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--project-tests", action="store_true")
args = parser.parse_args()
ROOT = Path(__file__).resolve().parents[3]
package_path = args.package_path.resolve()
sys.path[:0] = [str(package_path)] + [str(p.resolve()) for p in args.site]

import contextlib
import hashlib
import importlib.metadata
import io
import json
import os
import re
import time

import numpy as np
import pytest
import adaptenv
from adaptenv import ThermodynamicConfig, ThermodynamicLandscape

os.environ["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
args.output.mkdir(parents=True, exist_ok=True)
assert Path(adaptenv.__file__).is_relative_to(package_path)
assert adaptenv.__version__ == "0.7.0rc1"
paths = [str(ROOT / "external/Adapt-Env/tests")]
if args.project_tests:
    paths.append(str(ROOT / "tests/test_calibrate_architecture.py"))
start = time.perf_counter()
exit_code = int(pytest.main(paths + ["-q", "--junitxml", str(args.output / "pytest.xml")]))
report = dict(python=sys.version, executable=sys.executable,
              imported_package=adaptenv.__file__, tests_exit_code=exit_code,
              test_seconds=time.perf_counter() - start,
              dependencies={name: importlib.metadata.version(name) for name in
                            ("numpy", "scipy", "optuna", "pytest", "fastapi", "httpx")})
report["all_package_imports_from_selected_location"] = all(
    Path(module.__file__).is_relative_to(package_path)
    for name, module in sys.modules.items()
    if name.startswith("adaptenv") and getattr(module, "__file__", None))
assert report["all_package_imports_from_selected_location"]
if exit_code == 0:
    examples = []
    for relative in ("README.md", "docs/API_AND_SERVICE.md", "docs/MIGRATION_0_6.md", "docs/CALIBRATION.md"):
        scope = {}
        for number, block in enumerate(re.findall(r"```python\n(.*?)```",
                                (ROOT / "external/Adapt-Env" / relative).read_text(), re.S), 1):
            with contextlib.redirect_stdout(io.StringIO()):
                exec(compile(block, f"{relative}:{number}", "exec"), scope)
            examples.append(f"{relative}:{number}")
    report["documentation_examples"] = examples
    report["oracle_checks"] = []
    for length in (100, 500, 1000):
        rng = np.random.default_rng(712)
        matrix = rng.normal(size=(length, 20))
        x = rng.integers(0, 20, size=(256, length))
        start = time.perf_counter()
        model = ThermodynamicLandscape(ThermodynamicConfig(
            length=length, seed=77, interaction_concentration=.2,
            global_intercept=0., interaction_ratio=1.), additive_matrix=matrix)
        construction = time.perf_counter() - start
        start = time.perf_counter()
        y = model.evaluate_batch(x)
        scoring = time.perf_counter() - start
        assert np.isfinite(y).all()
        np.testing.assert_allclose(model.field.expected_mutation_energy(), model.sigma_a ** 2, rtol=1e-12)
        np.testing.assert_array_equal(y[:8], model.evaluate_batch(x[:8]))
        report["oracle_checks"].append(dict(length=length, queries=256,
            pair_tables=len(model.field.tables), table_bytes=model.field.tables.nbytes,
            construction_seconds=construction, scoring_seconds=scoring,
            exact_energy_check=True, batch_replay_check=True))
    report["oracle_timing_scope"] = "single engineering checks; not a comparative performance benchmark"
(args.output / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
raise SystemExit(exit_code)
