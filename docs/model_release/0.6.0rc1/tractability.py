"""Bounded local engineering smoke check, not an empirical/comparative experiment.

Run from any directory with Python containing the package's dependencies:
    python docs/model_release/0.6.0rc1/tractability.py --output /tmp/tractability.json

Each length runs in a fresh process. Timings are single measurements on this
host, not HYP-010 evidence. The synthetic frequency table exercises the profile
interface without implying a real biological alignment.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import resource
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "external/Adapt-Env"))


def worker(length):
    from adaptenv import ThermodynamicConfig, ThermodynamicLandscape, __version__
    from adaptenv.amino_acids import AA_ORDER

    rng = np.random.default_rng(20261119 + length)
    frequencies = rng.dirichlet(np.ones(20), size=length)
    sequences = rng.integers(20, size=(256, length))
    prior = dict(frequencies=frequencies, alphabet=AA_ORDER,
                 uniform_mixture=.01, source_precision="synthetic smoke fixture")
    start = time.perf_counter()
    model = ThermodynamicLandscape(ThermodynamicConfig(length=length, seed=19),
                                   profile_prior=prior)
    build_seconds = time.perf_counter() - start
    start = time.perf_counter()
    components = model.evaluate_batch_components(sequences)
    query_seconds = time.perf_counter() - start
    assert all(np.isfinite(values).all() for values in components.values())
    sites = np.arange(length)
    alleles = (sequences[0] + 1) % 20
    start = time.perf_counter()
    deltas = model.mutation_effects(sequences[0], sites, alleles)
    mutation_seconds = time.perf_counter() - start
    mutants = np.repeat(sequences[0][None], 16, axis=0)
    mutants[np.arange(16), sites[:16]] = alleles[:16]
    expected = model.evaluate_batch(mutants) - model.evaluate(sequences[0])
    error = float(np.max(np.abs(deltas[:16] - expected)))
    assert error < 1e-12
    tables_before = model.field.tables.nbytes
    model.evaluate_batch(sequences[::-1])
    assert model.field.tables.nbytes == tables_before
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return dict(
        version=__version__, length=length, queries=len(sequences), seed=19,
        target_degree=6, supports=model.summary()["n_supports"],
        table_bytes=tables_before, additive_bytes=model.additive_matrix.nbytes,
        profile_conditioned=model.summary()["profile_conditioned"],
        mean_adjacent_hamming_fraction=float(np.mean(sequences[:-1] != sequences[1:])),
        build_seconds=build_seconds, query_seconds=query_seconds,
        mutation_queries=len(sites), mutation_seconds_including_baseline=mutation_seconds,
        direct_mutation_check_max_absolute_error=error,
        peak_process_rss_bytes=int(peak if sys.platform == "darwin" else peak * 1024),
        fingerprint=model.fingerprint,
    )


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--length", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.length:
        print(json.dumps(worker(args.length)))
        return
    results = []
    for length in [100, 500, 1000]:
        result = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--length", str(length)],
            check=True, capture_output=True, text=True, timeout=120,
        )
        results.append(json.loads(result.stdout))
    report = dict(
        scope="Engineering smoke check only; no comparative performance or biological claim",
        platform=platform.platform(), python=platform.python_version(), numpy=np.__version__,
        query_distribution="independent uniform q20, not sampled close to a reference",
        profile_source="synthetic Dirichlet frequency fixture; biological accuracy not tested",
        peak_rss_scope="whole fresh process, including imports, inputs and validation scratch",
        results=results,
    )
    content = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(content)
    else:
        print(content, end="")


if __name__ == "__main__":
    main()
