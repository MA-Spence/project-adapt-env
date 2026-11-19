#!/usr/bin/env python3
"""EXP-065: a small, direct test of conditional aggregate statistical recovery.

No historical acceptance gates, fitted empirical interactions, seed selection,
or inferred missing genotypes. The scientific unit is a measured binary cube;
calibration and confirmation cubes have disjoint measured genotype sets.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
for path in (ROOT / "src", ROOT / "external/Adapt-Env"):
    sys.path.insert(0, str(path))

import numpy as np
import pandas as pd
import yaml
import scipy
from scipy.spatial.distance import cdist

from adaptenv import ThermodynamicConfig, ThermodynamicLandscape
from adaptenv.thermodynamic import sigmoid
from project_adapt_env.dense_descriptive import exact_summary
from project_adapt_env.summary_statistics import METRIC_VERSION, accessible_fraction, mutation_summary

# Deliberately small: the loss also uses their joint geometry, not only marginal bands.
FEATURES = (
    "fitness_sd", "walsh_epistatic_fraction", "walsh_high_order_fraction",
    "beneficial_diminishing_returns_slope", "reference_beneficial_fraction",
    "neutral_edge_fraction", "sign_epistasis_fraction",
    "reciprocal_sign_epistasis_fraction", "local_maxima_density",
    "accessible_path_fraction",
)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for part in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def write_json(path, value):
    Path(path).write_text(json.dumps(clean(value), indent=2, allow_nan=False) + "\n")


def measured_cubes(path, *, dimension, count, seed, missing_fitness="error"):
    """Uniform observed anchors/sites/alternative alleles; reject missing cubes.

    All retained cubes are genotype-disjoint, including across the eventual
    calibration/confirmation split. No activity filtering or low-index allele
    preference is used. Rejection sampling is conditional on measured completeness.
    """
    frame = pd.read_parquet(path) if Path(path).suffix == ".parquet" else pd.read_csv(path)
    if "fitness_is_higher_better" in frame:
        if not frame["fitness_is_higher_better"].astype(str).str.lower().isin(["true", "1"]).all():
            raise ValueError("input must declare higher-is-better fitness")
    columns = sorted([c for c in frame if c.startswith("pos") and c[3:].isdigit()], key=lambda c: int(c[3:]))
    if not columns:
        raise ValueError(f"{path}: explicit variable-site pos columns are required")
    if len(columns) < dimension:
        raise ValueError("input has fewer variable sites than the requested dimension")
    if frame[columns].isna().any().any():
        raise ValueError("input contains missing genotype coordinates")
    if missing_fitness not in ("error", "exclude"):
        raise ValueError("missing_fitness must be error or exclude")
    input_rows = len(frame)
    finite = np.isfinite(frame["fitness"].to_numpy(float))
    excluded = int((~finite).sum())
    if excluded and missing_fitness == "error":
        raise ValueError("input contains nonfinite fitness; declare exclusion explicitly")
    frame = frame.loc[finite].copy()
    if frame.empty:
        raise ValueError("input has no finite measured fitness")
    grouped = frame.groupby(columns, sort=True, dropna=False)["fitness"].mean()
    score = {tuple(map(str, key if isinstance(key, tuple) else (key,))): float(value)
             for key, value in grouped.items()}
    states = list(score)
    alleles = [sorted({s[i] for s in states}) for i in range(len(columns))]
    mutable = [i for i, a in enumerate(alleles) if len(a) > 1]
    if len(mutable) < dimension:
        raise ValueError("not enough polymorphic sites")
    rng = np.random.default_rng(seed)
    used, cubes, audit = set(), [], []
    for attempt in range(max(10000, count * 1000)):
        anchor = states[int(rng.integers(len(states)))]
        sites = sorted(rng.choice(mutable, dimension, replace=False).tolist())
        alternate = [rng.choice([a for a in alleles[i] if a != anchor[i]]) for i in sites]
        nodes = []
        for mask in range(1 << dimension):
            state = list(anchor)
            for bit, site in enumerate(sites):
                if mask & (1 << bit):
                    state[site] = alternate[bit]
            nodes.append(tuple(state))
        if any(n not in score or n in used for n in nodes):
            continue
        used.update(nodes)
        cubes.append(np.array([score[n] for n in nodes]))
        audit.append({"cube": len(cubes) - 1, "sites": sites, "anchor": anchor,
                      "alternatives": alternate, "node_keys": nodes, "attempt": attempt + 1,
                      "input_rows": input_rows, "excluded_nonfinite_rows": excluded,
                      "duplicate_finite_rows_averaged": len(frame) - len(score)})
        if len(cubes) == count:
            return cubes, audit
    raise ValueError(f"only {len(cubes)}/{count} complete genotype-disjoint cubes available in {path}")


def describe(cube, dimension, resolution):
    result = exact_summary(cube, dim=dimension, resolution=resolution)
    result.update(mutation_summary(cube, dim=dimension, resolution=resolution))
    result["accessible_path_fraction"] = accessible_fraction(cube, dim=dimension, resolution=resolution)
    return {k: v for k, v in result.items() if isinstance(v, (int, float, bool, np.number))}


class JointDistance:
    """Energy distance in a frozen, masked standardized descriptor embedding.

    Undefined values remain NaN in outputs. Internally each feature has a value
    coordinate and an availability indicator. Missing value coordinates are zero
    after empirical centering; their indicator prevents treating them as observed
    zeros. This compares missingness as well as the joint descriptor distribution.
    Scales are fitted on empirical calibration cubes only.
    """
    def __init__(self, rows):
        x = self.matrix(rows)
        self.center = np.array([np.mean(c[np.isfinite(c)]) if np.isfinite(c).any() else 0.0 for c in x.T])
        self.scale = np.array([max(float(np.std(c[np.isfinite(c)])), 0.05) if np.isfinite(c).any() else 0.05 for c in x.T])

    @staticmethod
    def matrix(rows):
        return np.array([[r.get(k, math.nan) for k in FEATURES] for r in rows], dtype=float)

    def embed(self, rows):
        x = self.matrix(rows)
        mask = np.isfinite(x)
        z = np.where(mask, (x - self.center) / self.scale, 0.0)
        return np.concatenate([z, mask.astype(float)], axis=1) / math.sqrt(len(FEATURES))

    def __call__(self, reference, generated):
        x, y = self.embed(reference), self.embed(generated)
        return max(0.0, float(2 * cdist(x, y).mean() - cdist(x, x).mean() - cdist(y, y).mean()))


def candidate_grid(config, family):
    mixtures = config["model"]["order_mixtures"]
    if family == "bare_gnk":
        return [{"family": family, "mixture": i, "scale": scale, "offset": offset}
                for i in range(len(mixtures)) for scale in config["fit"]["bare_scales"]
                for offset in config["fit"]["bare_offsets"]]
    return [{"family": family, "mixture": i, "beta": b, "p0": p, "ratio": r}
            for i, b, p, r in itertools.product(
                range(len(mixtures)) if family == "composite" else [0],
                config["fit"]["steepness"], config["fit"]["reference_occupancy"],
                config["fit"]["interaction_ratio"] if family == "composite" else [0.0])]


def latent_panel(config, *, seed, dimension, mixture, count):
    # Effective active-site model, not a claim that a 4-cube identifies the other
    # hundreds of residues. The separate length diagnostic tests scalable queries.
    model = ThermodynamicLandscape(ThermodynamicConfig(
        length=dimension, seed=seed, reference_rank=config["model"]["reference_rank"],
        representation_rank=config["model"]["representation_rank"],
        order_weights=tuple(config["model"]["order_mixtures"][mixture])))
    rng = np.random.default_rng(seed + 87111)
    nodes = []
    for _ in range(count):
        anchor = rng.integers(0, 20, dimension)
        alt = (anchor + rng.integers(1, 20, dimension)) % 20
        for mask in range(1 << dimension):
            nodes.append(np.where([(mask >> i) & 1 for i in range(dimension)], alt, anchor))
    a, g = model.latent_components(np.array(nodes))
    return a.reshape(count, -1), g.reshape(count, -1), model.fingerprint


def predictions(parameter, panels, *, dimension, resolution):
    rows, values = [], []
    for seed, (a, g, _) in panels.items():
        if parameter["family"] == "bare_gnk":
            cubes = parameter["offset"] + parameter["scale"] * g
        else:
            mu = math.log(parameter["p0"] / (1 - parameter["p0"]))
            cubes = sigmoid(mu + parameter["beta"] * (a + parameter["ratio"] * g))
        for index, cube in enumerate(cubes):
            rows.append({"seed": seed, "cube": index, **describe(cube, dimension, resolution)})
            values.append(cube)
    return rows, values


def scaling_diagnostic(config):
    rows = []
    for length in config["scaling"]["lengths"]:
        for seed in config["scaling"]["seeds"]:
            start = time.perf_counter()
            model = ThermodynamicLandscape(ThermodynamicConfig(length=length, seed=seed))
            build_time = time.perf_counter() - start
            rng = np.random.default_rng(seed + 4711)
            sources = rng.integers(0, 20, (config["scaling"]["mutations"], length))
            sites = rng.integers(0, length, len(sources))
            targets = sources.copy()
            targets[np.arange(len(sources)), sites] = (targets[np.arange(len(sources)), sites] + rng.integers(1, 20, len(sources))) % 20
            start = time.perf_counter()
            a0, g0 = model.latent_components(sources); a1, g1 = model.latent_components(targets)
            query_time = time.perf_counter() - start
            f0, f1 = model.evaluate_batch(sources), model.evaluate_batch(targets)
            rows.append({"length": length, "seed": seed, "build_seconds": build_time,
                         "paired_query_seconds": query_time, "n_mutations": len(sources),
                         "additive_mutation_rms": np.sqrt(np.mean((a1-a0)**2)),
                         "interaction_mutation_rms": np.sqrt(np.mean((g1-g0)**2)),
                         "observable_mutation_rms": np.sqrt(np.mean((f1-f0)**2)),
                         "n_supports": model.summary()["n_supports"], "fingerprint": model.fingerprint,
                         "scope": "uniform_q20_backgrounds_not_empirical_confirmation"})
    return rows


def run(config, output, *, quick=False, smoke=False, data_root=None):
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    # Copy before reducing budgets, so callers cannot accidentally mutate their config.
    config = json.loads(json.dumps(config))
    if data_root is not None:
        config["data"]["root"] = str(Path(data_root).resolve())
    if quick or smoke:
        config["sampling"]["calibration_cubes"] = 4
        config["sampling"]["confirmation_cubes"] = 4
        config["sampling"]["synthetic_cubes_per_seed"] = 3
        # Reserved engineering fixtures: never expose production confirmation seeds.
        config["sampling"]["seed"] = 69001
        config["fit"]["calibration_seeds"] = [69101, 69102]
        config["fit"]["confirmation_seeds"] = [69201, 69202]
        for key in ("steepness", "reference_occupancy", "interaction_ratio", "bare_scales", "bare_offsets"):
            config["fit"][key] = config["fit"][key][:2]
        config["model"]["order_mixtures"] = config["model"]["order_mixtures"][:1]
        config["scaling"] = {"lengths": [16, 32], "seeds": [811], "mutations": 64}
    train_seeds, test_seeds = config["fit"]["calibration_seeds"], config["fit"]["confirmation_seeds"]
    if set(train_seeds) & set(test_seeds):
        raise ValueError("generator calibration and confirmation seeds must be disjoint")
    dim = int(config["sampling"]["dimension"])
    ntrain, ntest = config["sampling"]["calibration_cubes"], config["sampling"]["confirmation_cubes"]
    resolution = float(config["measurement"]["resolution"])
    cache, fingerprints = {}, {}
    for mixture in range(len(config["model"]["order_mixtures"])):
        for seed in train_seeds + test_seeds:
            cache[mixture, seed] = latent_panel(config, seed=seed, dimension=dim, mixture=mixture,
                                               count=config["sampling"]["synthetic_cubes_per_seed"])
            fingerprints[f"{mixture}:{seed}"] = cache[mixture, seed][2]
    inputs, empirical_rows, calibration_rows, confirmation_rows, sensitivity_rows, outcomes, selections = {}, [], [], [], [], [], []
    datasets = ["synthetic_smoke"] if smoke else config["data"]["datasets"]
    for dataset in datasets:
        if smoke:
            rng = np.random.default_rng(18201)
            cubes = [sigmoid(rng.normal(size=1 << dim)) for _ in range(ntrain + ntest)]
            audit = [{"cube": i, "node_keys": [[str(i), str(n)] for n in range(1 << dim)]} for i in range(len(cubes))]
        else:
            path = ROOT / config["data"]["root"] / "datasets" / dataset / "variants.csv"
            if not path.is_file():
                path = path.with_suffix(".parquet")
            input_name = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
            inputs[input_name] = digest(path)
            cubes, audit = measured_cubes(path, dimension=dim, count=ntrain + ntest,
                seed=config["sampling"]["seed"], missing_fitness=config["data"]["missing_fitness"])
        raw = np.concatenate(cubes[:ntrain])
        low, high = float(raw.min()), float(raw.max())
        scale = high - low
        # A completely flat calibration panel has no estimated scale: preserve it
        # using a declared unit, label it, and let missingness/plateau loss report it.
        if scale <= 1e-15:
            scale = 1.0
        normalized = [(cube-low)/scale for cube in cubes]  # no clipping or inferred values
        empirical = [describe(c, dim, resolution) for c in normalized]
        train, test = empirical[:ntrain], empirical[ntrain:]
        for i, row in enumerate(empirical):
            empirical_rows.append({"dataset": dataset, "split": "calibration" if i < ntrain else "confirmation",
                                   "cube": i, "raw_offset": low, "raw_scale": scale,
                                   "node_keys": json.dumps(audit[i]["node_keys"]),
                                   "raw_values": json.dumps(cubes[i].tolist()),
                                   "sampling_audit": json.dumps(audit[i]), **row})
        distance = JointDistance(train)
        family_losses = {}
        for family in ("additive", "bare_gnk", "composite"):
            candidates = candidate_grid(config, family)
            scored = []
            for parameter in candidates:
                panels = {s: cache[parameter["mixture"], s] for s in train_seeds}
                rows, _ = predictions(parameter, panels, dimension=dim, resolution=resolution)
                loss = distance(train, rows)
                scored.append((loss, parameter))
                calibration_rows.append({"dataset": dataset, "family": family, "loss": loss,
                                         "parameters": json.dumps(parameter, sort_keys=True)})
            loss, selected = min(scored, key=lambda item: item[0])
            selections.append({"dataset": dataset, "family": family, "parameters": selected, "calibration_loss": loss})
            panels = {s: cache[selected["mixture"], s] for s in test_seeds}
            generated, generated_cubes = predictions(selected, panels, dimension=dim, resolution=resolution)
            heldout = distance(test, generated)
            family_losses[family] = heldout
            for row in generated:
                confirmation_rows.append({"dataset": dataset, "family": family, **row})
            seed_losses = [distance(test, [r for r in generated if r["seed"] == seed]) for seed in test_seeds]
            empirical_mean = {k: pd.Series([r[k] for r in test]).mean() for k in FEATURES}
            generated_mean = {k: pd.Series([r[k] for r in generated]).mean() for k in FEATURES}
            outcomes.append({"dataset": dataset, "family": family, "calibration_loss": loss,
                             "confirmation_joint_loss": heldout, "confirmation_seed_loss_mean": np.mean(seed_losses),
                             "confirmation_seed_loss_sd": np.std(seed_losses), "parameters": selected,
                             "empirical_confirmation_mean": empirical_mean,
                             "generated_confirmation_mean": generated_mean,
                             "descriptor_mean_error": {k: generated_mean[k] - empirical_mean[k] for k in FEATURES},
                             "empirical_resolved_fraction": {k: np.mean([np.isfinite(r[k]) for r in test]) for k in FEATURES},
                             "generated_resolved_fraction": {k: np.mean([np.isfinite(r[k]) for r in generated]) for k in FEATURES}})
            for tol in config["measurement"]["sensitivity_resolutions"]:
                empirical_s = [describe(c, dim, tol) for c in normalized[ntrain:]]
                generated_s = [describe(c, dim, tol) for c in generated_cubes]
                sensitivity_rows.append({"dataset": dataset, "family": family, "resolution": tol,
                                         "confirmation_joint_loss": distance(empirical_s, generated_s)})
        print(f"{dataset}: held-out joint distances {family_losses}", flush=True)
    pd.DataFrame(empirical_rows).to_csv(output / "empirical_cubes.csv", index=False)
    pd.DataFrame(calibration_rows).to_csv(output / "calibration.csv", index=False)
    pd.DataFrame(confirmation_rows).to_csv(output / "confirmation.csv", index=False)
    pd.DataFrame(sensitivity_rows).to_csv(output / "resolution_sensitivity.csv", index=False)
    pd.DataFrame(scaling_diagnostic(config)).to_csv(output / "scaling.csv", index=False)
    write_json(output / "selected_parameters.json", selections)
    sources = [Path(__file__), ROOT / "src/project_adapt_env/dense_descriptive.py",
               ROOT / "src/project_adapt_env/summary_statistics.py", *[(ROOT / "external/Adapt-Env/adaptenv" / name)
               for name in ("thermodynamic.py", "extensive_gnk.py", "ranked_categorical_gnk.py")]]
    write_json(output / "manifest.json", {"metric_version": METRIC_VERSION, "config": config,
               "input_sha256": inputs, "source_sha256": {str(p.relative_to(ROOT)): digest(p) for p in sources},
               "generator_fingerprints": fingerprints, "numpy": np.__version__, "pandas": pd.__version__,
               "scipy": scipy.__version__, "pyyaml": yaml.__version__,
               "python": sys.version, "quick": quick, "synthetic_smoke": smoke})
    result = {"experiment": "EXP-065", "status": "SMOKE_ONLY" if quick or smoke else "DESCRIPTIVE_CONFIRMATION_COMPLETE",
              "outcomes": outcomes, "features": FEATURES,
              "scope": "within-assay conditional aggregate recovery on measured complete subspaces",
              "interpretation": "Report joint and descriptor-specific held-out discrepancies; no inherited pass gate.",
              "diminishing_returns": "descriptive: inputs lack independent replicate-error correction",
              "generalization": "fresh seeds and genotype-disjoint cubes; not held-out proteins or full protein WHT"}
    write_json(output / "summary.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--data-root", type=Path, help="existing empirical bundle root; recorded in the manifest")
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--smoke", action="store_true", help="synthetic fixtures only; never scientific evidence")
    args = parser.parse_args()
    run(yaml.safe_load(args.config.read_text()), args.output_dir, quick=args.quick, smoke=args.smoke, data_root=args.data_root)
