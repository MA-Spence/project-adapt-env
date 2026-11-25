"""Compare old maps against frozen RUN-180 using the same Python interpreter."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
parser = argparse.ArgumentParser(__doc__)
parser.add_argument("--probe", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

if args.probe is not None:
    sys.path.insert(0, str(args.probe.resolve()))
    import numpy as np
    from adaptenv import ThermodynamicConfig, ThermodynamicLandscape, build_landscape_from_preset
    names = ["trpb_r0_sigmoid_independent_difference", "trpb_r0_sigmoid_shared_stop_sensitivity",
             "trpb_r4_sigmoid_independent_difference", "trpb_r4_sigmoid_shared_stop_sensitivity",
             "trpb_run175_r05_shared_stop_intrinsic", "trpb_run178_r05_profile_shared_stop_intrinsic"]
    results = {}
    for name in names + ["full_pairwise"]:
        for seed in (0, 1, 7, 42):
            model = (ThermodynamicLandscape(ThermodynamicConfig(length=4, seed=seed, additive_prior="synthetic"))
                     if name == "full_pairwise" else build_landscape_from_preset(name, seed=seed))
            sequences = np.random.default_rng(41).integers(20, size=(200, 4))
            arrays = [model.additive_matrix, model.reference, *model.latent_components(sequences),
                      model.evaluate_batch(sequences), model.mutation_effects(
                          sequences[0], sites=[0, 1, 2, 3], alleles=[1, 2, 3, 4])]
            digest = hashlib.sha256()
            for value in arrays:
                digest.update(np.asarray(value).tobytes())
            results[f"{name}/{seed}"] = dict(fingerprint=model.fingerprint,
                                            components_and_fitness_sha256=digest.hexdigest())
    print(json.dumps(results, sort_keys=True))
else:
    if args.output is None:
        parser.error("--output is required")
    current = ROOT / "external/Adapt-Env"
    reference = ROOT / "data/processed/aggregate-thermodynamic-gnk-calibration/RUN-180/source/external/Adapt-Env"
    def probe(path):
        return json.loads(subprocess.check_output([sys.executable, str(Path(__file__).resolve()),
                                                   "--probe", str(path)], text=True))
    expected, actual = probe(reference), probe(current)
    if expected != actual:
        raise SystemExit("Historical maps differ from frozen RUN-180 source")
    import numpy as np
    report = dict(status="passed", cases=len(actual), python=sys.version,
                  numpy=np.__version__, reference_source=str(reference.relative_to(ROOT)),
                  comparisons=actual)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"All {len(actual)} maps and fingerprints match the frozen source")
