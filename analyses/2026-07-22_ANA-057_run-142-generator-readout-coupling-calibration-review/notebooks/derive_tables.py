"""Derive ANA-057 tables from the collected RUN-142 (EXP-058) outputs.

Reads the DVC-tracked compact outputs and writes the analysis tables and
metrics.json used by ANA-057 / RES-057. Applies the preregistered EXP-058
decision rule only; adds no new success criteria.
"""
from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data/processed/generator-readout-coupling-calibration/RUN-142/outputs"
ANA = ROOT / "analyses/2026-07-22_ANA-057_run-142-generator-readout-coupling-calibration-review"
TAB = ANA / "tables"


def write(df, name):
    df.to_csv(TAB / f"{name}.csv", index=False)
    (TAB / f"{name}.md").write_text(df.to_markdown(index=False))


def main():
    g = pd.read_csv(OUT / "core_guardrail_checks.csv")
    sep = pd.read_csv(OUT / "gnk_separation_checks.csv")
    mono = pd.read_csv(OUT / "knob_monotonicity_checks.csv")

    core_rows = []
    for fx in ["A", "B", "C"]:
        s = g[g.fix_label == fx]
        core_rows.append(dict(
            fix_label=fx, n=len(s),
            low_plateau_pass=int(s.low_plateau_pass.sum()),
            high_plateau_pass=int(s.high_plateau_pass.sum()),
            beneficial_decline_pass=int(s.beneficial_decline_pass.sum()),
            core_guardrail_pass=int(s.core_guardrail_pass.sum()),
            obs_high_plateau_fitness_median_max=round(float(s.high_plateau_fitness_median.max()), 4),
        ))
    write(pd.DataFrame(core_rows), "run-142_core_guardrail_summary_by_fix")

    sep_rows = []
    for fx in ["A", "B", "C"]:
        s = sep[sep.fix_label == fx]
        sep_rows.append(dict(
            fix_label=fx, n=len(s),
            high_order_sep_pass=int(s.high_order_fraction_separation_pass.sum()),
            sign_epi_sep_pass=int(s.sign_epistasis_fraction_separation_pass.sum()),
            local_maxima_sep_pass=int(s.local_maxima_density_separation_pass.sum()),
            matched_core_guardrail_pass=int(s.matched_core_guardrail_pass.sum()),
            gnk_separation_pass=int(s.gnk_separation_pass.sum()),
            high_order_absdelta_max=round(float(s.high_order_fraction_full_minus_core.abs().max()), 4),
            high_order_absdelta_median=round(float(s.high_order_fraction_full_minus_core.abs().median()), 4),
            sign_epi_absdelta_max=round(float(s.sign_epistasis_fraction_full_minus_core.abs().max()), 4),
        ))
    write(pd.DataFrame(sep_rows), "run-142_gnk_separation_summary_by_fix")

    mm = mono.groupby(["fix_label", "metric", "expected_direction"]).agg(
        n=("direction_pass", "size"),
        direction_pass=("direction_pass", "sum"),
        pass_frac=("direction_pass", "mean"),
        median_rho=("spearman_rho", "median"),
        wrong_direction_frac=("wrong_direction", "mean"),
    ).round(3).reset_index()
    write(mm, "run-142_monotonicity_by_metric")

    piv = g[g.fix_label == "A"].pivot_table(
        index="threshold_steepness", columns="stability_margin",
        values="high_plateau_fitness_median", aggfunc="max").round(4).reset_index()
    write(piv, "run-142_high_plateau_collapse_fixA")


if __name__ == "__main__":
    main()
