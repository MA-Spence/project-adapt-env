#!/usr/bin/env python3
"""Download and normalize datasets for ensemble landscape validation.

This script intentionally does not touch labproj records or registries. It
writes external source mirrors under data/external and normalized analysis-ready
tables under data/processed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import textwrap
import zipfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import requests
import yaml


EXTERNAL_ROOT = Path("data/external/ensemble_landscape_validation_sources")
PROCESSED_ROOT = Path("data/processed/ensemble-landscape-validation-datasets")

GRAPHFLA_COMMIT = "main"
PROTEINGYM_VERSION = "v1.3"
HF_TRPB_COMMIT = "9603f13a275fb518182d3c3743049c6a3955ee3d"


@dataclass(frozen=True)
class SourceSpec:
    source_id: str
    url: str
    rel_path: str
    publication: str
    doi: str | None = None
    source_kind: str = "dataset"
    mirror_note: str | None = None
    expected_size: int | None = None
    tags: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class GraphFLADataset:
    dataset_id: str
    source_id: str
    rel_path: str
    landscape_name: str
    tier: str
    readout: str
    publication: str
    doi: str
    alphabet: str
    notes: str


def slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(spec: SourceSpec, root: Path, force: bool = False) -> dict[str, Any]:
    dest = root / spec.rel_path
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force:
        return file_entry(spec, dest, status="existing")

    tmp = dest.with_suffix(dest.suffix + ".tmp")
    headers = {"User-Agent": "project-adapt-env dataset preparation"}
    with requests.get(spec.url, headers=headers, stream=True, timeout=120) as resp:
        resp.raise_for_status()
        with tmp.open("wb") as handle:
            for chunk in resp.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    tmp.replace(dest)
    return file_entry(spec, dest, status="downloaded")


def file_entry(spec: SourceSpec, path: Path, status: str) -> dict[str, Any]:
    stat = path.stat()
    return {
        "source_id": spec.source_id,
        "status": status,
        "path": str(path),
        "size_bytes": stat.st_size,
        "sha256": sha256_file(path),
        "url": spec.url,
        "publication": spec.publication,
        "doi": spec.doi,
        "source_kind": spec.source_kind,
        "mirror_note": spec.mirror_note,
        "tags": spec.tags,
    }


def graphfla_url(path: str) -> str:
    return f"https://raw.githubusercontent.com/COLA-Laboratory/GraphFLA/{GRAPHFLA_COMMIT}/{path}"


def flip_url(path: str) -> str:
    return f"https://raw.githubusercontent.com/J-SNACKKB/FLIP/main/{path}"


GRAPHFLA_SOURCES = [
    ("graphfla_wu2016_gb1", "data/BioSequence/Wu2016_GB1.csv", "Wu et al. 2016 GB1 combinatorial landscape", "10.7554/eLife.16965"),
    ("graphfla_johnston2024_trpb4", "data/BioSequence/Johnston2024_TrpB4.csv", "Johnston et al. 2024 TrpB active-site landscape", "10.1073/pnas.2400439121"),
    ("graphfla_papkou2023_dhfr", "data/BioSequence/Papkou2023_DHFR.csv", "Papkou et al. 2023 DHFR/folA landscape", "10.1126/science.adh3860"),
    ("graphfla_phillips2021_cr9114_h1", "data/BioSequence/Phillips2021_CR9114_h1.csv", "Phillips et al. 2021 CR9114 H1 binding landscape", "10.7554/eLife.71393"),
    ("graphfla_phillips2021_cr9114_h3", "data/BioSequence/Phillips2021_CR9114_h3.csv", "Phillips et al. 2021 CR9114 H3 binding landscape", "10.7554/eLife.71393"),
    ("graphfla_phillips2021_cr9114_flub", "data/BioSequence/Phillips2021_CR9114_fluB.csv", "Phillips et al. 2021 CR9114 FluB binding landscape", "10.7554/eLife.71393"),
    ("graphfla_phillips2021_cr6261_h1", "data/BioSequence/Phillips2021_CR6261_h1.csv", "Phillips et al. 2021 CR6261 H1 binding landscape", "10.7554/eLife.71393"),
    ("graphfla_phillips2021_cr6261_h9", "data/BioSequence/Phillips2021_CR6261_h9.csv", "Phillips et al. 2021 CR6261 H9 binding landscape", "10.7554/eLife.71393"),
    ("graphfla_lite2020_pard2", "data/BioSequence/Lite2020_ParD2.csv", "Lite et al. 2020 ParD2-ParE specificity landscape", "10.7554/eLife.60924"),
    ("graphfla_lite2020_pard3", "data/BioSequence/Lite2020_ParD3.csv", "Lite et al. 2020 ParD3-ParE specificity landscape", "10.7554/eLife.60924"),
]


GRAPHFLA_DATASETS = [
    GraphFLADataset("gb1_wu2016_four_site_graphfla", "graphfla_wu2016_gb1", "graphfla/Wu2016_GB1.csv", "GB1 Wu 2016 four-site", "tier1_complete", "IgG Fc binding fitness", "Wu et al. 2016", "10.7554/eLife.16965", "20-letter protein segment", "Four mutable sites; GraphFLA compact mirror."),
    GraphFLADataset("trpb_johnston2024_four_site_graphfla", "graphfla_johnston2024_trpb4", "graphfla/Johnston2024_TrpB4.csv", "TrpB Johnston 2024 four-site", "tier1_complete", "enzyme activity", "Johnston et al. 2024", "10.1073/pnas.2400439121", "20-letter protein segment", "Four active-site residues."),
    GraphFLADataset("dhfr_papkou2023_nine_nt_graphfla", "graphfla_papkou2023_dhfr", "graphfla/Papkou2023_DHFR.csv", "DHFR/folA Papkou 2023", "tier1_near_complete", "trimethoprim resistance fitness", "Papkou et al. 2023", "10.1126/science.adh3860", "DNA 4-letter segment", "Nine nucleotide genotype; compact mirror avoids 2.9 GB Zenodo archive."),
    GraphFLADataset("cr9114_h1_phillips2021_graphfla", "graphfla_phillips2021_cr9114_h1", "graphfla/Phillips2021_CR9114_h1.csv", "CR9114 H1", "tier1_binary", "antibody binding", "Phillips et al. 2021", "10.7554/eLife.71393", "binary", "Sixteen heavy-chain mutations."),
    GraphFLADataset("cr9114_h3_phillips2021_graphfla", "graphfla_phillips2021_cr9114_h3", "graphfla/Phillips2021_CR9114_h3.csv", "CR9114 H3", "tier1_binary", "antibody binding", "Phillips et al. 2021", "10.7554/eLife.71393", "binary", "Sixteen heavy-chain mutations."),
    GraphFLADataset("cr9114_flub_phillips2021_graphfla", "graphfla_phillips2021_cr9114_flub", "graphfla/Phillips2021_CR9114_fluB.csv", "CR9114 FluB", "tier1_binary", "antibody binding", "Phillips et al. 2021", "10.7554/eLife.71393", "binary", "Sixteen heavy-chain mutations."),
    GraphFLADataset("cr6261_h1_phillips2021_graphfla", "graphfla_phillips2021_cr6261_h1", "graphfla/Phillips2021_CR6261_h1.csv", "CR6261 H1", "tier1_binary", "antibody binding", "Phillips et al. 2021", "10.7554/eLife.71393", "binary", "Eleven heavy-chain mutations."),
    GraphFLADataset("cr6261_h9_phillips2021_graphfla", "graphfla_phillips2021_cr6261_h9", "graphfla/Phillips2021_CR6261_h9.csv", "CR6261 H9", "tier1_binary", "antibody binding", "Phillips et al. 2021", "10.7554/eLife.71393", "binary", "Eleven heavy-chain mutations."),
    GraphFLADataset("pard2_lite2020_graphfla", "graphfla_lite2020_pard2", "graphfla/Lite2020_ParD2.csv", "ParD2-ParE specificity", "tier1_combinatorial", "toxin-antitoxin specificity", "Lite et al. 2020", "10.7554/eLife.60924", "20-letter protein segment", "Three mutable positions in GraphFLA representation."),
    GraphFLADataset("pard3_lite2020_graphfla", "graphfla_lite2020_pard3", "graphfla/Lite2020_ParD3.csv", "ParD3-ParE specificity", "tier1_combinatorial", "toxin-antitoxin specificity", "Lite et al. 2020", "10.7554/eLife.60924", "20-letter protein segment", "Three mutable positions in GraphFLA representation."),
]


PROTEINGYM_ASSAY_IDS = [
    "SPG1_STRSG_Olson_2014",
    "SPG1_STRSG_Wu_2016",
    "HIS7_YEAST_Pokusaeva_2019",
    "GFP_AEQVI_Sarkisyan_2016",
    "F7YBW8_MESOW_Aakre_2015",
    "F7YBW8_MESOW_Ding_2023",
    "PHOT_CHLRE_Chen_2023",
]


SOURCES = [
    *[
        SourceSpec(
            source_id=source_id,
            url=graphfla_url(path),
            rel_path=f"graphfla/{Path(path).name}",
            publication=publication,
            doi=doi,
            source_kind="compact_mirror_csv",
            mirror_note="GraphFLA normalized mirror of the cited empirical landscape.",
            tags=["GraphFLA", "tier1"],
        )
        for source_id, path, publication, doi in GRAPHFLA_SOURCES
    ],
    SourceSpec(
        source_id="proteingym_v1_3_reference",
        url="https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/main/reference_files/DMS_substitutions.csv",
        rel_path="proteingym/v1.3/reference/DMS_substitutions.csv",
        publication="Notin et al. 2023 ProteinGym substitution benchmark",
        doi="10.48550/arXiv.2307.01183",
        source_kind="benchmark_metadata_csv",
        tags=["ProteinGym", "tier3"],
    ),
    SourceSpec(
        source_id="proteingym_v1_3_substitutions",
        url="https://proteingym.s3.us-east-2.amazonaws.com/DMS_substitutions.parquet",
        rel_path="proteingym/v1.3/DMS_substitutions.parquet",
        publication="Notin et al. 2023 ProteinGym substitution benchmark",
        doi="10.48550/arXiv.2307.01183",
        source_kind="benchmark_parquet",
        tags=["ProteinGym", "tier2", "tier3"],
    ),
    SourceSpec(
        source_id="wchnicholas_protein_g_four_mutants_repo",
        url="https://github.com/wchnicholas/ProteinGFourMutants/archive/refs/heads/master.zip",
        rel_path="gb1_wu2016/ProteinGFourMutants-master.zip",
        publication="Wu et al. 2016 GB1 four-mutant source repository",
        doi="10.7554/eLife.16965",
        source_kind="source_repository_zip",
        tags=["GB1", "Wu2016", "Olson2014"],
    ),
    SourceSpec(
        source_id="hf_saprothub_trpb_dataset",
        url=f"https://huggingface.co/datasets/SaProtHub/Dataset-TrpB_fitness_landsacpe/resolve/{HF_TRPB_COMMIT}/dataset.csv",
        rel_path="trpb_johnston2024/hf_saprothub_dataset.csv",
        publication="SaProtHub mirror of Johnston et al. 2024 TrpB dataset",
        doi="10.1073/pnas.2400439121",
        source_kind="full_sequence_mirror_csv",
        mirror_note="Hugging Face mirror with full Tm9D8* sequences and train/valid/test split labels.",
        tags=["TrpB", "tier1", "split"],
    ),
    SourceSpec(
        source_id="dhfr_papkou2023_zenodo_archive_metadata",
        url="https://zenodo.org/api/records/8228920",
        rel_path="dhfr_papkou2023/zenodo_8228920_metadata.json",
        publication="Papkou et al. 2023 analysis/code Zenodo record",
        doi="10.5281/zenodo.8228920",
        source_kind="zenodo_metadata_json",
        mirror_note="The full Zenodo archive is 2.9 GB; compact CSV is retrieved via GraphFLA for immediate processing.",
        tags=["DHFR", "Papkou2023"],
    ),
    SourceSpec(
        source_id="wagner_fitness_landscape_sampling_papkou2023_tsv",
        url="https://raw.githubusercontent.com/andreas-wagner-uzh/fitness_landscape_sampling/main/fitness_data_science_papkou2023.tsv",
        rel_path="dhfr_papkou2023/fitness_data_science_papkou2023.tsv",
        publication="Papkou et al. 2023 DHFR/folA landscape, Wagner-lab sampling repository",
        doi="10.1126/science.adh3860",
        source_kind="full_landscape_tsv",
        mirror_note="Compact public TSV used by Wagner-lab follow-up analyses; contains DNA genotype, amino-acid genotype, fitness, and Wald SE.",
        tags=["DHFR", "Papkou2023", "tier1"],
    ),
    *[
        SourceSpec(
            source_id=f"flip_{slug(path)}",
            url=flip_url(path),
            rel_path=f"flip/{path.removeprefix('splits/')}",
            publication="Dallago et al. 2021 FLIP benchmark",
            doi="10.48550/arXiv.2111.06772",
            source_kind="benchmark_split_or_data",
            tags=["FLIP", "tier3"],
        )
        for path in [
            "splits/gb1/README.md",
            "splits/gb1/5LDE_1.fasta",
            "splits/gb1/four_mutations_full_data.csv.zip",
            "splits/gb1/splits.zip",
            "splits/aav/README.md",
            "splits/aav/P03135.fasta",
            "splits/aav/full_data.csv.zip",
            "splits/aav/splits.zip",
            "splits/meltome/README.md",
            "splits/meltome/full_dataset.json.zip",
            "splits/meltome/full_dataset_sequences.fasta.zip",
            "splits/meltome/meltome_PIDE20_clusters.tsv",
            "splits/meltome/splits.zip",
        ]
    ],
]


def normalize_graphfla(spec: GraphFLADataset, source_root: Path, processed_root: Path) -> dict[str, Any]:
    source_path = source_root / spec.rel_path
    df = pd.read_csv(source_path)
    pos_cols = [c for c in df.columns if re.fullmatch(r"pos\d+", str(c))]
    if "sequences" not in df.columns or "fitness" not in df.columns:
        raise ValueError(f"Unexpected GraphFLA schema for {source_path}: {df.columns.tolist()}")

    out = pd.DataFrame(
        {
            "dataset_id": spec.dataset_id,
            "source_id": spec.source_id,
            "tier": spec.tier,
            "landscape_name": spec.landscape_name,
            "readout": spec.readout,
            "sequence": df["sequences"].astype(str),
            "genotype": df["sequences"].astype(str),
            "fitness": pd.to_numeric(df["fitness"], errors="coerce"),
            "fitness_is_higher_better": True,
            "alphabet": spec.alphabet,
            "n_variable_sites": len(pos_cols),
            "source_row": range(len(df)),
        }
    )
    for c in pos_cols:
        out[c] = df[c].astype(str)

    return write_dataset(
        processed_root,
        spec.dataset_id,
        out,
        {
            "dataset_id": spec.dataset_id,
            "landscape_name": spec.landscape_name,
            "source_ids": [spec.source_id],
            "publication": spec.publication,
            "doi": spec.doi,
            "tier": spec.tier,
            "readout": spec.readout,
            "notes": spec.notes,
        },
    )


def write_dataset(root: Path, dataset_id: str, df: pd.DataFrame, metadata: dict[str, Any]) -> dict[str, Any]:
    dest = root / "datasets" / dataset_id
    dest.mkdir(parents=True, exist_ok=True)
    csv_path = dest / "variants.csv"
    parquet_path = dest / "variants.parquet"
    metadata_path = dest / "metadata.yaml"
    df.to_csv(csv_path, index=False)
    df.to_parquet(parquet_path, index=False)
    with metadata_path.open("w") as handle:
        yaml.safe_dump({**metadata, "n_rows": int(len(df)), "columns": list(df.columns)}, handle, sort_keys=False)
    return {
        "dataset_id": dataset_id,
        "n_rows": int(len(df)),
        "csv_path": str(csv_path),
        "parquet_path": str(parquet_path),
        "metadata_path": str(metadata_path),
        "sha256_csv": sha256_file(csv_path),
        "sha256_parquet": sha256_file(parquet_path),
        **metadata,
    }


def normalize_hf_trpb(source_root: Path, processed_root: Path) -> dict[str, Any]:
    source_path = source_root / "trpb_johnston2024/hf_saprothub_dataset.csv"
    df = pd.read_csv(source_path)
    required = {"protein", "label", "stage"}
    if not required.issubset(df.columns):
        raise ValueError(f"Unexpected TrpB HF schema: {df.columns.tolist()}")
    out = pd.DataFrame(
        {
            "dataset_id": "trpb_johnston2024_full_sequence_hf",
            "source_id": "hf_saprothub_trpb_dataset",
            "tier": "tier1_complete_split",
            "landscape_name": "TrpB Johnston 2024 full-sequence mirror",
            "readout": "enzyme activity",
            "sequence": df["protein"].astype(str),
            "genotype": pd.NA,
            "fitness": pd.to_numeric(df["label"], errors="coerce"),
            "fitness_is_higher_better": True,
            "split": df["stage"].astype(str),
            "alphabet": "20-letter full protein sequence",
            "n_variable_sites": 4,
            "source_row": range(len(df)),
        }
    )
    return write_dataset(
        processed_root,
        "trpb_johnston2024_full_sequence_hf",
        out,
        {
            "source_ids": ["hf_saprothub_trpb_dataset"],
            "publication": "Johnston et al. 2024; SaProtHub mirror",
            "doi": "10.1073/pnas.2400439121",
            "tier": "tier1_complete_split",
            "readout": "enzyme activity",
            "notes": "Full-sequence mirror with split labels; use GraphFLA TrpB4 for compact exact 4-site descriptor work.",
        },
    )


def normalize_papkou_full_tsv(source_root: Path, processed_root: Path) -> dict[str, Any]:
    source_path = source_root / "dhfr_papkou2023/fitness_data_science_papkou2023.tsv"
    df = pd.read_csv(source_path, sep="\t")
    required = {"dna", "aa", "fitness", "wald_se"}
    if not required.issubset(df.columns):
        raise ValueError(f"Unexpected Papkou TSV schema: {df.columns.tolist()}")
    out = pd.DataFrame(
        {
            "dataset_id": "dhfr_papkou2023_full_wagner_tsv",
            "source_id": "wagner_fitness_landscape_sampling_papkou2023_tsv",
            "tier": "tier1_near_complete",
            "landscape_name": "DHFR/folA Papkou 2023 full TSV",
            "readout": "trimethoprim resistance fitness",
            "sequence": df["aa"].astype(str),
            "genotype": df["dna"].astype(str),
            "fitness": pd.to_numeric(df["fitness"], errors="coerce"),
            "fitness_wald_se": pd.to_numeric(df["wald_se"], errors="coerce"),
            "fitness_is_higher_better": True,
            "alphabet": "DNA 4-letter segment with translated amino-acid state",
            "n_variable_sites": 9,
            "source_row": range(len(df)),
        }
    )
    return write_dataset(
        processed_root,
        "dhfr_papkou2023_full_wagner_tsv",
        out,
        {
            "source_ids": ["wagner_fitness_landscape_sampling_papkou2023_tsv", "dhfr_papkou2023_zenodo_archive_metadata"],
            "publication": "Papkou et al. 2023; Wagner-lab compact TSV mirror",
            "doi": "10.1126/science.adh3860",
            "tier": "tier1_near_complete",
            "readout": "trimethoprim resistance fitness",
            "notes": "Use this as the primary DHFR/folA empirical landscape. The GraphFLA mirror is retained separately and has fewer rows.",
        },
    )


def normalize_proteingym(source_root: Path, processed_root: Path) -> list[dict[str, Any]]:
    reference_path = source_root / "proteingym/v1.3/reference/DMS_substitutions.csv"
    parquet_path = source_root / "proteingym/v1.3/DMS_substitutions.parquet"
    ref = pd.read_csv(reference_path)
    tsuboyama_ids = sorted(ref.loc[ref["DMS_id"].astype(str).str.contains("Tsuboyama_2023", na=False), "DMS_id"].astype(str))
    selected_ids = sorted(set(PROTEINGYM_ASSAY_IDS + tsuboyama_ids))
    assay_df = pd.read_parquet(parquet_path, filters=[("DMS_id", "in", selected_ids)])

    entries: list[dict[str, Any]] = []
    id_to_ref = ref.set_index("DMS_id", drop=False)
    assay_index_rows = []
    for dms_id, sub in assay_df.groupby("DMS_id", sort=True):
        sub = sub.copy()
        meta = id_to_ref.loc[dms_id].to_dict() if dms_id in id_to_ref.index else {}
        dataset_id = f"proteingym_{slug(str(dms_id))}"
        tier = "tier2_stability_core" if "Tsuboyama_2023" in str(dms_id) else "tier3_order1_or_sparse_multimutant"
        out = pd.DataFrame(
            {
                "dataset_id": dataset_id,
                "source_id": "proteingym_v1_3_substitutions",
                "tier": tier,
                "landscape_name": str(dms_id),
                "readout": str(meta.get("selection_type", "")),
                "sequence": sub["mutated_sequence"].astype(str),
                "target_sequence": sub["target_seq"].astype(str),
                "genotype": sub["mutant"].astype(str),
                "fitness": pd.to_numeric(sub["DMS_score"], errors="coerce"),
                "fitness_binary": pd.to_numeric(sub["DMS_score_bin"], errors="coerce"),
                "fitness_is_higher_better": True,
                "alphabet": "20-letter protein sequence",
                "source_row": sub.index.astype(int),
            }
        )
        entries.append(
            write_dataset(
                processed_root,
                dataset_id,
                out,
                {
                    "source_ids": ["proteingym_v1_3_reference", "proteingym_v1_3_substitutions"],
                    "publication": "ProteinGym v1.3 substitution benchmark",
                    "doi": "10.48550/arXiv.2307.01183",
                    "tier": tier,
                    "readout": str(meta.get("selection_type", "")),
                    "notes": f"Extracted assay {dms_id} from ProteinGym. Original reference metadata retained in assay index.",
                },
            )
        )
        assay_index_rows.append(
            {
                "dataset_id": dataset_id,
                "DMS_id": dms_id,
                "n_rows": len(out),
                "tier": tier,
                "selection_type": meta.get("selection_type"),
                "taxon": meta.get("taxon"),
                "DMS_number_single_mutants": meta.get("DMS_number_single_mutants"),
                "DMS_number_multiple_mutants": meta.get("DMS_number_multiple_mutants"),
            }
        )

    assay_index = pd.DataFrame(assay_index_rows).sort_values(["tier", "DMS_id"])
    index_dir = processed_root / "proteingym"
    index_dir.mkdir(parents=True, exist_ok=True)
    assay_index.to_csv(index_dir / "selected_assay_index.csv", index=False)
    return entries


def first_member(zip_path: Path, suffix: str) -> str:
    with zipfile.ZipFile(zip_path) as zf:
        matches = [name for name in zf.namelist() if name.lower().endswith(suffix.lower()) and not name.endswith("/")]
        if not matches:
            raise ValueError(f"No {suffix} file found in {zip_path}")
        return matches[0]


def read_first_csv_zip(zip_path: Path) -> pd.DataFrame:
    member = first_member(zip_path, ".csv")
    with zipfile.ZipFile(zip_path) as zf:
        with zf.open(member) as handle:
            return pd.read_csv(handle)


def normalize_flip(source_root: Path, processed_root: Path) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []

    gb1_zip = source_root / "flip/gb1/four_mutations_full_data.csv.zip"
    if gb1_zip.exists():
        df = read_first_csv_zip(gb1_zip)
        seq_col = "sequence" if "sequence" in df.columns else df.columns[0]
        fitness_col = "Fitness" if "Fitness" in df.columns else ("fitness" if "fitness" in df.columns else "Fittness")
        out = pd.DataFrame(
            {
                "dataset_id": "flip_gb1_four_mutations_full_data",
                "source_id": "flip_splits_gb1_four_mutations_full_data_csv_zip",
                "tier": "tier3_ml_split",
                "landscape_name": "FLIP GB1 four-mutation benchmark",
                "readout": "GB1 binding fitness",
                "sequence": df[seq_col].astype(str),
                "fitness": pd.to_numeric(df[fitness_col], errors="coerce"),
                "fitness_is_higher_better": True,
                "source_row": range(len(df)),
            }
        )
        for col in df.columns:
            if "split" in col.lower() or col.lower() in {"set", "stage"}:
                out[col] = df[col].astype(str)
        entries.append(write_dataset(processed_root, "flip_gb1_four_mutations_full_data", out, {"source_ids": ["flip_splits_gb1_four_mutations_full_data_csv_zip"], "publication": "Dallago et al. 2021 FLIP benchmark", "doi": "10.48550/arXiv.2111.06772", "tier": "tier3_ml_split", "readout": "GB1 binding fitness", "notes": "FLIP full GB1 data with any split columns preserved."}))

    aav_zip = source_root / "flip/aav/full_data.csv.zip"
    if aav_zip.exists():
        df = read_first_csv_zip(aav_zip)
        seq_col = "mutated_region" if "mutated_region" in df.columns else ("sequence" if "sequence" in df.columns else df.columns[0])
        fitness_col = "score" if "score" in df.columns else ("viral_selection" if "viral_selection" in df.columns else ("label" if "label" in df.columns else None))
        out = pd.DataFrame(
            {
                "dataset_id": "flip_aav_full_data",
                "source_id": "flip_splits_aav_full_data_csv_zip",
                "tier": "tier3_ml_split",
                "landscape_name": "FLIP AAV capsid benchmark",
                "readout": "AAV viability/selection",
                "sequence": df[seq_col].astype(str),
                "fitness": pd.to_numeric(df[fitness_col], errors="coerce") if fitness_col else pd.NA,
                "fitness_is_higher_better": True,
                "source_row": range(len(df)),
            }
        )
        for col in df.columns:
            if col not in out.columns and col.lower() in {"mutation_mask", "category", "is_viable", "split", "set"}:
                out[col] = df[col].astype(str)
        entries.append(write_dataset(processed_root, "flip_aav_full_data", out, {"source_ids": ["flip_splits_aav_full_data_csv_zip"], "publication": "Dallago et al. 2021 FLIP benchmark", "doi": "10.48550/arXiv.2111.06772", "tier": "tier3_ml_split", "readout": "AAV viability/selection", "notes": "FLIP full AAV data; raw split zip retained under data/external."}))

    meltome_zip = source_root / "flip/meltome/full_dataset.json.zip"
    if meltome_zip.exists():
        member = first_member(meltome_zip, ".json")
        with zipfile.ZipFile(meltome_zip) as zf:
            with zf.open(member) as handle:
                data = json.load(handle)
        rows = []
        for item in data:
            rows.append(
                {
                    "dataset_id": "flip_meltome_stability_summary",
                    "source_id": "flip_splits_meltome_full_dataset_json_zip",
                    "tier": "tier3_ml_split",
                    "landscape_name": "FLIP meltome thermostability benchmark",
                    "readout": "melting temperature",
                    "protein_id": item.get("proteinId"),
                    "uniprot_accession": item.get("uniprotAccession"),
                    "run_name": item.get("runName"),
                    "sequence": item.get("sequence"),
                    "melting_point": item.get("meltingPoint"),
                    "quant_norm_melting_point": item.get("quantNormMeltingPoint"),
                    "n_melting_curve_points": len(item.get("meltingBehaviour") or []),
                }
            )
        out = pd.DataFrame(rows)
        entries.append(write_dataset(processed_root, "flip_meltome_stability_summary", out, {"source_ids": ["flip_splits_meltome_full_dataset_json_zip"], "publication": "Dallago et al. 2021 FLIP benchmark", "doi": "10.48550/arXiv.2111.06772", "tier": "tier3_ml_split", "readout": "melting temperature", "notes": "Protein-level summary extracted from FLIP meltome JSON; raw curves remain in external zip."}))

    inventory_rows = []
    for zip_path in sorted((source_root / "flip").rglob("*.zip")):
        with zipfile.ZipFile(zip_path) as zf:
            for info in zf.infolist():
                if not info.is_dir():
                    inventory_rows.append({"zip_path": str(zip_path), "member": info.filename, "size_bytes": info.file_size})
    if inventory_rows:
        inv = pd.DataFrame(inventory_rows)
        inv_dir = processed_root / "flip"
        inv_dir.mkdir(parents=True, exist_ok=True)
        inv.to_csv(inv_dir / "zip_inventory.csv", index=False)

    return entries


def write_manifest(entries: list[dict[str, Any]], dataset_entries: list[dict[str, Any]]) -> None:
    EXTERNAL_ROOT.mkdir(parents=True, exist_ok=True)
    PROCESSED_ROOT.mkdir(parents=True, exist_ok=True)
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "external_root": str(EXTERNAL_ROOT),
        "processed_root": str(PROCESSED_ROOT),
        "scope": "Science plan part 1 ensemble landscape validation dataset preparation.",
        "record_boundary": "Generated as agent-side data preparation; no labproj record, registry, or claim files are modified.",
        "sources": entries,
        "processed_datasets": dataset_entries,
    }
    with (EXTERNAL_ROOT / "manifest.yaml").open("w") as handle:
        yaml.safe_dump(manifest, handle, sort_keys=False)
    with (PROCESSED_ROOT / "dataset_index.csv").open("w") as handle:
        pd.DataFrame(dataset_entries).to_csv(handle, index=False)
    with (PROCESSED_ROOT / "manifest.yaml").open("w") as handle:
        yaml.safe_dump(manifest, handle, sort_keys=False)


def write_readme(dataset_entries: list[dict[str, Any]]) -> None:
    lines = [
        "# Ensemble Landscape Validation Datasets",
        "",
        "Generated by `scripts/prepare_ensemble_landscape_validation_datasets.py`.",
        "This bundle is agent-side data preparation for the scratch science plan; it is not a labproj record entry.",
        "",
        "## Contents",
        "",
        "- `datasets/<dataset_id>/variants.csv`: portable CSV normalized variant table.",
        "- `datasets/<dataset_id>/variants.parquet`: same table in columnar form for descriptor pipelines.",
        "- `datasets/<dataset_id>/metadata.yaml`: per-dataset provenance and schema.",
        "- `dataset_index.csv`: row counts and file paths across all processed datasets.",
        "- `manifest.yaml`: source download provenance and checksums.",
        "",
        "## Normalized Core Columns",
        "",
        "`dataset_id`, `source_id`, `tier`, `landscape_name`, `readout`, `sequence`, `genotype`, `fitness`, and `fitness_is_higher_better` are standardized where the source permits it. Source-specific columns are retained when they are needed for coordinates or split labels.",
        "",
        "## Dataset Count",
        "",
        f"Processed datasets: {len(dataset_entries)}",
        "",
    ]
    PROCESSED_ROOT.mkdir(parents=True, exist_ok=True)
    (PROCESSED_ROOT / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="Re-download files even when they already exist.")
    parser.add_argument("--skip-download", action="store_true", help="Only process files already present under data/external.")
    args = parser.parse_args()

    EXTERNAL_ROOT.mkdir(parents=True, exist_ok=True)
    PROCESSED_ROOT.mkdir(parents=True, exist_ok=True)

    source_entries: list[dict[str, Any]] = []
    if args.skip_download:
        for spec in SOURCES:
            path = EXTERNAL_ROOT / spec.rel_path
            if path.exists():
                source_entries.append(file_entry(spec, path, status="existing"))
            else:
                source_entries.append({"source_id": spec.source_id, "status": "missing", "path": str(path), "url": spec.url})
    else:
        for spec in SOURCES:
            print(f"retrieving {spec.source_id}")
            source_entries.append(download(spec, EXTERNAL_ROOT, force=args.force))

    dataset_entries: list[dict[str, Any]] = []
    for spec in GRAPHFLA_DATASETS:
        print(f"processing {spec.dataset_id}")
        dataset_entries.append(normalize_graphfla(spec, EXTERNAL_ROOT, PROCESSED_ROOT))

    print("processing TrpB full-sequence mirror")
    dataset_entries.append(normalize_hf_trpb(EXTERNAL_ROOT, PROCESSED_ROOT))

    print("processing Papkou DHFR full TSV")
    dataset_entries.append(normalize_papkou_full_tsv(EXTERNAL_ROOT, PROCESSED_ROOT))

    print("processing ProteinGym selected assays")
    dataset_entries.extend(normalize_proteingym(EXTERNAL_ROOT, PROCESSED_ROOT))

    print("processing FLIP datasets")
    dataset_entries.extend(normalize_flip(EXTERNAL_ROOT, PROCESSED_ROOT))

    write_manifest(source_entries, dataset_entries)
    write_readme(dataset_entries)
    print(f"wrote {len(dataset_entries)} processed datasets to {PROCESSED_ROOT}")


if __name__ == "__main__":
    main()
