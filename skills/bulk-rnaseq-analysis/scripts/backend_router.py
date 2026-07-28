#!/usr/bin/env python3
"""Discover and route between the RNAseq-Templates and TCGA backends."""

from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path
from typing import Iterable


RNASEQ_SENTINELS = (
    "RNAseq_lib",
    "notebooks",
    "references/TEMPLATE_SELECTION.md",
)
TCGA_SENTINELS = (
    "tcga_toolkit/VERSION",
    "tcga_toolkit/scripts/run_task.R",
    "tcga_toolkit/specs/README.md",
)

TCGA_SOURCES = {"tcga", "target", "gtex"}
TCGA_ANALYSES = {
    "audit-data",
    "clinical-association",
    "cnv",
    "cohort-qc",
    "drug-response",
    "external-validation",
    "gtex-compare",
    "immune-phenotype",
    "maf",
    "methylation",
    "mutation",
    "mutation-survival",
    "pan-cancer",
    "prepare-bulk-rna",
    "prognostic-model",
    "stage",
    "subtype",
    "survival-map",
    "tmb",
}
GENERIC_ANALYSES = {
    "deg",
    "deseq2",
    "enrichment",
    "geo",
    "gsea",
    "gsva",
    "limma-voom",
    "ora",
    "report",
    "survival",
    "time-course",
    "tme",
    "visualization",
    "wgcna",
}


def unique_paths(paths: Iterable[Path]) -> list[Path]:
    seen: set[str] = set()
    answer: list[Path] = []
    for path in paths:
        resolved = path.expanduser().resolve()
        key = str(resolved)
        if key not in seen:
            seen.add(key)
            answer.append(resolved)
    return answer


def candidates(repo_name: str, env_name: str, start: Path) -> list[Path]:
    values: list[Path] = []
    configured = os.environ.get(env_name)
    if configured:
        configured_path = Path(configured)
        values.append(configured_path)
        if repo_name == "TCGA" and configured_path.name == "tcga_toolkit":
            values.append(configured_path.parent)

    start = start.resolve()
    for parent in (start, *start.parents):
        if parent.name == repo_name:
            values.append(parent)
        values.append(parent / repo_name)

    home = Path.home()
    values.extend(
        (
            home / "Projects" / repo_name,
            home / "Documents" / "Projects" / repo_name,
            home
            / "Library"
            / "Mobile Documents"
            / "com~apple~CloudDocs"
            / "Projects"
            / repo_name,
        )
    )
    return unique_paths(values)


def validate_root(path: Path, sentinels: tuple[str, ...]) -> bool:
    return path.is_dir() and all((path / sentinel).exists() for sentinel in sentinels)


def discover_one(
    repo_name: str, env_name: str, sentinels: tuple[str, ...], start: Path
) -> dict[str, object]:
    checked = candidates(repo_name, env_name, start)
    found = next((path for path in checked if validate_root(path, sentinels)), None)
    return {
        "found": found is not None,
        "root": str(found) if found else None,
        "configured_by": env_name if os.environ.get(env_name) else "auto-discovery",
        "required_sentinels": list(sentinels),
        "checked": [str(path) for path in checked],
    }


def repository_version(root: str | None, backend: str) -> str | None:
    if not root:
        return None
    repo = Path(root)
    if backend == "tcga-toolkit":
        version_file = repo / "tcga_toolkit" / "VERSION"
        return version_file.read_text(encoding="utf-8").strip()

    readme = repo / "README.md"
    if readme.exists():
        match = re.search(r"version-([0-9][0-9.]*)-", readme.read_text(encoding="utf-8"))
        if match:
            return match.group(1)
    return None


def discover(start: Path) -> dict[str, object]:
    rnaseq = discover_one(
        "RNAseq-Templates", "RNASEQ_TEMPLATES_ROOT", RNASEQ_SENTINELS, start
    )
    tcga = discover_one("TCGA", "TCGA_TOOLKIT_ROOT", TCGA_SENTINELS, start)
    rnaseq["version"] = repository_version(rnaseq["root"], "rnaseq-templates")
    tcga["version"] = repository_version(tcga["root"], "tcga-toolkit")
    return {"rnaseq-templates": rnaseq, "tcga-toolkit": tcga}


def normalize_analyses(value: str) -> set[str]:
    return {
        item.strip().lower().replace("_", "-").replace(" ", "-")
        for item in value.split(",")
        if item.strip()
    }


def route(source: str, input_type: str, analyses: set[str]) -> dict[str, object]:
    source = source.lower()
    input_type = input_type.lower()
    tcga_hits = sorted(analyses & TCGA_ANALYSES)
    generic_hits = sorted(analyses & GENERIC_ANALYSES)
    warnings: list[str] = []

    if {"deseq2", "deg"} & analyses and input_type in {
        "tpm",
        "vst",
        "rlog",
        "normalized",
    }:
        warnings.append(
            "DESeq2 requires integer-like raw counts; select a compatible method "
            "or provide raw counts."
        )
    if "tme" in analyses and input_type in {"vst", "rlog"}:
        warnings.append(
            "TME input cannot use VST/rlog as TPM; provide TPM or raw counts plus "
            "gene lengths."
        )

    if source in TCGA_SOURCES:
        backends = ["tcga-toolkit"]
        reason = f"{source.upper()} is handled by the cohort-aware TCGA toolkit."
    elif tcga_hits and generic_hits:
        backends = ["rnaseq-templates", "tcga-toolkit"]
        reason = (
            "The request combines generic expression analysis with "
            f"TCGA-specialized tasks: {', '.join(tcga_hits)}."
        )
    elif tcga_hits:
        backends = ["tcga-toolkit"]
        reason = f"Specialized cancer-cohort tasks require the TCGA toolkit: {', '.join(tcga_hits)}."
    elif source in {"local", "geo", "external"} or generic_hits:
        backends = ["rnaseq-templates"]
        reason = "The request is a generic matrix-based bulk RNA-seq workflow."
    else:
        backends = []
        reason = "Source and requested analyses are insufficient for deterministic routing."

    return {
        "primary_backend": backends[0] if backends else None,
        "backends": backends,
        "source": source,
        "input_type": input_type,
        "analyses": sorted(analyses),
        "reason": reason,
        "warnings": warnings,
        "needs_clarification": not backends,
    }


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    commands = root.add_subparsers(dest="command", required=True)

    discover_parser = commands.add_parser("discover", help="Find both repositories")
    discover_parser.add_argument("--start", type=Path, default=Path.cwd())
    discover_parser.add_argument("--json", action="store_true")

    route_parser = commands.add_parser("route", help="Select one or both backends")
    route_parser.add_argument(
        "--source",
        default="unknown",
        choices=("local", "geo", "tcga", "target", "gtex", "external", "unknown"),
    )
    route_parser.add_argument(
        "--input-type",
        default="unknown",
        choices=(
            "raw-counts",
            "tpm",
            "vst",
            "rlog",
            "normalized",
            "maf",
            "cnv",
            "methylation",
            "clinical",
            "unknown",
        ),
    )
    route_parser.add_argument("--analyses", default="")
    route_parser.add_argument("--json", action="store_true")
    return root


def emit(payload: dict[str, object], as_json: bool) -> None:
    if as_json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
        return
    for key, value in payload.items():
        print(f"{key}: {value}")


def main() -> int:
    args = parser().parse_args()
    if args.command == "discover":
        emit(discover(args.start), args.json)
        return 0
    result = route(args.source, args.input_type, normalize_analyses(args.analyses))
    emit(result, args.json)
    return 2 if result["needs_clarification"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
