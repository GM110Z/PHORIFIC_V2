#!/usr/bin/env python3
"""Create the MGE hotspot view used for downstream annotation."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path


STANDARD_FIELDS = [
    "hotspot_id", "locus_id", "genome", "contig", "start", "end",
    "left_boundary_id", "right_boundary_id", "left_boundary_cluster",
    "right_boundary_cluster", "defence_system", "defence_present",
    "architecture_members", "architecture_type", "left_anchor_role", "right_anchor_role",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare MGE hotspot tables for downstream annotation")
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--kind", choices=("plasmid", "prophage"), required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--include-edge", action="store_true",
                        help="include experimental one-edge hotspots (excluded by default)")
    args = parser.parse_args()

    results = args.results_dir.resolve()
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    standard = read_tsv(results / "hotspots.tsv")
    edge = (read_tsv(results / f"{args.kind}_edge_hotspot.tsv")
            if args.include_edge else [])
    merged = []
    for row in standard:
        copy = dict(row)
        copy.update({
            "architecture_type": "two_sided",
            "left_anchor_role": "non_defence_boundary",
            "right_anchor_role": "non_defence_boundary",
        })
        merged.append(copy)
    for row in edge:
        edge_side = row["edge_anchor_side"]
        terminal_id = row.get("terminal_edge_gene_id", row.get("defence_anchor_id", ""))
        terminal_cluster = row.get("terminal_edge_gene_cluster", row.get("defence_anchor_cluster", ""))
        if edge_side == "right":
            left_id, right_id = row["conserved_boundary_id"], terminal_id
            left_cluster, right_cluster = row["conserved_boundary_cluster"], terminal_cluster
            left_role, right_role = "non_defence_boundary", "open_element_edge"
        else:
            left_id, right_id = terminal_id, row["conserved_boundary_id"]
            left_cluster, right_cluster = terminal_cluster, row["conserved_boundary_cluster"]
            left_role, right_role = "open_element_edge", "non_defence_boundary"
        merged.append({
            "hotspot_id": row["hotspot_id"], "locus_id": row["locus_id"],
            "genome": row["mge_id"], "contig": row["mge_id"],
            "start": row["start"], "end": row["end"],
            "left_boundary_id": left_id, "right_boundary_id": right_id,
            "left_boundary_cluster": left_cluster, "right_boundary_cluster": right_cluster,
            "defence_system": row["defence_system"],
            "defence_present": "0" if row["defence_system"].lower() == "no defence detected" else "1",
            "architecture_members": row["architecture_members"],
            "architecture_type": "one_sided_edge",
            "left_anchor_role": left_role, "right_anchor_role": right_role,
        })
    with (out / "hotspots.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=STANDARD_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in STANDARD_FIELDS} for row in merged)
    target = results / "loci"
    link = out / "loci"
    if link.is_symlink() and link.resolve() != target:
        link.unlink()
    elif link.exists() and not link.is_symlink():
        raise SystemExit(f"error: refusing to replace non-symlink path: {link}")
    if not link.exists():
        os.symlink(target, link, target_is_directory=True)
    print(f"Prepared {len(standard)} standard hotspot instances"
          f" and {len(edge)} optional edge instances in {out}")


if __name__ == "__main__":
    main()
