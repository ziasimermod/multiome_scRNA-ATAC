#!/usr/bin/env python3
"""Prepare 08L WT HFD-open/closed candidate BEDs and their frozen neutral pairs.

Usage: python3 analysis/08s_prepare_wt_622_homer_inputs.py \
  STEP08L_DIR OUTPUT_INPUTS_DIR

The candidate universe contains 622 peaks total (387 opening, 235 closing).
All 622 are HOMER targets. Their 08L distinct neutral controls form the
background; 14 targets have no qualifying matched control, recorded in audit.
"""
import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path


def require(test, message):
    if not test:
        raise SystemExit(message)


def rows(path):
    require(path.is_file(), f"Missing Step 08L input: {path}")
    with (gzip.open(path, "rt", newline="") if path.suffix == ".gz"
          else path.open(newline="")) as handle:
        return list(csv.DictReader(handle))


def digest(path):
    checksum = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            checksum.update(block)
    return checksum.hexdigest()


def bed(row):
    start, end = int(row["start_1based"]), int(row["end_1based"])
    require(start >= 1 and end >= start,
            f"Invalid 1-based peak coordinates: {row['feature']}")
    return (row["chrom"], start - 1, end, row["feature"])


def write_bed(path, records):
    with path.open("w") as handle:
        for record in records:
            handle.write("\t".join(map(str, record)) + "\n")


require(len(sys.argv) == 3,
        "Usage: 08s_prepare_wt_622_homer_inputs.py "
        "STEP08L_DIR OUTPUT_INPUTS_DIR")
source, output = Path(sys.argv[1]), Path(sys.argv[2])
candidate_file = source / "01_all_622_open_close_genotype_and_support.csv.gz"
pairs_file = source / "02_all_622_target_neutral_matches.csv"
universe_file = source / "02_evaluable_WT_universe_sequence_covariates.csv.gz"
candidate_rows = rows(candidate_file)
pair_rows = rows(pairs_file)
universe_rows = rows(universe_file)
require(len(candidate_rows) == len(pair_rows) == 622 and
        len(universe_rows) == 2028, "Frozen 08L candidate/pair/universe counts changed")
candidates = {r["feature"]: r for r in candidate_rows}
pairs = {r["target_feature"]: r for r in pair_rows}
universe = {r["feature"]: r for r in universe_rows}
require(len(candidates) == len(pairs) == 622 and len(universe) == 2028 and
        candidates.keys() == pairs.keys(),
        "Duplicate or missing 08L candidate, pair, or universe IDs")
require({r["WT_direction"] for r in candidate_rows} ==
        {"HFD_open", "HFD_closed"}, "Unexpected WT HFD direction")
require(sum(r["WT_direction"] == "HFD_open" for r in candidate_rows) == 387 and
        sum(r["WT_direction"] == "HFD_closed" for r in candidate_rows) == 235,
        "The frozen 387 opening / 235 closing partition changed")
output.mkdir(parents=True, exist_ok=True)
audit = []
all_control_ids = set()
for direction in ("open", "closed"):
    key = f"WT_HFD_{direction}"
    selected = [r for r in candidate_rows if r["WT_direction"] == f"HFD_{direction}"]
    all_targets, matched_targets, controls, pair_records, missing = [], [], [], [], []
    for row in selected:
        feature = row["feature"]
        pair = pairs[feature]
        target = universe.get(feature)
        require(target is not None and pair["WT_direction"] == f"HFD_{direction}",
                f"Lost WT direction or universe record: {feature}")
        require(bed(row) == bed(target) and
                target["membership"] == f"HFD_{direction}",
                f"WT candidate BED or membership differs from 08L: {feature}")
        all_targets.append(bed(target))
        background_feature = pair["background_feature"]
        if not background_feature or background_feature == "NA":
            missing.append(feature)
            continue
        background = universe.get(background_feature)
        require(background is not None and background["membership"] == "neutral",
                f"Missing/non-neutral matched background: {feature}")
        require(pair["target_context"] == target["context"] ==
                pair["background_context"] == background["context"],
                f"08L target/control genomic context changed: {feature}")
        require(background_feature not in all_control_ids,
                f"A background peak was reused: {background_feature}")
        all_control_ids.add(background_feature)
        matched_targets.append(bed(target))
        controls.append(bed(background))
        pair_records.append((feature, background_feature, f"HFD_{direction}",
                             target["context"]))
    require(len(all_targets) == len(selected) and
            len(matched_targets) + len(missing) == len(selected) and
            len(matched_targets) > 10 and len(matched_targets) == len(controls),
            f"Matching audit failed: {key}")
    write_bed(output / f"{key}_all_targets.bed", all_targets)
    write_bed(output / f"{key}_matched_targets.bed", matched_targets)
    write_bed(output / f"{key}_matched_neutral_background.bed", controls)
    with (output / f"{key}_paired_intervals.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("target_feature", "background_feature", "WT_direction",
                         "context"))
        writer.writerows(pair_records)
    audit.append({"set": key, "WT_candidates": len(selected),
                  "all_targets_in_primary_HOMER": len(all_targets),
                  "matched_target_control_pairs": len(matched_targets),
                  "targets_without_matched_control": len(missing),
                  "unmatched_features": missing})
require(sum(a["WT_candidates"] for a in audit) == 622 and
        len(all_control_ids) == sum(a["matched_target_control_pairs"] for a in audit),
        "The two direction groups no longer partition the 622 WT candidates")
with (output / "00_matching_audit.json").open("w") as handle:
    json.dump({"source_files": {
        str(p): digest(p) for p in (candidate_file, pairs_file, universe_file)},
        "sets": audit}, handle, indent=2)
    handle.write("\n")
for group in audit:
    print(group["set"], "WT candidates", group["WT_candidates"],
          "matched pairs", group["matched_target_control_pairs"],
          "targets lacking matched control", group["targets_without_matched_control"])
