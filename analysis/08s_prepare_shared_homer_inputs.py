#!/usr/bin/env python3
"""Prepare exact matched target/background BEDs for two shared-loss HOMER runs.

Usage: python3 analysis/08s_prepare_shared_homer_inputs.py \
  /path/to/08n_wt_ISC_interval_loss_pathways \
  /path/to/08l_wt_ISC_open_close_motifs \
  /path/to/output_directory
"""
import csv
import gzip
import json
import sys
from pathlib import Path


def fail(message):
    raise SystemExit(message)


def read_csv(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle))


def read_bed(path):
    records = {}
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                fail(f"Expected four-column BED: {path}")
            chrom, start, end, feature = parts[:4]
            if feature in records or int(start) < 0 or int(end) <= int(start):
                fail(f"Duplicate/invalid BED interval {feature}: {path}")
            records[feature] = (chrom, int(start), int(end), feature)
    return records


def write_bed(path, rows):
    with path.open("w") as handle:
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")


if len(sys.argv) != 4:
    fail("Usage: 08s_prepare_shared_homer_inputs...py STEP08N_DIR STEP08L_DIR OUT_DIR")
step_n, step_l, out = (Path(x) for x in sys.argv[1:])
out.mkdir(parents=True, exist_ok=True)
pair_file = step_l / "02_all_622_target_neutral_matches.csv"
universe_file = step_l / "02_evaluable_WT_universe_sequence_covariates.csv.gz"
if not pair_file.is_file() or not universe_file.is_file():
    fail("The Step 08L matched-pair or evaluable-universe file is missing.")
pairs = read_csv(pair_file)
universe = read_csv(universe_file)
pair_by_target = {row["target_feature"]: row for row in pairs}
universe_by_feature = {row["feature"]: row for row in universe}
if len(pairs) != 622 or len(pair_by_target) != 622 or len(universe_by_feature) != 2028:
    fail("Frozen 08L matching/universe identity changed.")
manifest = []
for direction in ("open", "closed"):
    key = f"shared_HFD_{direction}_lost"
    target_bed = read_bed(step_n / "sets" / f"{key}.bed")
    expected_direction = f"HFD_{direction}"
    matched_targets, matched_controls, paired_rows = [], [], []
    unmatched = []
    for feature, original_bed in target_bed.items():
        pair = pair_by_target.get(feature)
        if pair is None or pair["WT_direction"] != expected_direction:
            fail(f"08N target {feature} does not match the 08L direction.")
        t = universe_by_feature[feature]
        expected_bed = (t["chrom"], int(t["start_1based"]) - 1,
                        int(t["end_1based"]), feature)
        if original_bed != expected_bed:
            fail(f"08N/08L coordinates differ at {feature}.")
        control = pair["background_feature"]
        if not control or control == "NA":
            unmatched.append(feature)
            continue
        b = universe_by_feature.get(control)
        if b is None or b["membership"] != "neutral":
            fail(f"08L matched control is absent/non-neutral: {control}")
        if pair["target_context"] != pair["background_context"] or \
                pair["target_context"] != t["context"] or \
                pair["background_context"] != b["context"]:
            fail(f"Genomic contexts disagree for {feature}/{control}.")
        matched_targets.append(original_bed)
        matched_controls.append((b["chrom"], int(b["start_1based"]) - 1,
                                 int(b["end_1based"]), control))
        paired_rows.append({"target_feature": feature,
                            "background_feature": control,
                            "WT_direction": expected_direction,
                            "context": t["context"]})
    if len(matched_targets) < 10 or len({x[3] for x in matched_controls}) != len(matched_controls):
        fail(f"Insufficient/duplicated matched controls for {key}.")
    if len(matched_targets) + len(unmatched) != len(target_bed):
        fail(f"Matched + unmatched target accounting failed for {key}.")
    write_bed(out / f"{key}_matched_targets.bed", matched_targets)
    write_bed(out / f"{key}_matched_neutral_background.bed", matched_controls)
    with (out / f"{key}_paired_intervals.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=paired_rows[0].keys())
        writer.writeheader()
        writer.writerows(paired_rows)
    manifest.append({"set": key, "original_targets": len(target_bed),
                     "matched_targets": len(matched_targets),
                     "matched_controls": len(matched_controls),
                     "unmatched_excluded": len(unmatched),
                     "unmatched_features": unmatched})
with (out / "00_matching_audit.json").open("w") as handle:
    json.dump(manifest, handle, indent=2)
    handle.write("\n")
for row in manifest:
    print(row["set"], "targets", row["original_targets"], "matched",
          row["matched_targets"], "unmatched", row["unmatched_excluded"])
