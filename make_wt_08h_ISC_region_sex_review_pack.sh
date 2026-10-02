#!/usr/bin/env bash
# Package Step 08H's frozen WT peak list and sex-marker sensitivity audit.
set -euo pipefail

usage() {
  cat >&2 <<'EOF'
Usage: bash make_wt_08h_ISC_region_sex_review_pack.sh [WT_OUTPUT_DIR] [PACK_DIR]

Run from the repository root after finishing Step 08H. Defaults:
  WT_OUTPUT_DIR=/scratch/dsaiz/Yesenia_scData2026/results/v2_resequenced/independent/wt
  PACK_DIR=WT_OUTPUT_DIR/review_packs

The 08H notebook may live in analysis/ or in the earlier deliverables bundle.
EOF
}
if (( $# > 2 )); then usage; exit 2; fi
if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then usage; exit 0; fi
command -v tar >/dev/null || { echo 'tar is not on PATH.' >&2; exit 1; }

wt_root="${1:-/scratch/dsaiz/Yesenia_scData2026/results/v2_resequenced/independent/wt}"
pack_dir="${2:-$wt_root/review_packs}"
step_dir="$wt_root/08h_wt_ISC_region_freeze_sex_audit"
prior_dir="$wt_root/08g_C0_C2_and_C0W0_ISC"
repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
notebook="$repo_root/analysis/08h_freeze_wt_ISC_candidate_regions_and_audit_sex.Rmd"
if [[ ! -s "$notebook" ]]; then
  notebook="$repo_root/deliverables/wt_crypt_followup/analysis/08h_freeze_wt_ISC_candidate_regions_and_audit_sex.Rmd"
fi

[[ -s "$step_dir/STEP_08H_COMPLETE.txt" ]] || {
  echo "Step 08H has not completed: $step_dir/STEP_08H_COMPLETE.txt" >&2; exit 1;
}
[[ -s "$prior_dir/STEP_08G_COMPLETE.txt" ]] || {
  echo "Missing completed Step 08G: $prior_dir/STEP_08G_COMPLETE.txt" >&2; exit 1;
}
[[ -s "$notebook" ]] || {
  echo "Cannot find the 08H provenance notebook in analysis/ or deliverables/." >&2; exit 1;
}

# Check every essential input before making a staging area or output archive.
required_08h=(
  STEP_08H_COMPLETE.txt
  input_manifest.csv
  sessionInfo.txt
  00_WT_pooled_mouse_sex_composition.csv
  00_matrix_alignment.csv
  01_sex_marker_and_depth_balance_by_group_diet.csv
  01_sex_marker_proxy_strata_by_group_diet.csv
  01_nucleus_level_sex_proxy_and_gate_audit.csv.gz
  02_all_622_top15_ATAC_peaks_with_evidence_flags.csv.gz
  02_C0_supported_37_peak_candidates_full_evidence.csv
  02_higher_count_C0_supported_sensitivity.csv
  02_six_nested_C0_and_top10_HFD_open_peaks.csv
  02_candidate_support_summary.csv
  03_Xist_stratum_eligibility.csv
)
for item in "${required_08h[@]}"; do
  [[ -s "$step_dir/$item" ]] || {
    echo "Missing/empty required 08H output: $step_dir/$item" >&2; exit 1;
  }
done
for item in 02_C0_supported_37_GRCm39.bed \
            02_higher_count_C0_supported_GRCm39.bed \
            03_C0_supported_region_effects_by_Xist_detection.csv; do
  [[ -f "$step_dir/$item" ]] || {
    echo "Missing 08H output: $step_dir/$item" >&2; exit 1;
  }
done
[[ -s "$step_dir/plots/01_Xist_detection_by_group_diet.pdf" ]] || {
  echo "Missing 08H plot: $step_dir/plots/01_Xist_detection_by_group_diet.pdf" >&2; exit 1;
}
required_08g=(
  01_core_marker_detection_by_state_diet.csv
  02_C0W0_gate_size_and_score_sensitivity.csv
  02_C0W0_selected_identity_and_depth_by_diet.csv
  03_contrast_definitions_and_diet_balance.csv
  04_descriptive_effect_counts.csv
  05_candidate_peak_support_tiers.csv
  06_C0W0_top15_candidates_with_depth_check.csv.gz
)
for item in "${required_08g[@]}"; do
  [[ -s "$prior_dir/$item" ]] || {
    echo "Missing/empty 08G context: $prior_dir/$item" >&2; exit 1;
  }
done

mkdir -p -- "$pack_dir"
stage="$(mktemp -d "$pack_dir/.wt08h_review.XXXXXXXX")"
temp_archive="$(mktemp "$pack_dir/.wt08h_archive.XXXXXXXX")"
trap 'rm -rf -- "$stage"; rm -f -- "$temp_archive"' EXIT
content="$stage/wt_08h_ISC_region_sex_review"
mkdir -p -- "$content/08h/plots" "$content/08g_context" "$content/scripts"

for item in "${required_08h[@]}" \
            02_C0_supported_37_GRCm39.bed \
            02_higher_count_C0_supported_GRCm39.bed \
            03_C0_supported_region_effects_by_Xist_detection.csv; do
  cp -- "$step_dir/$item" "$content/08h/$item"
done
# This summary is conditional: 08H writes it only if a stratum was eligible.
if [[ -f "$step_dir/03_Xist_stratum_peak_summary.csv" ]]; then
  cp -- "$step_dir/03_Xist_stratum_peak_summary.csv" "$content/08h/"
fi
cp -- "$step_dir/plots/01_Xist_detection_by_group_diet.pdf" "$content/08h/plots/"
for item in "${required_08g[@]}"; do
  cp -- "$prior_dir/$item" "$content/08g_context/$item"
done
cp -- "$notebook" "$content/scripts/"

cat > "$content/README_review_first.md" <<'EOF'
# WT ISC-enriched candidate peaks and sex-marker review

1. Read `08h/00_WT_pooled_mouse_sex_composition.csv`, then
   `08h/01_sex_marker_and_depth_balance_by_group_diet.csv` and
   `08h/plots/01_Xist_detection_by_group_diet.pdf`. Compare CON and HFD
   Xist/Y detection and RNA/ATAC depth within C0_all and C0W0_top15.
2. Read `08h/01_sex_marker_proxy_strata_by_group_diet.csv` and
   `08h/03_Xist_stratum_eligibility.csv`. Xist non-detection is not a male
   call. The per-nucleus data are in `08h/01_nucleus_level_sex_proxy_and_gate_audit.csv.gz`.
3. Inspect `08h/02_candidate_support_summary.csv`, then the full evidence
   for the 37 C0-supported regions in
   `08h/02_C0_supported_37_peak_candidates_full_evidence.csv`.
   The entire 622-peak top-15 set is also preserved. The higher-count and
   nested-top10 files are sensitivity subsets; their BEDs give exact GRCm39
   intervals for later cohort comparison.
4. Read `08h/03_C0_supported_region_effects_by_Xist_detection.csv` and,
   when present, `08h/03_Xist_stratum_peak_summary.csv`. Compare each
   region's effect sign, counts, detection and evaluability within the
   eligible Xist strata. Check depth and strata sizes before interpreting
   changes in sign. A missing summary means no stratum produced rows.
5. Use `08g_context/` for the original marker gate, group sizes and prior
   ATAC effect definitions. `scripts/` contains the source 08H notebook.

These are pooled descriptive effects from one WT library per diet. Xist and
Y transcript detection are imperfect proxies; stratification is a
sensitivity diagnostic and cannot adjust for mouse sex or yield biological
replicate-level significance. The list was selected in WT, so agreement
within overlapping subsets is not independent validation.
EOF

stamp="$(date -u +%Y-%m-%d_%H%M%S_UTC)"
archive="$pack_dir/wt_08h_ISC_region_sex_review_${stamp}.tar.gz"
[[ ! -e "$archive" ]] || { echo "Archive already exists: $archive" >&2; exit 1; }
tar -C "$stage" -czf "$temp_archive" wt_08h_ISC_region_sex_review
mv -- "$temp_archive" "$archive"
if command -v sha256sum >/dev/null; then sha256sum "$archive"; fi
echo "Review pack: $archive"
echo 'Start with: wt_08h_ISC_region_sex_review/README_review_first.md'
