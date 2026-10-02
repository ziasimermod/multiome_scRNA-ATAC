#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: bash make_wt_crypt_followup_review_pack.sh WT_OUTPUT_DIR [PACK_DIR]" >&2
  exit 2
fi

wt_root="${1%/}"
pack_dir="${2:-$PWD}"
release='wt_crypt_followup'
review_dir="05_wnn/$release/05f_review"
effect_dir='08e_refined_wt_crypt_diet'
need_files=(
  "$review_dir/STEP_05F_COMPLETE.txt"
  "$review_dir/input_manifest.csv"
  "$review_dir/refined_working_labels.csv"
  "$review_dir/matrix_cell_alignment_audit.csv"
  "$review_dir/refined_cell_assignments.csv.gz"
  "$review_dir/focused_composition_by_diet.csv"
  "$review_dir/qc_and_phase_by_working_unit_diet.csv"
  "$review_dir/r0.3_graph_and_boundary_stability.csv"
  "$review_dir/focused_reduction_depth_correlations.csv"
  "$review_dir/fetal_repair_sentinel_feature_availability.csv"
  "$review_dir/sentinel_detection_by_cluster_diet.csv"
  "$review_dir/fetal_revival_sentinel_co_detection.csv"
  "$review_dir/marker_program_audit_summary.csv"
  "$review_dir/RNA_only_r0.3_positive_markers.csv.gz"
  "$review_dir/WNN_r0.3_positive_markers_separately_by_diet.csv.gz"
  "$review_dir/targeted_promoter_evidence_with_availability.csv"
  "$review_dir/plots/01_fetal_repair_sentinel_detection.pdf"
  "$effect_dir/STEP_08E_COMPLETE.txt"
  "$effect_dir/input_manifest.csv"
  "$effect_dir/matrix_cell_alignment_audit.csv"
  "$effect_dir/refined_diet_contrast_definition_and_counts.csv"
  "$effect_dir/persistent_reference_gene_canonicalization.csv"
  "$effect_dir/RNA_top_descriptive_effects.csv"
  "$effect_dir/ATAC_top_descriptive_effects.csv"
  "$effect_dir/persistent_genes_refined_WT_concordance.csv"
  "$effect_dir/persistent_gene_concordance_counts.csv"
  "$effect_dir/equal_cell_and_depth_matched_resampling.csv.gz"
  "$effect_dir/phase_stratum_eligibility.csv"
  "$effect_dir/RNA_equal_phase_weight_effect_sensitivity.csv"
  "$effect_dir/target_gene_promoter_overlapping_peak_map.csv.gz"
  "$effect_dir/target_gene_promoter_assay_availability.csv"
  "$effect_dir/refined_C0_C2_RNA_promoter_peak_evidence.csv.gz"
  "$effect_dir/plots/01_source_d21_vs_refined_WT_d21_RNA.pdf"
  "$effect_dir/plots/02_refined_RNA_depth_sensitivity.pdf"
)

[[ -d "$wt_root" ]] || { echo "WT output directory absent: $wt_root" >&2; exit 1; }
for relative in "${need_files[@]}"; do
  [[ -s "$wt_root/$relative" ]] || {
    echo "Missing required review file: $wt_root/$relative" >&2
    exit 1
  }
done
mkdir -p "$pack_dir"
manifest_file="$(mktemp)"
trap 'rm -f "$manifest_file"' EXIT
printf '%s\n' "${need_files[@]}" > "$manifest_file"
stamp="$(date +%Y-%m-%d_%H%M%S)"
archive="$pack_dir/wt_crypt_refined_followup_review_v3_${stamp}.tar.gz"
tar -C "$wt_root" -czf "$archive" -T "$manifest_file"
sha256sum "$archive"
echo "Review pack: $archive"
