#!/usr/bin/env bash
# Run on SOL from any directory. Creates one 05D/05E/05F/08E review archive.
set -euo pipefail

wt_output='/scratch/dsaiz/Yesenia_scData2026/results/v2_resequenced/independent/wt'
pack_dir="$wt_output/review_packs"
[[ -d "$wt_output" ]] || { printf 'WT output directory missing: %s\n' "$wt_output" >&2; exit 1; }

step5d='05_wnn/gourab_d21_rna_reference_mapping/05d_symmetric_pseudobulk_validation'
step5e='05_wnn/focused_wt_crypt_stem_progenitor'
step5f='05_wnn/wt_crypt_followup/05f_review'
step8e='08e_refined_wt_crypt_diet'

required=(
  "$step5d/STEP_05D_COMPLETE.txt"
  "$step5d/tables/primary_similarity_with_bootstrap_support.csv"
  "$step5d/tables/primary_mutual_top1_matches.csv"
  "$step5d/tables/focused_WT_crypt_secretory_primary_matches.csv"
  "$step5d/tables/curated_lineage_program_scores.csv"
  "$step5d/plots/03_primary_diet_balanced_similarity.pdf"
  "$step5e/STEP_05E_COMPLETE.txt"
  "$step5e/tables/primary_r0.3_cluster_summary.csv"
  "$step5e/tables/primary_r0.3_parent_composition.csv"
  "$step5e/tables/primary_r0.3_targeted_RNA_promoter_summary_by_cluster.csv"
  "$step5e/tables/core_WNN_top_markers_review_resolutions.csv"
  "$step5e/plots/09_core_marker_audit_resolution_0.3.pdf"
  "$step5e/plots/12_primary_r0.3_program_heatmap.pdf"
  "$step5f/STEP_05F_COMPLETE.txt"
  "$step5f/input_manifest.csv"
  "$step5f/refined_working_labels.csv"
  "$step5f/matrix_cell_alignment_audit.csv"
  "$step5f/refined_cell_assignments.csv.gz"
  "$step5f/focused_composition_by_diet.csv"
  "$step5f/qc_and_phase_by_working_unit_diet.csv"
  "$step5f/r0.3_graph_and_boundary_stability.csv"
  "$step5f/focused_reduction_depth_correlations.csv"
  "$step5f/fetal_repair_sentinel_feature_availability.csv"
  "$step5f/sentinel_detection_by_cluster_diet.csv"
  "$step5f/fetal_revival_sentinel_co_detection.csv"
  "$step5f/marker_program_audit_summary.csv"
  "$step5f/RNA_only_r0.3_positive_markers.csv.gz"
  "$step5f/WNN_r0.3_positive_markers_separately_by_diet.csv.gz"
  "$step5f/targeted_promoter_evidence_with_availability.csv"
  "$step5f/plots/01_fetal_repair_sentinel_detection.pdf"
  "$step8e/STEP_08E_COMPLETE.txt"
  "$step8e/input_manifest.csv"
  "$step8e/matrix_cell_alignment_audit.csv"
  "$step8e/refined_diet_contrast_definition_and_counts.csv"
  "$step8e/persistent_reference_gene_canonicalization.csv"
  "$step8e/RNA_top_descriptive_effects.csv"
  "$step8e/ATAC_top_descriptive_effects.csv"
  "$step8e/persistent_genes_refined_WT_concordance.csv"
  "$step8e/persistent_gene_concordance_counts.csv"
  "$step8e/equal_cell_and_depth_matched_resampling.csv.gz"
  "$step8e/phase_stratum_eligibility.csv"
  "$step8e/RNA_equal_phase_weight_effect_sensitivity.csv"
  "$step8e/target_gene_promoter_overlapping_peak_map.csv.gz"
  "$step8e/target_gene_promoter_assay_availability.csv"
  "$step8e/refined_C0_C2_RNA_promoter_peak_evidence.csv.gz"
  "$step8e/plots/01_source_d21_vs_refined_WT_d21_RNA.pdf"
  "$step8e/plots/02_refined_RNA_depth_sensitivity.pdf"
)

missing=0
for relative in "${required[@]}"; do
  if [[ ! -s "$wt_output/$relative" ]]; then
    printf 'Missing required review file: %s\n' "$wt_output/$relative" >&2
    missing=1
  fi
done
(( missing == 0 )) || exit 1

mkdir -p "$pack_dir"
stamp="$(date +%Y-%m-%d_%H%M%S)"
archive="$pack_dir/wt_crypt_source_to_diet_review_${stamp}.tar.gz"
[[ ! -e "$archive" ]] || { printf 'Archive already exists: %s\n' "$archive" >&2; exit 1; }
temporary="$(mktemp "$pack_dir/.wt_crypt_pack_XXXXXXXX")"
trap 'rm -f "$temporary"' EXIT
tar -C "$wt_output" -czf "$temporary" "${required[@]}"
tar -tzf "$temporary" >/dev/null
mv -- "$temporary" "$archive"
trap - EXIT
sha256sum "$archive"
printf 'Review pack: %s\n' "$archive"
printf 'Files included: %d\n' "${#required[@]}"
printf 'Full genome-wide RNA/ATAC effect tables remain in %s\n' "$wt_output/$step8e"
