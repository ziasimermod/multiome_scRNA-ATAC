#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 3 ]]; then
  echo "Usage: $0 <wt_output_dir> [repo_root] [destination_dir]" >&2
  exit 2
fi

wt_output="$(cd "$1" && pwd)"
repo_root="${2:-$(pwd)}"
destination_dir="${3:-$repo_root}"

step5e_dir="$wt_output/05_wnn/focused_wt_crypt_stem_progenitor"
step5d_dir="$wt_output/05_wnn/gourab_d21_rna_reference_mapping/05d_symmetric_pseudobulk_validation"
step6_dir="$wt_output/06_annotation"

required_files=(
  "$step5e_dir/STEP_05E_COMPLETE.txt"
  "$step5e_dir/STEP_05E_READY_FOR_REVIEW.txt"
  "$step5e_dir/step5e_summary.csv"
  "$step5e_dir/tables/parent_cluster_annotation_audit.csv"
  "$step5e_dir/tables/focused_resolution_cluster_summary.csv"
  "$step5e_dir/tables/focused_resolution_parent_composition.csv"
  "$step5e_dir/tables/adjacent_resolution_crosswalk.csv"
  "$step5e_dir/tables/core_WNN_vs_RNA_crosswalk.csv"
  "$step5e_dir/tables/core_vs_boundary_WNN_crosswalk.csv"
  "$step5e_dir/tables/crypt_program_definitions.csv"
  "$step5e_dir/tables/crypt_program_gene_availability.csv"
  "$step5e_dir/tables/core_WNN_top_markers_review_resolutions.csv"
  "$step5e_dir/tables/core_WNN_marker_count_review_resolutions.csv"
  "$step5e_dir/tables/primary_r0.3_cluster_summary.csv"
  "$step5e_dir/tables/primary_r0.3_parent_composition.csv"
  "$step5e_dir/tables/primary_r0.3_targeted_RNA_promoter_summary_by_cluster.csv"
  "$step5e_dir/tables/primary_r0.3_targeted_RNA_promoter_summary_by_cluster_diet.csv"
  "$step5e_dir/tables/primary_r0.3_targeted_promoter_peak_map.csv"
  "$step5e_dir/tables/core_cell_assignments.csv.gz"
  "$step5e_dir/tables/boundary_cell_assignments.csv.gz"
  "$step5e_dir/plots/01_core_WNN_parent_clusters.pdf"
  "$step5e_dir/plots/02_core_WNN_by_diet.pdf"
  "$step5e_dir/plots/03_core_WNN_modality_weights.pdf"
  "$step5e_dir/plots/06_core_curated_program_support.pdf"
  "$step5e_dir/plots/07_core_source_program_support.pdf"
  "$step5e_dir/plots/08_core_cell_cycle_phase.pdf"
  "$step5e_dir/plots/09_core_marker_audit_resolution_0.3.pdf"
  "$step5e_dir/plots/10_primary_r0.3_targeted_RNA_evidence.pdf"
  "$step5e_dir/plots/11_primary_r0.3_targeted_promoter_ATAC_evidence.pdf"
  "$step5e_dir/plots/12_primary_r0.3_program_heatmap.pdf"
)

for required_file in "${required_files[@]}"; do
  if [[ ! -s "$required_file" ]]; then
    echo "Missing required review file: $required_file" >&2
    exit 1
  fi
done

mkdir -p "$destination_dir"

timestamp="$(date +%Y-%m-%d_%H%M%S)"
pack_name="wt_step05e_crypt_refinement_review_${timestamp}"
staging_root="$(mktemp -d)"
pack_dir="$staging_root/$pack_name"

cleanup() {
  rm -rf "$staging_root"
}
trap cleanup EXIT

mkdir -p \
  "$pack_dir/step05e/tables" \
  "$pack_dir/step05e/plots" \
  "$pack_dir/step05d_context/tables" \
  "$pack_dir/step05d_context/plots" \
  "$pack_dir/step06_context" \
  "$pack_dir/notebooks"

cp "$step5e_dir/STEP_05E_COMPLETE.txt" "$pack_dir/step05e/"
cp "$step5e_dir/STEP_05E_READY_FOR_REVIEW.txt" "$pack_dir/step05e/"
cp "$step5e_dir/step5e_summary.csv" "$pack_dir/step05e/"

find "$step5e_dir/tables" -maxdepth 1 -type f \
  \( -name '*.csv' -o -name '*.csv.gz' \) \
  -exec cp {} "$pack_dir/step05e/tables/" \;

find "$step5e_dir/plots" -maxdepth 1 -type f -name '*.pdf' \
  -exec cp {} "$pack_dir/step05e/plots/" \;

step5d_context_tables=(
  "source_integrated_cluster_top_markers.csv"
  "primary_similarity_with_bootstrap_support.csv"
  "primary_mutual_top1_matches.csv"
  "focused_WT_crypt_secretory_primary_matches.csv"
  "focused_WT_crypt_secretory_panel_and_diet_sensitivity.csv"
  "curated_lineage_program_scores.csv"
  "epithelial_cluster_cell_counts_by_diet.csv"
)

for table_name in "${step5d_context_tables[@]}"; do
  source_path="$step5d_dir/tables/$table_name"
  if [[ -s "$source_path" ]]; then
    cp "$source_path" "$pack_dir/step05d_context/tables/"
  fi
done

step5d_context_plots=(
  "03_primary_diet_balanced_similarity.pdf"
  "04_feature_panel_sensitivity.pdf"
  "05_diet_stratum_sensitivity.pdf"
  "06_curated_lineage_program_specificity.pdf"
  "07_focused_match_bootstrap_stability.pdf"
)

for plot_name in "${step5d_context_plots[@]}"; do
  source_path="$step5d_dir/plots/$plot_name"
  if [[ -s "$source_path" ]]; then
    cp "$source_path" "$pack_dir/step05d_context/plots/"
  fi
done

step6_context_files=(
  "$step6_dir/cell_state_annotations_applied.csv"
  "$step6_dir/cell_state_annotation_summary.csv"
  "$step6_dir/markers/RNA_top_annotation_candidates.csv"
  "$step6_dir/markers/RNA_all_positive_markers.csv"
  "$step6_dir/STEP_06_COMPLETE.txt"
)

for source_path in "${step6_context_files[@]}"; do
  if [[ -s "$source_path" ]]; then
    cp "$source_path" "$pack_dir/step06_context/"
  fi
done

step5e_notebook="$repo_root/analysis/05e_refine_wt_crypt_stem_progenitor_states.Rmd"
step5d_notebook="$repo_root/analysis/05d_validate_gourab_wt_annotations_with_symmetric_pseudobulk.Rmd"

for notebook_path in "$step5e_notebook" "$step5d_notebook"; do
  if [[ -s "$notebook_path" ]]; then
    cp "$notebook_path" "$pack_dir/notebooks/"
  fi
done

printf '%s\n' \
  'WT Step 05E focused crypt stem/progenitor review pack' \
  '' \
  'Review in this order:' \
  '1. step05e/STEP_05E_READY_FOR_REVIEW.txt' \
  '2. step05e/step5e_summary.csv' \
  '3. step05e/plots/01_core_WNN_parent_clusters.pdf' \
  '4. step05e/plots/04_core_WNN_vs_RNA_resolution_*.pdf' \
  '5. step05e/tables/adjacent_resolution_crosswalk.csv' \
  '6. step05e/tables/core_WNN_vs_RNA_crosswalk.csv' \
  '7. step05e/plots/05_boundary_WNN_resolution_*.pdf' \
  '8. step05e/tables/core_vs_boundary_WNN_crosswalk.csv' \
  '9. step05e/plots/06_core_curated_program_support.pdf' \
  '10. step05e/plots/07_core_source_program_support.pdf' \
  '11. step05e/plots/08_core_cell_cycle_phase.pdf' \
  '12. step05e/plots/09_core_marker_audit_resolution_0.3.pdf' \
  '13. step05e/tables/core_WNN_top_markers_review_resolutions.csv' \
  '14. step05e/tables/core_WNN_marker_count_review_resolutions.csv' \
  '15. step05e/plots/12_primary_r0.3_program_heatmap.pdf' \
  '16. step05e/plots/10_primary_r0.3_targeted_RNA_evidence.pdf' \
  '17. step05e/plots/11_primary_r0.3_targeted_promoter_ATAC_evidence.pdf' \
  '18. step05e/tables/primary_r0.3_targeted_RNA_promoter_summary_by_cluster_diet.csv' \
  '19. step05e/tables/primary_r0.3_targeted_promoter_peak_map.csv' \
  '20. step05e/tables/focused_resolution_cluster_summary.csv' \
  '21. step05d_context/ for the symmetric source comparison' \
  '' \
  'Interpretation guardrails:' \
  '- Step 05E is review-only and does not overwrite Step 6 annotations.' \
  '- A UMAP island or high relative module score is not sufficient evidence.' \
  '- Require adjacent-resolution stability, RNA sensitivity, boundary' \
  '  robustness, coherent detected markers, acceptable QC, and explicit' \
  '  cell-cycle and diet composition.' \
  '- RNA and exact-promoter ATAC evidence are evaluated separately; promoter' \
  '  accessibility supports but does not prove expression or identity.' \
  '- Diet comparisons remain descriptive because diet and pooled library are' \
  '  confounded.' \
  '- Large Seurat RDS files are deliberately excluded from this pack.' \
  > "$pack_dir/README.txt"

(
  cd "$pack_dir"
  find . -type f ! -name 'SHA256SUMS.txt' -print0 \
    | sort -z \
    | xargs -0 sha256sum > SHA256SUMS.txt
)

archive_path="$destination_dir/${pack_name}.tar.gz"
tar -C "$staging_root" -czf "$archive_path" "$pack_name"

echo "Created review pack: $archive_path"
echo "Contents:"
tar -tzf "$archive_path"
