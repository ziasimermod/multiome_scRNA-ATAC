#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 3 ]]; then
  echo "Usage: $0 <wt_output_dir> [repo_root] [destination_dir]" >&2
  exit 2
fi

wt_output="$(cd "$1" && pwd)"
repo_root="${2:-$(pwd)}"
destination_dir="${3:-$repo_root}"

step5c_dir="$wt_output/05_wnn/gourab_d21_rna_reference_mapping"
step5d_dir="$step5c_dir/05d_symmetric_pseudobulk_validation"

required_files=(
  "$step5d_dir/STEP_05D_COMPLETE.txt"
  "$step5d_dir/STEP_05D_READY_FOR_REVIEW.txt"
  "$step5d_dir/step5d_summary.csv"
  "$step5d_dir/tables/step5c_inverse_qc_disposition.csv"
  "$step5d_dir/tables/primary_similarity_with_bootstrap_support.csv"
  "$step5d_dir/tables/primary_mutual_top1_matches.csv"
  "$step5d_dir/tables/focused_WT_crypt_secretory_primary_matches.csv"
  "$step5d_dir/tables/focused_WT_crypt_secretory_panel_and_diet_sensitivity.csv"
  "$step5d_dir/plots/03_primary_diet_balanced_similarity.pdf"
  "$step5d_dir/plots/06_curated_lineage_program_specificity.pdf"
  "$step5d_dir/plots/07_focused_match_bootstrap_stability.pdf"
)

for required_file in "${required_files[@]}"; do
  if [[ ! -s "$required_file" ]]; then
    echo "Missing required review file: $required_file" >&2
    exit 1
  fi
done

mkdir -p "$destination_dir"

timestamp="$(date +%Y-%m-%d_%H%M%S)"
pack_name="wt_step05d_symmetric_pseudobulk_review_${timestamp}"
staging_root="$(mktemp -d)"
pack_dir="$staging_root/$pack_name"

cleanup() {
  rm -rf "$staging_root"
}
trap cleanup EXIT

mkdir -p \
  "$pack_dir/step05d/tables" \
  "$pack_dir/step05d/plots" \
  "$pack_dir/step05c_context/tables" \
  "$pack_dir/step05c_context/plots" \
  "$pack_dir/step06_context" \
  "$pack_dir/notebooks"

cp "$step5d_dir/STEP_05D_COMPLETE.txt" "$pack_dir/step05d/"
cp "$step5d_dir/STEP_05D_READY_FOR_REVIEW.txt" "$pack_dir/step05d/"
cp "$step5d_dir/step5d_summary.csv" "$pack_dir/step05d/"

find "$step5d_dir/tables" -maxdepth 1 -type f -name '*.csv' \
  -exec cp {} "$pack_dir/step05d/tables/" \;

find "$step5d_dir/plots" -maxdepth 1 -type f -name '*.pdf' \
  -exec cp {} "$pack_dir/step05d/plots/" \;

step5c_context_tables=(
  "source_label_counts_by_diet.csv"
  "source_RNA_top_annotation_candidates.csv"
  "WT_focus_clusters_0_6_8_and_secretory_crosswalk.csv"
  "Gourab_state_by_predicted_WT_cluster_crosswalk.csv"
  "Gourab_state_by_predicted_WT_state_crosswalk.csv"
  "inverse_source_label_specificity_summary.csv"
  "reciprocal_focus_cluster_Gourab_state_agreement.csv"
  "Harmony_sensitivity_audit.csv"
)

for table_name in "${step5c_context_tables[@]}"; do
  source_path="$step5c_dir/tables/$table_name"
  if [[ -s "$source_path" ]]; then
    cp "$source_path" "$pack_dir/step05c_context/tables/"
  fi
done

step5c_context_plots=(
  "03_WT_WNN_reference_mapping_overview.pdf"
  "05_inverse_Gourab_to_WT_mapping_overview.pdf"
  "06_inverse_Gourab_to_WT_state_crosswalk.pdf"
  "08_reciprocal_focus_cluster_agreement.pdf"
  "09_Harmony_platform_alignment_sensitivity.pdf"
)

for plot_name in "${step5c_context_plots[@]}"; do
  source_path="$step5c_dir/plots/$plot_name"
  if [[ -s "$source_path" ]]; then
    cp "$source_path" "$pack_dir/step05c_context/plots/"
  fi
done

step6_files=(
  "$wt_output/06_annotation/cell_state_annotations_applied.csv"
  "$wt_output/06_annotation/cell_state_annotation_summary.csv"
  "$wt_output/06_annotation/markers/RNA_top_annotation_candidates.csv"
  "$wt_output/06_annotation/STEP_06_COMPLETE.txt"
)

for source_path in "${step6_files[@]}"; do
  if [[ -s "$source_path" ]]; then
    cp "$source_path" "$pack_dir/step06_context/"
  fi
done

notebook_path="$repo_root/analysis/05d_validate_gourab_wt_annotations_with_symmetric_pseudobulk.Rmd"
if [[ -s "$notebook_path" ]]; then
  cp "$notebook_path" "$pack_dir/notebooks/"
fi

cat > "$pack_dir/README.txt" <<'EOF'
WT Step 05D symmetric pseudobulk annotation-validation review pack

Review in this order:
1. step05d/tables/step5c_inverse_qc_disposition.csv
2. step05d/plots/03_primary_diet_balanced_similarity.pdf
3. step05d/tables/primary_similarity_with_bootstrap_support.csv
4. step05d/plots/07_focused_match_bootstrap_stability.pdf
5. step05d/tables/focused_WT_crypt_secretory_primary_matches.csv
6. step05d/plots/06_curated_lineage_program_specificity.pdf
7. step05d/plots/04_feature_panel_sensitivity.pdf
8. step05d/plots/05_diet_stratum_sensitivity.pdf

The pseudobulk RDS and Seurat objects are deliberately excluded because they
are large and reproducible from the notebook. Step 05D does not overwrite any
WT annotation and does not perform biological diet inference.
EOF

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
