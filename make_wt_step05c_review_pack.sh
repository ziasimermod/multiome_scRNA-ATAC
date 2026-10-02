#!/usr/bin/env bash
set -euo pipefail

# Review-pack builder version 2.0.1.
# Cache-safe release identifier: step05c-review-pack-v2.0.1-20260921
# Package the completed bidirectional Gourab d21 <-> WT annotation audit,
# together with the WT marker, annotation, and WNN context needed to diagnose
# asymmetric or potentially miscalibrated inverse mapping.
#
# Usage from the repository root:
#   bash make_wt_step05c_review_pack.sh "$wt_output"
#
# Optional arguments:
#   1. WT cohort output directory
#   2. archive destination directory
#   3. project repository directory

wt_review_output="${1:-/scratch/dsaiz/Yesenia_scData2026/results/v2_resequenced/independent/wt}"
wt_review_destination="${2:-${wt_review_output}/review_packs}"
wt_review_project="${3:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"

step5c_dir="${wt_review_output}/05_wnn/gourab_d21_rna_reference_mapping"
step5_dir="${wt_review_output}/05_wnn"
step6_dir="${wt_review_output}/06_annotation"
wt_config_dir="${wt_review_project}/config/datasets/v2_resequenced/cohorts/wt"
dataset_config_dir="${wt_review_project}/config/datasets/v2_resequenced"
source_notebook_v2="${wt_review_project}/analysis/05c_map_gourab_d21_scRNA_to_wt_multiome_v2.Rmd"
source_notebook_original="${wt_review_project}/analysis/05c_map_gourab_d21_scRNA_to_wt_multiome.Rmd"
package_setup_v2="${wt_review_project}/setup/install_packages_v2.R"
package_setup_original="${wt_review_project}/setup/install_packages.R"

if [[ -f "$source_notebook_v2" ]]; then
  source_notebook="$source_notebook_v2"
else
  source_notebook="$source_notebook_original"
fi

if [[ -f "$package_setup_v2" ]]; then
  package_setup="$package_setup_v2"
else
  package_setup="$package_setup_original"
fi

for required_path in \
  "$step5c_dir" \
  "$source_notebook" \
  "$package_setup"; do
  if [[ ! -e "$required_path" ]]; then
    echo "Missing required Step 05C review input: $required_path" >&2
    exit 1
  fi
done

mkdir -p "$wt_review_destination"

timestamp="$(date +%Y-%m-%d_%H%M%S)"
pack_name="wt_step05c_bidirectional_annotation_review_v2_${timestamp}"
archive_path="${wt_review_destination}/${pack_name}.tar.gz"
staging_root="$(mktemp -d "${TMPDIR:-/tmp}/wt_step05c_review.XXXXXX")"
pack_root="${staging_root}/${pack_name}"

cleanup() {
  rm -rf -- "$staging_root"
}
trap cleanup EXIT

mkdir -p "$pack_root"

copy_required() {
  local source_root="$1"
  local relative_path="$2"
  local destination_prefix="$3"
  local source_path="${source_root}/${relative_path}"
  local destination_path="${pack_root}/${destination_prefix}/${relative_path}"

  if [[ ! -f "$source_path" ]]; then
    echo "Missing required review file: $source_path" >&2
    echo "Finish or rerun the updated Step 05C notebook first." >&2
    exit 1
  fi

  mkdir -p "$(dirname "$destination_path")"
  cp -- "$source_path" "$destination_path"
}

copy_optional() {
  local source_root="$1"
  local relative_path="$2"
  local destination_prefix="$3"
  local source_path="${source_root}/${relative_path}"
  local destination_path="${pack_root}/${destination_prefix}/${relative_path}"

  if [[ ! -f "$source_path" ]]; then
    echo "Optional review file not present; skipping: $source_path" >&2
    return 0
  fi

  mkdir -p "$(dirname "$destination_path")"
  cp -- "$source_path" "$destination_path"
}

copy_optional_by_basename() {
  local source_root="$1"
  local requested_basename="$2"
  local destination_prefix="$3"
  local source_path
  local relative_path
  local destination_path
  local match_count=0

  if [[ ! -d "$source_root" ]]; then
    echo "Optional context directory not present; skipping: $source_root" >&2
    return 0
  fi

  while IFS= read -r -d '' source_path; do
    relative_path="${source_path#${source_root}/}"
    destination_path="${pack_root}/${destination_prefix}/${relative_path}"
    mkdir -p "$(dirname "$destination_path")"
    cp -- "$source_path" "$destination_path"
    match_count=$((match_count + 1))
  done < <(
    find "$source_root" \
      -type f \
      -name "$requested_basename" \
      -print0 | sort -z
  )

  if [[ "$match_count" -eq 0 ]]; then
    echo "Optional context file not found recursively: ${source_root}/**/${requested_basename}" >&2
  fi
}

copy_reviewable_tree() {
  local source_root="$1"
  local destination_prefix="$2"
  local source_path
  local relative_path
  local destination_path

  while IFS= read -r -d '' source_path; do
    relative_path="${source_path#${source_root}/}"
    destination_path="${pack_root}/${destination_prefix}/${relative_path}"
    mkdir -p "$(dirname "$destination_path")"
    cp -- "$source_path" "$destination_path"
  done < <(
    find "$source_root" \
      -type f \
      \( \
        -name '*.csv' -o \
        -name '*.csv.gz' -o \
        -name '*.tsv' -o \
        -name '*.tsv.gz' -o \
        -name '*.txt' -o \
        -name '*.pdf' -o \
        -name '*.log' \
      \) \
      -print0 | sort -z
  )
}

step5c_required_files=(
  "STEP_05C_COMPLETE.txt"
  "STEP_05C_READY_FOR_REVIEW.txt"
  "step5c_summary.csv"
  "WT_Gourab_d21_label_transfer_assignments.csv"
  "Gourab_d21_WT_inverse_label_transfer_assignments.csv"
  "tables/source_integration_suffix_audit.csv"
  "tables/source_label_counts_by_diet.csv"
  "tables/source_cell_recovery_audit.csv"
  "tables/source_count_metadata_concordance.csv"
  "tables/source_marker_panel_availability.csv"
  "tables/source_RNA_all_positive_markers.csv"
  "tables/source_RNA_top_annotation_candidates.csv"
  "tables/reference_query_feature_overlap_audit.csv"
  "tables/WT_label_transfer_score_summary.csv"
  "tables/WT_WNN_cluster_by_Gourab_d21_state_crosswalk.csv"
  "tables/WT_WNN_cluster_mapping_confidence.csv"
  "tables/WT_current_cell_state_by_Gourab_d21_state_crosswalk.csv"
  "tables/WT_focus_clusters_0_6_8_and_secretory_crosswalk.csv"
  "tables/WT_inverse_reference_composition.csv"
  "tables/inverse_mapping_score_summary.csv"
  "tables/Gourab_state_by_predicted_WT_state_crosswalk.csv"
  "tables/Gourab_state_by_predicted_WT_cluster_crosswalk.csv"
  "tables/inverse_source_label_specificity_summary.csv"
  "tables/reciprocal_WT_state_Gourab_state_agreement.csv"
  "tables/reciprocal_WT_cluster_Gourab_state_agreement.csv"
  "tables/reciprocal_focus_cluster_Gourab_state_agreement.csv"
  "tables/Harmony_sensitivity_audit.csv"
  "plots/01_source_reference_UMAP.pdf"
  "plots/02_source_marker_DotPlot.pdf"
  "plots/03_WT_WNN_reference_mapping_overview.pdf"
  "plots/04_focused_cluster_mapping_crosswalk.pdf"
  "plots/05_inverse_Gourab_to_WT_mapping_overview.pdf"
  "plots/06_inverse_Gourab_to_WT_state_crosswalk.pdf"
  "plots/07_reciprocal_focused_cluster_agreement.pdf"
)

for relative_path in "${step5c_required_files[@]}"; do
  copy_required "$step5c_dir" "$relative_path" "05c_bidirectional_mapping"
done

# Copy any additional reviewable outputs introduced by a later Step 05C
# revision. Large RDS objects and raw matrices are deliberately excluded.
copy_reviewable_tree "$step5c_dir" "05c_bidirectional_mapping"

# Preserve the applied WT annotations used to build the inverse reference.
copy_required \
  "$step6_dir" \
  "cell_state_annotations_applied.csv" \
  "context/06_annotation"

step6_context_files=(
  "STEP_06_COMPLETE.txt"
  "cell_state_annotation_summary.csv"
  "RNA_top_annotation_candidates.csv"
  "RNA_all_positive_markers.csv"
  "14_WNN_UMAP_broad_compartment.pdf"
  "15_WNN_UMAP_cell_state.pdf"
  "16_cell_state_sample_composition_descriptive.pdf"
)

for relative_path in "${step6_context_files[@]}"; do
  copy_optional_by_basename \
    "$step6_dir" \
    "$relative_path" \
    "context/06_annotation"
done

step5_context_files=(
  "chosen_cluster_sizes_resolution_0.8.csv"
  "cluster_sample_composition_descriptive.csv"
  "modality_weight_summary.csv"
  "02_WNN_UMAP_by_sample.pdf"
  "05_colon_marker_DotPlot_resolution_0.8.pdf"
  "06_chosen_WNN_clusters_resolution_0.8.pdf"
)

for relative_path in "${step5_context_files[@]}"; do
  copy_optional_by_basename \
    "$step5_dir" \
    "$relative_path" \
    "context/05_wnn"
done

wt_config_files=(
  "clustering_decision.csv"
  "cell_state_annotations.csv"
)

for relative_path in "${wt_config_files[@]}"; do
  copy_optional "$wt_config_dir" "$relative_path" "source/config/wt"
done

copy_optional \
  "$dataset_config_dir" \
  "samples.csv" \
  "source/config/dataset"

mkdir -p "$pack_root/source/analysis"
mkdir -p "$pack_root/source/setup"
cp -- "$source_notebook" "$pack_root/source/analysis/"
cp -- "$package_setup" "$pack_root/source/setup/"

mkdir -p "$pack_root/provenance"

{
  printf 'relative_path\tbytes\n'
  find "$step5c_dir" \
    -type f \
    -printf '%P\t%s\n' | sort
} > "$pack_root/provenance/step05c_complete_output_inventory.tsv"

{
  printf 'relative_path\tbytes\treason_excluded\n'
  find "$step5c_dir" \
    -type f \
    \( -name '*.rds' -o -name '*.h5' -o -name '*.h5ad' \) \
    -printf '%P\t%s\tlarge_binary_object\n' | sort
} > "$pack_root/provenance/excluded_large_objects.tsv"

{
  printf 'Created: %s\n' "$(date --iso-8601=seconds)"
  printf 'Repository: %s\n' "$wt_review_project"
  printf 'Branch: %s\n' "$(git -C "$wt_review_project" branch --show-current 2>/dev/null || printf 'unavailable')"
  printf 'Commit: %s\n' "$(git -C "$wt_review_project" rev-parse HEAD 2>/dev/null || printf 'unavailable')"
  printf '\nGit status at packaging:\n'
  git -C "$wt_review_project" status --short --branch 2>/dev/null || true
} > "$pack_root/provenance/git_provenance.txt"

cat > "$pack_root/README_REVIEW_PACK.md" <<EOF
# WT Step 05C bidirectional annotation review pack v2

Created: $(date --iso-8601=seconds)

This pack tests annotation correspondence in both directions:

1. Gourab d21 labels transferred onto reviewed WT snMultiome RNA;
2. reviewed WT cell states and WNN clusters transferred onto Gourab d21 RNA;
3. pairwise forward/reverse percentages and conservative reciprocal support;
4. source-label specificity, fragmentation, and per-cell mapping uncertainty;
5. independent marker evidence recalculated from the full source H5 counts.
6. WT de novo marker and annotation evidence used to construct the inverse
   reference, when present in the Step 6 output directory.

## Review order

1. Confirm STEP_05C_COMPLETE.txt and step5c_summary.csv.
2. Inspect source cell recovery and count-metadata concordance.
3. Review source_RNA_top_annotation_candidates.csv and the source marker
   DotPlot before trusting the inherited labels.
4. Inspect both mapping-overview PDFs on their original embeddings.
5. Review inverse_source_label_specificity_summary.csv. A broad source label
   that maps to several reviewed WT states is not automatically wrong, but it
   is less specific than its name may imply.
6. Review reciprocal_focus_cluster_Gourab_state_agreement.csv and its heatmap.
   Strong one-way mapping with weak reverse mapping indicates containment or
   granularity mismatch, not one-to-one equivalence.
7. Use the full assignment tables to inspect low-score, high-entropy, or
   low-margin cells behind any surprising aggregate result.
8. Compare the inverse assignments against WT_inverse_reference_composition.csv
   and the WT Step 6 marker tables. Unequal reference-state sizes, different
   label granularity, and cycling programs can all distort maximum transfer
   scores without implying that either annotation is simply correct or wrong.

## Questions this pack is designed to resolve

- Does the inverse mapping fail because individual source cells have weak or
  ambiguous WT support, or because a coherent Gourab label spans several more
  finely resolved WT states?
- Are maximum inverse scores inflated toward large WT reference classes?
- Do cell-state and raw WNN-cluster transfers tell the same story?
- Does the source Stem label contain ISC-like, crypt-progenitor, and cycling
  states in proportions consistent with an intentionally broad label?
- Do DSC and Goblet labels map onto biologically adjacent WT secretory states,
  or across unrelated compartments?
- Are apparent disagreements concentrated by diet, RNA complexity, or source
  integrated cluster?

## Interpretation limits

- Transfer scores are algorithmic mapping support, not biological identity
  probabilities and not substitutes for annotation confidence.
- Reciprocal agreement exposes asymmetry but is not an independent validation
  experiment; the WT review used related intestinal marker knowledge.
- Whole-cell scRNA and snMultiome RNA differ in preparation and transcript
  capture, so incomplete concordance is expected.
- Neither source nor WT has independent mouse-level diet replicates here.
  No diet-effect inference or abundance test is supported.
- Harmony, when available, is supplementary visualization only.
- No WT or Gourab annotation is changed by this notebook.

Large Seurat objects and raw H5 matrices are intentionally excluded. Source
output directory:

    $step5c_dir
EOF

manifest_path="$pack_root/MANIFEST.tsv"

{
  printf 'relative_path\tbytes\tsha256\n'

  while IFS= read -r -d '' packaged_file; do
    relative_path="${packaged_file#${pack_root}/}"
    file_bytes="$(stat -c '%s' "$packaged_file")"
    file_sha256="$(sha256sum "$packaged_file" | awk '{print $1}')"
    printf '%s\t%s\t%s\n' \
      "$relative_path" \
      "$file_bytes" \
      "$file_sha256"
  done < <(
    find "$pack_root" \
      -type f \
      ! -name 'MANIFEST.tsv' \
      -print0 | sort -z
  )
} > "$manifest_path"

tar \
  -C "$staging_root" \
  -czf "$archive_path" \
  "$pack_name"

printf 'Created Step 05C review pack:\n%s\n' "$archive_path"
