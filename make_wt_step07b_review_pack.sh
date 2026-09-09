#!/usr/bin/env bash
set -euo pipefail

# Package the WT Step 7B manual-review outputs generated through Section 7.
#
# Usage from the repository root:
#   bash make_wt_step07b_review_pack.sh "$wt_output"
#
# Optional arguments:
#   1. WT output directory
#   2. archive destination directory
#   3. project repository directory

wt_review_output="${1:-/scratch/dsaiz/Yesenia_scData2026/results/v2_resequenced/independent/wt}"
wt_review_destination="${2:-${wt_review_output}/review_packs}"
wt_review_project="${3:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"

trajectory_dir="${wt_review_output}/07b_trajectories"
step6_dir="${wt_review_output}/06_annotation"
step7_dir="${wt_review_output}/07_diet_response"
source_notebook="${wt_review_project}/analysis/07b_model_wt_epithelial_trajectories.Rmd"

if [[ ! -d "$trajectory_dir" ]]; then
  echo "Missing Step 7B output directory: $trajectory_dir" >&2
  exit 1
fi

if [[ ! -f "$source_notebook" ]]; then
  echo "Missing Step 7B source notebook: $source_notebook" >&2
  exit 1
fi

mkdir -p "$wt_review_destination"

timestamp="$(date +%Y-%m-%d_%H%M%S)"
pack_name="wt_step07b_trajectory_review_${timestamp}"
archive_path="${wt_review_destination}/${pack_name}.tar.gz"
staging_root="$(mktemp -d "${TMPDIR:-/tmp}/wt_step07b_review.XXXXXX")"
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
    echo "Optional file not present; skipping: $source_path" >&2
    return 0
  fi

  mkdir -p "$(dirname "$destination_path")"
  cp -- "$source_path" "$destination_path"
}

trajectory_required_files=(
  "trajectory_decision_proposed.csv"
  "trajectory_analysis_unit_resolution.csv"
  "trajectory_cell_cycle_scoring_audit.csv"
  "trajectory_matrix_cell_alignment_audit.csv"
  "tables/candidate_unit_nuclei_by_library.csv"
  "tables/trajectory_program_feature_audit.csv"
  "tables/trajectory_program_summary_by_unit_and_library.csv"
  "tables/trajectory_program_scores_by_nucleus.csv.gz"
  "tables/candidate_unit_cell_cycle_by_library.csv"
  "tables/proposed_root_marker_expression.csv"
  "plots/01_absorptive_secretory_program_mixing.pdf"
  "plots/02_candidate_RNA_topology_unit_diet_phase.pdf"
  "plots/03_proposed_root_marker_audit.pdf"
)

for relative_path in "${trajectory_required_files[@]}"; do
  copy_required "$trajectory_dir" "$relative_path" "07b_trajectories"
done

copy_optional \
  "$trajectory_dir" \
  "STEP_07B_READY_FOR_TRAJECTORY_REVIEW.txt" \
  "07b_trajectories"

# Include the annotation and eligibility context needed to interpret the
# candidate units without repackaging the large Seurat object.
copy_required \
  "$step6_dir" \
  "cell_state_annotations_applied.csv" \
  "context/06_annotation"

copy_required \
  "$step7_dir" \
  "step7_summary.csv" \
  "context/07_diet_response"

copy_required \
  "$step7_dir" \
  "tables/analysis_population_eligibility.csv" \
  "context/07_diet_response"

copy_required \
  "$step7_dir" \
  "tables/analysis_unit_QC_by_diet.csv" \
  "context/07_diet_response"

mkdir -p "$pack_root/source"
cp -- "$source_notebook" "$pack_root/source/"

cat > "$pack_root/REVIEW_NOTES.txt" <<EOF
WT Step 7B trajectory review pack
Created: $(date --iso-8601=seconds)
Dataset: v2_resequenced
Cohort: wt
Included scope: outputs generated through Section 7
Post-approval Slingshot modeling included: FALSE
Source output directory: $wt_review_output
Source repository: $wt_review_project
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

tar -tzf "$archive_path" >/dev/null

echo "Created review pack:"
echo "  $archive_path"
echo
echo "Archive size:"
du -h "$archive_path"
echo
echo "Packaged files:"
tar -tzf "$archive_path" | sort

