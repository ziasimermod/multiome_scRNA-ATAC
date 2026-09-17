#!/usr/bin/env bash
set -euo pipefail

# Package the completed WT Step 08C composition-decomposition review.
#
# Usage from the repository root:
#   bash make_wt_step08c_review_pack.sh "$wt_output"
#
# Optional arguments:
#   1. WT cohort output directory
#   2. archive destination directory
#   3. project repository directory

wt_review_output="${1:-/scratch/dsaiz/Yesenia_scData2026/results/v2_resequenced/independent/wt}"
wt_review_destination="${2:-${wt_review_output}/review_packs}"
wt_review_project="${3:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"

step8b_dir="${wt_review_output}/08b_source_matched_persistent_hfd"
step8c_dir="${wt_review_output}/08c_broad_stem_decomposition"

source_notebook="${wt_review_project}/analysis/08c_decompose_wt_broad_stem_and_prioritize_regulatory_loci.Rmd"
source_guide="${wt_review_project}/analysis/README_STEP08_PERSISTENT_HFD.md"
motif_decision="${wt_review_project}/config/datasets/v2_resequenced/cohorts/wt/motif_peak_set_decision.csv"
reference_program="${wt_review_project}/config/datasets/v2_resequenced/reference_programs/persistent_hfd_genes.csv"

for required_path in \
  "$step8c_dir" \
  "$source_notebook" \
  "$motif_decision" \
  "$reference_program"; do
  if [[ ! -e "$required_path" ]]; then
    echo "Missing required Step 08C review input: $required_path" >&2
    exit 1
  fi
done

mkdir -p "$wt_review_destination"

timestamp="$(date +%Y-%m-%d_%H%M%S)"
pack_name="wt_step08c_broad_stem_decomposition_review_${timestamp}"
archive_path="${wt_review_destination}/${pack_name}.tar.gz"
staging_root="$(mktemp -d "${TMPDIR:-/tmp}/wt_step08c_review.XXXXXX")"
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
    echo "Finish or rerun Step 08C before constructing the review pack." >&2
    exit 1
  fi

  mkdir -p "$(dirname "$destination_path")"
  cp -- "$source_path" "$destination_path"
}

step8c_required_files=(
  "STEP_08C_COMPLETE.txt"
  "step8c_parameters.csv"
  "matrix_cell_alignment_audit.csv"
  "step8c_unit_annotation_validation.csv"
  "step8c_classification_summary.csv"
  "step8c_summary.csv"
  "tables/broad_stem_component_counts_and_proportions.csv"
  "tables/RNA_linear_component_means.csv.gz"
  "tables/RNA_broad_stem_composition_decomposition.csv.gz"
  "tables/persistent_gene_RNA_composition_decomposition.csv"
  "tables/persistent_gene_RNA_mechanism_summary.csv"
  "tables/ATAC_detection_component_means.csv.gz"
  "tables/ATAC_binary_detection_downsampling_stability.csv.gz"
  "tables/ATAC_broad_stem_composition_decomposition.csv.gz"
  "tables/motif_ready_peak_sets.csv.gz"
  "tables/motif_peak_set_summary.csv"
  "tables/motif_background_universe.csv.gz"
  "tables/persistent_gene_exact_promoter_decomposition.csv.gz"
  "tables/persistent_gene_prioritized_promoter_candidates.csv"
  "plots/01_broad_stem_component_composition.pdf"
  "plots/02_persistent_RNA_observed_vs_standardized.pdf"
  "plots/03_persistent_RNA_unit8_vs_unit6.pdf"
  "plots/04_ATAC_unit8_vs_unit6_response_classes.pdf"
  "plots/05_motif_ready_peak_set_inventory.pdf"
)

for relative_path in "${step8c_required_files[@]}"; do
  copy_required "$step8c_dir" "$relative_path" "08c_broad_stem_decomposition"
done

# The promoter plot is absent only when no candidate passes the explicit
# prioritization rules; preserve that outcome without making packaging fail.
if [[ -f "$step8c_dir/plots/06_prioritized_persistent_gene_promoter_candidates.pdf" ]]; then
  copy_required \
    "$step8c_dir" \
    "plots/06_prioritized_persistent_gene_promoter_candidates.pdf" \
    "08c_broad_stem_decomposition"
fi

# Compact Step 08B context needed to compare observed pooled effects with the
# new standardized/component-resolved effects.
step8b_context_files=(
  "STEP_08B_COMPLETE.txt"
  "tables/matched_population_counts_by_diet.csv"
  "tables/persistent_gene_RNA_evidence_by_matched_population.csv.gz"
  "tables/ATAC_effect_summary_by_population_direction_context.csv"
  "tables/persistent_gene_exact_promoter_summary.csv"
  "tables/persistent_gene_integrated_primary_candidates.csv"
)

for relative_path in "${step8b_context_files[@]}"; do
  copy_required "$step8b_dir" "$relative_path" "context/08b_source_matched"
done

mkdir -p "$pack_root/source/analysis"
mkdir -p "$pack_root/source/config"
cp -- "$source_notebook" "$pack_root/source/analysis/"
cp -- "$motif_decision" "$pack_root/source/config/"
cp -- "$reference_program" "$pack_root/source/config/"

if [[ -f "$source_guide" ]]; then
  cp -- "$source_guide" "$pack_root/source/analysis/"
fi

mkdir -p "$pack_root/provenance"

{
  printf 'Created: %s\n' "$(date --iso-8601=seconds)"
  printf 'Repository: %s\n' "$wt_review_project"
  printf 'Branch: %s\n' "$(git -C "$wt_review_project" branch --show-current 2>/dev/null || printf 'unavailable')"
  printf 'Commit: %s\n' "$(git -C "$wt_review_project" rev-parse HEAD 2>/dev/null || printf 'unavailable')"
  printf '\nGit status at packaging:\n'
  git -C "$wt_review_project" status --short --branch 2>/dev/null || true
} > "$pack_root/provenance/git_provenance.txt"

cat > "$pack_root/README_REVIEW_PACK.md" <<EOF
# WT Step 08C broad-stem decomposition review pack

Created: $(date --iso-8601=seconds)

This pack supports review of four separate questions:

1. How much did the unit-8/unit-6 composition change between pooled diet
   libraries?
2. Which persistent-gene RNA effects remain after fixing component weights?
3. Which effects are shared, unit-6-specific, unit-8-specific, opposing, or
   composition dominant?
4. Which stable ATAC peak sets and exact-promoter candidates are appropriate
   to approve for motif enrichment?

## Required review order

1. Read STEP_08C_COMPLETE.txt and step8c_summary.csv.
2. Inspect the composition plot and component-count table.
3. Review persistent_gene_RNA_composition_decomposition.csv.
4. Review motif_peak_set_summary.csv and the ATAC component scatter.
5. Review persistent_gene_prioritized_promoter_candidates.csv.
6. Edit motif_peak_set_decision.csv only after those checks. Step 08D will
   refuse to run until at least one set has approved=TRUE with reviewer and
   review_date completed.

## Interpretation limits

- CON and HFD each represent one pooled WT library; all diet effects remain
  descriptive and library confounded.
- The exact midpoint decomposition separates observed mixture and within-state
  terms but does not create independent replicates.
- Cell downsampling evaluates sampling stability only.
- Peak opening is accessibility, not necessarily enhancer activation.
- Promoter overlap and motif enrichment do not prove regulatory causality.

The large Step 6 Seurat object is intentionally excluded. Source output:

    $wt_review_output
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

printf 'Created Step 08C review pack:\n%s\n' "$archive_path"
