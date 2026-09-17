#!/usr/bin/env bash
set -euo pipefail

# Package the completed WT Step 08B source-matched persistent-HFD review.
#
# Usage from the repository root:
#   bash make_wt_step08b_review_pack.sh "$wt_output"
#
# Optional arguments:
#   1. WT cohort output directory
#   2. archive destination directory
#   3. project repository directory

wt_review_output="${1:-/scratch/dsaiz/Yesenia_scData2026/results/v2_resequenced/independent/wt}"
wt_review_destination="${2:-${wt_review_output}/review_packs}"
wt_review_project="${3:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"

step6_dir="${wt_review_output}/06_annotation"
step7_dir="${wt_review_output}/07_diet_response"
step8_dir="${wt_review_output}/08_persistent_hfd_programs"
step8b_dir="${wt_review_output}/08b_source_matched_persistent_hfd"

source_notebook="${wt_review_project}/analysis/08b_refine_persistent_hfd_programs_in_matched_WT_populations.Rmd"
source_guide="${wt_review_project}/analysis/README_STEP08_PERSISTENT_HFD.md"
reference_program="${wt_review_project}/config/datasets/v2_resequenced/reference_programs/persistent_hfd_genes.csv"

if [[ ! -d "$step8b_dir" ]]; then
  echo "Missing Step 08B output directory: $step8b_dir" >&2
  exit 1
fi

if [[ ! -f "$source_notebook" ]]; then
  echo "Missing Step 08B source notebook: $source_notebook" >&2
  exit 1
fi

if [[ ! -f "$reference_program" ]]; then
  echo "Missing persistent-HFD reference program: $reference_program" >&2
  exit 1
fi

mkdir -p "$wt_review_destination"

timestamp="$(date +%Y-%m-%d_%H%M%S)"
pack_name="wt_step08b_source_matched_review_${timestamp}"
archive_path="${wt_review_destination}/${pack_name}.tar.gz"
staging_root="$(mktemp -d "${TMPDIR:-/tmp}/wt_step08b_review.XXXXXX")"
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
    echo "Finish or rerun Step 08B before constructing the review pack." >&2
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
    echo "Optional context file not present; skipping: $source_path" >&2
    return 0
  fi

  mkdir -p "$(dirname "$destination_path")"
  cp -- "$source_path" "$destination_path"
}

step8b_required_files=(
  "STEP_08B_COMPLETE.txt"
  "step8b_parameters.csv"
  "matrix_cell_alignment_audit.csv"
  "step8b_run_inventory.csv"
  "step8b_unit_annotation_validation.csv"
  "matched_population_definitions.csv"
  "msigdb_input_audit.csv"
  "step8b_primary_population_summary.csv"
  "tables/matched_population_counts_by_diet.csv"
  "tables/matched_population_support.csv"
  "tables/RNA_pooled_library_effects_by_matched_population.csv.gz"
  "tables/RNA_top_diet_candidates_by_matched_population.csv"
  "tables/RNA_balanced_downsampling_stability.csv"
  "tables/persistent_gene_RNA_evidence_by_matched_population.csv.gz"
  "tables/strict_Step8_unit8_reproduction_audit.csv"
  "tables/persistent_program_ranked_RNA_enrichment.csv"
  "tables/persistent_gene_population_sensitivity.csv.gz"
  "tables/population_sensitivity_summary.csv"
  "tables/persistent_gene_primary_RNA_candidate_inventory.csv"
  "tables/persistent_gene_concordant_RNA_candidates.csv"
  "tables/persistent_gene_discordant_RNA_candidates.csv"
  "tables/targeted_ORA_input_audit.csv"
  "tables/targeted_MSigDB_ORA_results.csv"
  "tables/ATAC_feature_context_GRCm39.csv.gz"
  "tables/ATAC_pooled_library_effects_by_matched_population.csv.gz"
  "tables/ATAC_top_diet_candidates_by_matched_population.csv"
  "tables/ATAC_balanced_downsampling_stability.csv"
  "tables/ATAC_effect_summary_by_population_direction_context.csv"
  "tables/persistent_gene_exact_promoter_evidence_by_population_peak.csv.gz"
  "tables/persistent_gene_exact_promoter_summary.csv"
  "tables/persistent_gene_integrated_primary_candidates.csv"
  "plots/01_source_d21_vs_WT_d21_broad_stem_RNA.pdf"
  "plots/S01_source_d50_vs_WT_d21_broad_stem_RNA.pdf"
  "plots/02_persistent_program_ranked_RNA_enrichment.pdf"
  "plots/03_persistent_gene_population_sensitivity.pdf"
  "plots/04_primary_population_ATAC_effect_inventory.pdf"
  "plots/05_persistent_gene_RNA_exact_promoter_ATAC.pdf"
)

for relative_path in "${step8b_required_files[@]}"; do
  copy_required "$step8b_dir" "$relative_path" "08b_source_matched_persistent_hfd"
done

# Context needed to interpret analysis units and trace Step 08B back to the
# strict unit-8 audit. The large multiome RDS is intentionally excluded.
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

copy_optional \
  "$step7_dir" \
  "tables/analysis_unit_QC_by_diet.csv" \
  "context/07_diet_response"

copy_required \
  "$step8_dir" \
  "STEP_08_COMPLETE.txt" \
  "context/08_strict_unit8"

copy_required \
  "$step8_dir" \
  "tables/persistent_gene_primary_WT_stem_RNA_summary.csv" \
  "context/08_strict_unit8"

copy_required \
  "$step8_dir" \
  "tables/persistent_gene_exact_promoter_peak_overlaps.csv" \
  "context/08_strict_unit8"

mkdir -p "$pack_root/source/analysis"
mkdir -p "$pack_root/source/config"
cp -- "$source_notebook" "$pack_root/source/analysis/"
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
# WT Step 08B source-matched persistent-HFD review pack

Created: $(date --iso-8601=seconds)

This pack supports joint review of:

1. the strict unit-8 provenance check;
2. broad-stem versus strict-ISC and crypt-progenitor sensitivity;
3. all-goblet versus core-goblet and individual goblet components;
4. persistent-gene RNA concordance and cell-sampling stability;
5. independent ATAC diet effects in the same matched populations;
6. exact promoter RNA-ATAC agreement, opposition, or absence.

## Primary review questions

- Does broad stem (units 8 + 6) recover the source-derived persistent program
  more consistently than strict ISC-like unit 8 alone?
- Is that gain broadly distributed or driven almost entirely by crypt
  progenitor unit 6?
- Does immature goblet unit 13 weaken goblet concordance relative to the core
  goblet composite, and which genes account for that disagreement?
- Which broad-stem persistent genes are source-concordant, clear, and stable?
- For those genes, which exact promoter peaks show a stable diet effect in the
  matching direction, the opposing direction, or no evaluable effect?
- Are independent ATAC effects concentrated in promoters, gene bodies, or
  intergenic regions before any RNA-based filtering?

## Interpretation limits

- CON and HFD each represent one pooled WT library, so diet and library remain
  confounded and all effect estimates are descriptive.
- Cell downsampling evaluates sensitivity to cell composition and abundance;
  it does not create biological replication.
- Population-sensitivity comparisons reuse cells and are not independent
  statistical comparisons.
- Promoter overlap and nearest-gene context do not establish regulatory
  causality.
- Motif enrichment and cross-genotype comparisons are intentionally deferred
  until the broad-stem RNA and ATAC evidence is reviewed.

The multi-gigabyte Seurat object is not included. Source output directory:

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

tar -tzf "$archive_path" >/dev/null

echo "Created WT Step 08B review pack:"
echo "  $archive_path"
echo
echo "Archive size:"
du -h "$archive_path"
echo
echo "Packaged files:"
tar -tzf "$archive_path" | sort
