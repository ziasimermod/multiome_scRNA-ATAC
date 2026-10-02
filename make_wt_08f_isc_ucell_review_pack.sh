#!/usr/bin/env bash
# Create a compact 08F review pack from completed WT SOL output.
set -euo pipefail

usage() {
  echo 'Usage: bash make_wt_08f_isc_ucell_review_pack.sh [WT_OUTPUT_DIR] [PACK_DIR]' >&2
  echo 'Defaults: WT_OUTPUT_DIR=/scratch/dsaiz/Yesenia_scData2026/results/v2_resequenced/independent/wt' >&2
  echo '          PACK_DIR=WT_OUTPUT_DIR/review_packs' >&2
}
if (( $# > 2 )); then usage; exit 2; fi
wt_root="${1:-/scratch/dsaiz/Yesenia_scData2026/results/v2_resequenced/independent/wt}"
pack_dir="${2:-$wt_root/review_packs}"
effect_dir="$wt_root/08f_isc_ucell_wt_crypt"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
summarizer="$script_dir/analysis/08f_summarize_isc_ucell_review.R"
notebook="$script_dir/analysis/08f_score_wt_crypt_isc_ucell_and_accessibility.Rmd"

[[ -d "$effect_dir" ]] || { echo "Missing Step 08F output: $effect_dir" >&2; exit 1; }
[[ -s "$effect_dir/STEP_08F_COMPLETE.txt" ]] || {
  echo "Step 08F is not complete: $effect_dir/STEP_08F_COMPLETE.txt" >&2; exit 1;
}
[[ -s "$summarizer" && -s "$notebook" ]] || {
  echo "Run this script from the checkout with the 08F notebook and summarizer." >&2; exit 1;
}
command -v Rscript >/dev/null || { echo 'Rscript is missing from PATH.' >&2; exit 1; }
command -v tar >/dev/null || { echo 'tar is missing from PATH.' >&2; exit 1; }
mkdir -p "$pack_dir"

stage="$(mktemp -d "$pack_dir/.wt08f_review.XXXXXXXX")"
trap 'rm -rf -- "$stage"' EXIT
content="$stage/wt_08f_isc_ucell_review"
mkdir -p "$content/tables" "$content/08f/plots" "$content/scripts"
Rscript "$summarizer" "$effect_dir" "$content/tables"

required=(
  STEP_08F_COMPLETE.txt
  input_manifest.csv
  sessionInfo.txt
  00_matrix_alignment.csv
  01_signature_assay_availability.csv
  01_marker_detection_by_parent_diet.csv
  01_scores_and_preserved_annotations.csv.gz
  02_all_gate_memberships_with_05f_labels.csv.gz
  02_gate_counts_and_depth_by_parent_diet.csv
  02_gate_boundary_and_ties.csv
  02_maxRank400_vs_600_gate_overlap.csv
  02_contrast_definitions_and_eligibility.csv
  03_selected_05f_subclusters_by_parent_diet.csv
  06_which_parent_effects_are_available.csv
  06_peak_comparison_counts.csv
  06_parent_identity_and_depth_review.csv
  07_ATAC_equal_cell_and_depth_sensitivity.csv.gz
  08_ATAC_candidates_with_posthoc_promoter_context.csv.gz
  08_Gourab_gene_assay_audit.csv
  08_posthoc_Gourab_RNA_effects_by_population.csv
  08_source_gene_promoter_peak_annotations.csv.gz
  09_prespecified_ATAC_locus_examples.csv
  09_ATAC_first_examples_with_RNA_context.csv
  09_coverage_plot_status.csv
)
for item in "${required[@]}"; do
  [[ -s "$effect_dir/$item" ]] || {
    echo "Missing/nonempty required output: $effect_dir/$item" >&2; exit 1;
  }
  cp -- "$effect_dir/$item" "$content/08f/$item"
done

# Source PDFs plus coverage, if plotted. Failed coverage attempts are recorded
# in 09_coverage_plot_status.csv and do not prevent a useful pack.
plots=(
  01_ISC_UCell_on_global_WNN.pdf
  02_nested_gate_locations.pdf
  03_parent_score_distributions.pdf
  04_preserved_subcluster_composition.pdf
  05_proxy_vs_union_peak_effects.pdf
)
for item in "${plots[@]}"; do
  [[ -s "$effect_dir/plots/$item" ]] || {
    echo "Missing required figure: $effect_dir/plots/$item" >&2; exit 1;
  }
  cp -- "$effect_dir/plots/$item" "$content/08f/plots/$item"
done
shopt -s nullglob
for file in "$effect_dir"/plots/08_locus_ATAC_and_RNA_context.pdf \
            "$effect_dir"/plots/09_*_locus_*.pdf \
            "$effect_dir"/09_*_review_only_GRCm39.bed; do
  if [[ -s "$file" ]]; then
    case "$file" in
      */plots/*) cp -- "$file" "$content/08f/plots/" ;;
      *) cp -- "$file" "$content/08f/" ;;
    esac
  fi
done
shopt -u nullglob

cp -- "$summarizer" "$notebook" "$content/scripts/"
stamp="$(date -u +%Y-%m-%d_%H%M%S_UTC)"
archive="$pack_dir/wt_08f_isc_ucell_review_${stamp}.tar.gz"
tar -C "$stage" -czf "$archive" wt_08f_isc_ucell_review
if command -v sha256sum >/dev/null; then sha256sum "$archive"; fi
echo "Review pack: $archive"
echo "First table inside pack: wt_08f_isc_ucell_review/tables/01_proxy_and_union_15_10_05_viable_peak_counts.csv"
