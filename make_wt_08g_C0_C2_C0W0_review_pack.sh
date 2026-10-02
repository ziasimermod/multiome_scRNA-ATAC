#!/usr/bin/env bash
# Package completed WT Step 08G for interpretation and reproducible handoff.
set -euo pipefail

usage() {
  cat >&2 <<'EOF'
Usage: bash make_wt_08g_C0_C2_C0W0_review_pack.sh [WT_OUTPUT_DIR] [PACK_DIR]
Defaults: WT_OUTPUT_DIR=/scratch/dsaiz/Yesenia_scData2026/results/v2_resequenced/independent/wt
          PACK_DIR=WT_OUTPUT_DIR/review_packs
Place this script in the repository root and the companion R script in analysis/.
EOF
}
if (( $# > 2 )); then usage; exit 2; fi
if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then usage; exit 0; fi
wt_root="${1:-/scratch/dsaiz/Yesenia_scData2026/results/v2_resequenced/independent/wt}"
pack_dir="${2:-$wt_root/review_packs}"
step_dir="$wt_root/08g_C0_C2_and_C0W0_ISC"
prior_dir="$wt_root/08f_isc_ucell_wt_crypt"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
summarizer="$script_dir/analysis/08g_summarize_C0_C2_C0W0_review.R"
notebook="$script_dir/analysis/08g_disentangle_C0_C2_and_C0W0_ISC_accessibility.Rmd"

[[ -s "$step_dir/STEP_08G_COMPLETE.txt" ]] || {
  echo "Missing completed Step 08G output: $step_dir/STEP_08G_COMPLETE.txt" >&2; exit 1;
}
[[ -s "$prior_dir/08_Gourab_gene_assay_audit.csv" ]] || {
  echo "Missing post hoc source audit: $prior_dir/08_Gourab_gene_assay_audit.csv" >&2; exit 1;
}
[[ -s "$summarizer" ]] || {
  echo "Missing companion R script: $summarizer" >&2; exit 1;
}
command -v Rscript >/dev/null || { echo 'Rscript is not on PATH.' >&2; exit 1; }
command -v tar >/dev/null || { echo 'tar is not on PATH.' >&2; exit 1; }
mkdir -p -- "$pack_dir"
stage="$(mktemp -d "$pack_dir/.wt08g_review.XXXXXXXX")"
temp_archive="$(mktemp "$pack_dir/.wt08g_archive.XXXXXXXX")"
trap 'rm -rf -- "$stage"; rm -f -- "$temp_archive"' EXIT
content="$stage/wt_08g_C0_C2_C0W0_review"
mkdir -p "$content/tables" "$content/08g/plots" "$content/08f_context" "$content/scripts"

Rscript "$summarizer" "$step_dir" "$prior_dir" "$content/tables"

required=(
  STEP_08G_COMPLETE.txt
  input_manifest.csv
  sessionInfo.txt
  00_matrix_alignment.csv
  01_UCell_score_scale_explanation.csv
  01_core_marker_detection_by_state_diet.csv
  01_score_by_detected_genes_and_state.csv
  02_C0W0_gate_size_and_score_sensitivity.csv
  02_all_nuclei_fixed_labels_and_C0W0_gates.csv.gz
  02_C0W0_selected_identity_and_depth_by_diet.csv
  03_contrast_definitions_and_diet_balance.csv
  04_RNA_evaluable_pooled_effects.csv.gz
  04_ATAC_evaluable_pooled_effects.csv.gz
  04_descriptive_effect_counts.csv
  05_C0W0_top15_peak_candidates_with_all_contrasts.csv.gz
  05_candidate_peak_support_tiers.csv
  05_C0_C2_at_least_one_clear_peak.csv.gz
  05_C0_vs_C2_peak_overlap_summary.csv
  05_ATAC_peak_counts_and_diet_balance.csv
  06_top200_ATAC_equal_cell_and_depth_sensitivity.csv.gz
  06_C0W0_top15_candidates_with_depth_check.csv.gz
  07_peak_promoter_RNA_and_Gourab_overlay.csv.gz
  07_example_loci_selected_without_Gourab_overlay.csv
  07_example_gene_per_nucleus_normalized_RNA.csv.gz
  08_example_coverage_status.csv
  09_motif_input_manifest_not_enrichment.csv
)
for item in "${required[@]}"; do
  [[ -s "$step_dir/$item" ]] || {
    echo "Missing/empty required 08G output: $step_dir/$item" >&2; exit 1;
  }
  cp -- "$step_dir/$item" "$content/08g/$item"
done

plots=(
  01_score_vs_detected_ISC_genes.pdf
  05_evaluable_peak_counts.pdf
  07_examples_RNA_and_ATAC_effects.pdf
  07_example_gene_RNA_violins.pdf
  08_C0_C2_and_C0W0_on_global_WNN.pdf
)
for item in "${plots[@]}"; do
  [[ -s "$step_dir/plots/$item" ]] || {
    echo "Missing required 08G figure: $step_dir/plots/$item" >&2; exit 1;
  }
  cp -- "$step_dir/plots/$item" "$content/08g/plots/$item"
done

# Coverage PDFs are conditional. Their successes/failures remain explicit
# in 08_example_coverage_status.csv. BED files may be zero length.
shopt -s nullglob
for file in "$step_dir"/plots/08_*_example_*.pdf; do
  [[ -s "$file" ]] && cp -- "$file" "$content/08g/plots/"
done
shopt -u nullglob
for item in 09_HFD_higher_candidate_regions.bed \
            09_CON_higher_candidate_regions.bed \
            09_evaluable_unmatched_universe.bed; do
  [[ -f "$step_dir/$item" ]] || {
    echo "Missing 08G BED input (an empty file is valid): $step_dir/$item" >&2; exit 1;
  }
  cp -- "$step_dir/$item" "$content/08g/$item"
done

cp -- "$prior_dir/08_Gourab_gene_assay_audit.csv" "$content/08f_context/"
cp -- "$summarizer" "$content/scripts/"
if [[ -s "$notebook" ]]; then cp -- "$notebook" "$content/scripts/"; fi
stamp="$(date -u +%Y-%m-%d_%H%M%S_UTC)"
archive="$pack_dir/wt_08g_C0_C2_C0W0_review_${stamp}.tar.gz"
[[ ! -e "$archive" ]] || { echo "Archive already exists: $archive" >&2; exit 1; }
tar -C "$stage" -czf "$temp_archive" wt_08g_C0_C2_C0W0_review
mv -- "$temp_archive" "$archive"
if command -v sha256sum >/dev/null; then sha256sum "$archive"; fi
echo "Review pack: $archive"
echo 'Start with: wt_08g_C0_C2_C0W0_review/tables/README_review_first.md'
