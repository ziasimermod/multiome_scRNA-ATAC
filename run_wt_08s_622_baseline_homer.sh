#!/usr/bin/env bash
# Run on an allocated SOL compute node from the repository root after 08L.
# Uses all frozen WT HFD opening/closing candidates. Distinct 08L matched
# neutral peaks provide the background. HOMER tests independent sets, not pairs.
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/scratch/dsaiz/multiome_scRNA-ATAC}"
STEP08L="${STEP08L:-/scratch/dsaiz/Yesenia_scData2026/results/v2_resequenced/independent/wt/08l_wt_ISC_open_close_motifs}"
OUT="${OUT:-/scratch/dsaiz/Yesenia_scData2026/results/v2_resequenced/independent/wt/08s_wt_622_baseline_HOMER}"
FASTA="${FASTA:-/scratch/dsaiz/Yesenia_scData2026/reference/GRCm39_2024-A/fasta/genome.fa}"
HOMER_BIN="${HOMER_BIN:-/scratch/dsaiz/tools/homer/bin/findMotifsGenome.pl}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

if [[ "$(hostname -s)" == *login* ]]; then
  echo "Run this analysis on an allocated SOL compute node." >&2
  exit 1
fi
[[ -r "$FASTA" && -r "${FASTA}.fai" ]] || {
  echo "Missing indexed GRCm39 FASTA: $FASTA" >&2; exit 1;
}
[[ -x "$HOMER_BIN" ]] || {
  echo "HOMER executable missing: $HOMER_BIN" >&2; exit 1;
}
[[ "$FASTA" == *GRCm39* || "$FASTA" == *mm39* ]] || {
  echo "Expected the matching GRCm39/mm39 reference, not mm10: $FASTA" >&2; exit 1;
}
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
  echo "THREADS must be a positive integer: $THREADS" >&2; exit 1;
}
[[ -f "$REPO_ROOT/analysis/08s_prepare_wt_622_homer_inputs.py" ]] || {
  echo "The versioned 08S v2 Python companion is missing under $REPO_ROOT/analysis" >&2
  exit 1
}

# HOMER calls these programs by name from findMotifsGenome.pl. An absolute
# path for that entry point does not make the helpers visible to its shell.
homer_bin_dir="$(cd -- "$(dirname -- "$HOMER_BIN")" && pwd -P)"
export PATH="$homer_bin_dir:$PATH"
for helper in bed2pos.pl checkPeakFile.pl cleanUpPeakFile.pl mergePeaks \
    homerTools cleanUpSequences.pl removePoorSeq.pl freq2group.pl \
    makeBinaryFile.pl homer2 findKnownMotifs.pl compareMotifs.pl; do
  command -v "$helper" >/dev/null 2>&1 || {
    echo "Missing HOMER helper: $helper (expected under $homer_bin_dir)" >&2
    exit 1
  }
done

mkdir -p "$OUT/inputs"
python3 "$REPO_ROOT/analysis/08s_prepare_wt_622_homer_inputs.py" \
  "$STEP08L" "$OUT/inputs"

cat > "$OUT/00_run_parameters.txt" <<EOF
Assembly: GRCm39/mm39 (BED coordinates and custom FASTA)
FASTA: $FASTA
FASTA index SHA256: $(sha256sum "${FASTA}.fai" | cut -d' ' -f1)
FASTA size (bytes): $(wc -c < "$FASTA")
HOMER executable: $HOMER_BIN
Known and de novo motif comparison set: vertebrates (-mset vertebrates)
Background: Step 08L matched neutral accessible peaks
EOF

check_bed_coordinates() {
  local bed="$1"
  awk -v file="$bed" 'BEGIN { FS="\t"; bad=0 }
    NR==FNR { chrlen[$1]=$2; next }
    !($1 in chrlen) || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ ||
      $2 < 0 || $3 <= $2 || $3 > chrlen[$1] {
        print "Invalid GRCm39 BED interval in " file " at line " FNR ": " $0 > "/dev/stderr";
        bad=1
      }
    END { exit bad }' "${FASTA}.fai" "$bed"
}

check_motif_taxonomy() {
  local results="$1"
  # Check motif labels, not sequence patterns: cross-species PWMs can match
  # mouse DNA, but yeast/plant/fly models are outside the approved library.
  if awk -F '\t' 'NR>1 && tolower($1) ~ /yeast|saccharomyces|s[.]cerevisiae|arabidopsis|a[.]thaliana|drosophila|c[.]elegans|plant/ {
    print "Nonvertebrate motif model: " $1 > "/dev/stderr"; bad=1
  } END { exit bad }' "$results/knownResults.txt"; then
    return 0
  fi
  echo "Motif taxonomy check failed: $results/knownResults.txt" >&2
  return 1
}

for direction in open closed; do
  key="WT_HFD_${direction}"
  target="$OUT/inputs/${key}_all_targets.bed"
  background="$OUT/inputs/${key}_matched_neutral_background.bed"
  results="$OUT/${key}"
  check_bed_coordinates "$target"
  check_bed_coordinates "$background"
  mkdir -p "$results"
  rm -f -- "$results/knownResults.txt"
  echo "HOMER WT ${direction}: $(wc -l < "$target") full-set targets / $(wc -l < "$background") matched neutral controls"
  # -mset restricts HOMER's motif collection. The genome argument selects
  # sequence coordinates; it does not impose motif taxonomy.
  "$HOMER_BIN" "$target" "$FASTA" "$results" \
    -bg "$background" -size given -len 8,10,12 \
    -noweight -mset vertebrates -p "$THREADS" \
    > "$results/run.log" 2>&1
  [[ -s "$results/knownResults.txt" ]] || {
    echo "HOMER did not create knownResults.txt; inspect $results/run.log" >&2
    exit 1
  }
  awk 'END { exit NR < 2 }' "$results/knownResults.txt" || {
    echo "HOMER produced no known motif rows: $results/knownResults.txt" >&2
    exit 1
  }
  if grep -En 'command not found|No such file or directory' "$results/run.log"; then
    echo "HOMER reported missing programs or files in $results/run.log" >&2
    exit 1
  fi
  check_motif_taxonomy "$results"
  cp "$results/knownResults.txt" "$OUT/01_${key}_HOMER_knownResults.txt"
  echo "Known results: $OUT/01_${key}_HOMER_knownResults.txt"
done

cat > "$OUT/README_interpretation.txt" <<'EOF'
The source is the frozen Step 08L set of 622 WT ISC-enriched candidate peaks:
387 HFD-opening and 235 HFD-closing. All 622 enter the foreground. The
background contains distinct Step 08L context/GC/width/WT-CON-detection
matched neutral peaks (378 opening, 230 closing); 9 opening and 5 closing
foreground peaks have no qualifying matched neutral control. The matching
audit identifies them. The custom GRCm39/mm39 FASTA supplies the matching
sequences; -mset vertebrates restricts HOMER's motif library. Compare this
version to the version 3 shared-loss run. This establishes motif enrichment for the
WT candidate sequence sets relative to neutral WT accessible peaks. HOMER
does not test paired discordance, TF occupancy, diet P values, or whether a
motif specifically predicts loss of WT effects in PPAR or IL17. Compare motif
prevalence and model/family, not P values across the overlapping sets.
EOF
echo "WT 622-source HOMER v3 complete: $OUT"
