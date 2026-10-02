#!/usr/bin/env bash
# Run from the repository root on an allocated SOL compute node, after 08L/08N.
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/scratch/dsaiz/multiome_scRNA-ATAC}"
STEP08L="${STEP08L:-/scratch/dsaiz/Yesenia_scData2026/results/v2_resequenced/independent/wt/08l_wt_ISC_open_close_motifs}"
STEP08N="${STEP08N:-/scratch/dsaiz/Yesenia_scData2026/results/v2_resequenced/independent/wt/08n_wt_ISC_interval_loss_pathways}"
OUT="${OUT:-/scratch/dsaiz/Yesenia_scData2026/results/v2_resequenced/independent/wt/08s_shared_loss_HOMER}"
FASTA="${FASTA:-/scratch/dsaiz/Yesenia_scData2026/reference/GRCm39_2024-A/fasta/genome.fa}"
HOMER_BIN="${HOMER_BIN:-/scratch/dsaiz/tools/homer/bin/findMotifsGenome.pl}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

if [[ "$(hostname -s)" == *login* ]]; then
  echo "Run this analysis on an allocated SOL compute node, not a login node." >&2
  exit 1
fi
[[ -r "$FASTA" && -r "${FASTA}.fai" ]] || {
  echo "Missing indexed GRCm39 FASTA: $FASTA" >&2; exit 1;
}
[[ -x "$HOMER_BIN" ]] || {
  echo "HOMER is not installed at: $HOMER_BIN" >&2; exit 1;
}
[[ "$FASTA" == *GRCm39* || "$FASTA" == *mm39* ]] || {
  echo "Expected the matching GRCm39/mm39 reference, not mm10: $FASTA" >&2; exit 1;
}
# findMotifsGenome.pl invokes its companion HOMER programs by name via the
# shell. Calling this one script by absolute path is insufficient unless the
# entire HOMER bin directory is on PATH for its subprocesses as well.
homer_bin_dir="$(cd -- "$(dirname -- "$HOMER_BIN")" && pwd -P)"
export PATH="$homer_bin_dir:$PATH"
for helper in bed2pos.pl checkPeakFile.pl cleanUpPeakFile.pl mergePeaks \
    homerTools cleanUpSequences.pl removePoorSeq.pl freq2group.pl \
    makeBinaryFile.pl homer2 findKnownMotifs.pl compareMotifs.pl; do
  if ! command -v "$helper" >/dev/null 2>&1; then
    echo "Missing HOMER helper: $helper (expected under $homer_bin_dir)" >&2
    echo "Check the HOMER installation before running the enrichment." >&2
    exit 1
  fi
done
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
  echo "THREADS must be a positive integer: $THREADS" >&2; exit 1;
}
mkdir -p "$OUT/inputs"
python3 "$REPO_ROOT/analysis/08s_prepare_shared_homer_inputs.py" \
  "$STEP08N" "$STEP08L" "$OUT/inputs"

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
  if awk -F '\t' 'NR>1 && tolower($1) ~ /yeast|saccharomyces|s[.]cerevisiae|arabidopsis|a[.]thaliana|drosophila|c[.]elegans|plant/ {
    print "Nonvertebrate motif model: " $1 > "/dev/stderr"; bad=1
  } END { exit bad }' "$results/knownResults.txt"; then
    return 0
  fi
  echo "Motif taxonomy check failed: $results/knownResults.txt" >&2
  return 1
}

for direction in open closed; do
  key="shared_HFD_${direction}_lost"
  target="$OUT/inputs/${key}_matched_targets.bed"
  background="$OUT/inputs/${key}_matched_neutral_background.bed"
  results="$OUT/${key}"
  check_bed_coordinates "$target"
  check_bed_coordinates "$background"
  mkdir -p "$results"
  # A previous HOMER attempt may have returned zero despite missing helper
  # commands. Do not let an old knownResults.txt satisfy the success check.
  rm -f -- "$results/knownResults.txt"
  echo "HOMER: $key (target $(wc -l < "$target"), background $(wc -l < "$background"))"
  # -size given preserves the consensus interval lengths; -noweight avoids
  # reweighting controls already matched by 08L on context/GC/width/detection.
  # HOMER still treats the two sets as independent, unlike 08O's paired test.
  "$HOMER_BIN" "$target" "$FASTA" "$results" \
    -bg "$background" -size given -len 8,10,12 \
    -noweight -mset vertebrates -p "$THREADS" \
    > "$results/run.log" 2>&1
  [[ -s "$results/knownResults.txt" ]] || {
    echo "HOMER did not produce knownResults.txt; inspect $results/run.log" >&2
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
  cp "$results/knownResults.txt" \
    "$OUT/01_${key}_HOMER_knownResults.txt"
  if [[ -s "$results/homerResults.html" ]]; then
    echo "De novo HTML: $results/homerResults.html"
  fi
  echo "Known motifs: $OUT/01_${key}_HOMER_knownResults.txt"
done

cat > "$OUT/README_interpretation.txt" <<'EOF'
Target/background coordinates are the Step 08L matched pairs restricted to
Step 08N shared-loss peaks. The custom GRCm39/mm39 FASTA supplies the
matching sequences; -mset vertebrates restricts HOMER's motif collection.
Compare against the version 3 WT 622 baseline run. HOMER known motif statistics use its own library
and enrichment test; they are not equivalent to Step 08O's paired JASPAR
discordance P values. These are DNA sequence tests, not TF occupancy or
genotype-by-diet interaction tests. Examine matched target/background
fractions and motif prevalence before interpreting any TF name.
EOF
echo "HOMER run complete: $OUT"
