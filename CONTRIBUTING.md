# Contributing to the Mana Lab Multiome workflow

The goal is a readable analysis history that a lab member can follow in RStudio.

## Branch and review convention

1. Create a short-lived branch for one bounded analysis change, for example `qc-review`, `common-peaks`, or `colon-annotation`.
2. Keep each pull request focused on one numbered workflow stage or one documented correction.
3. Explain the biological reason for a changed threshold or marker list, not only the line of code that changed.
4. Have at least one lab member review changes to sample metadata, QC thresholds, cluster resolution, or cell-type annotation.
5. Merge only after the corresponding notebook runs from its documented input checkpoint.

## What to commit

- R/R Markdown code and documentation;
- `config/datasets/<dataset_version>/samples.csv` for the current analysis;
- reviewed `config/datasets/<dataset_version>/qc_thresholds.csv` and its decision notes;
- future manual annotation worksheets;
- small summary tables when they are intentionally part of the analysis record.

## What not to commit

- FASTQ, BAM, fragment, H5, or Cell Ranger output directories;
- RDS/Seurat checkpoints;
- rendered notebook caches;
- Slurm logs containing large console output;
- secrets, tokens, or personal credentials.

Before every commit, inspect:

```bash
git status
git diff
```

If a large biological file appears, stop and update `.gitignore` before staging anything.

## Reproducibility note

An analysis-changing decision should leave three records:

1. the configuration or annotation file that stores the chosen value;
2. a notebook or plot that supports the decision;
3. a commit/PR message that explains why it changed.

