# Yesenia snMultiome handoff

The GitHub teaching notebooks are the source to follow. Existing results are
reference provenance, not Yesenia's new analysis or mandatory threshold choices.

## Main workflow

Run analysis/00 through analysis/07 in order, reviewing QC, clustering and
annotation decisions for each cohort. analysis/07d_whole_wnn_cluster_diet_response.Rmd
is the independent whole-WNN-cluster diet-response companion for WT, PPAR and IL17.
The pipeline supplies descriptive effects; pooled nuclei do not create biological
replicates. Cluster IDs and QC decisions are specific to each cohort.

## Optional work and presentation

05B is the IL17 annotation review. 05C-05F, 07B-07C, and 08/08B-08T are optional
reference mapping, refinement, trajectory or locus followups. They remain available
under analysis/ and are indexed in docs/HANDOFF_SCRIPT_REGISTRY.tsv. Follow each
notebook's input requirements; the lettered filenames do not define one mandatory
linear teaching sequence. 07B includes the d21 reference lineage audit.

analysis/presentation_multicohort_diet_program_figures.Rmd is presentation output,
not step 09. The meeting-figure and review-pack scripts are also support tools.
No new scientific step 08 or step 09 is invented by this cleanup.

## Run location

config/project_config.R accepts MULTIOME_RESULTS_DIR for an independent output root:

```r
Sys.setenv(MULTIOME_RESULTS_DIR = "/scratch/ybarrer1/Yesenia_scData2026/results",
           MULTIOME_DATASET_VERSION = "v2_resequenced", MULTIOME_COHORT_ID = "wt")
```

Raw input sample paths and references still point to their shared recorded locations.
Use reviewed config CSVs as starting references, review each decision and record your
own changes. Restart R before switching roots/cohorts. By default, fresh runs write
to the original project's now-empty results/ tree. Existing checkpoints are not
silently used as the starting point of a new run.

Some optional shell and Python tools use explicit path arguments or named environment
overrides rather than project_config.R; supply their input/output paths as documented
by their own --help/usage. RDS objects can retain original absolute fragment or cached
paths internally; source raw inputs remain in place, and no serialized object is rewritten.

## Provenance

Reference location: `/scratch/dsaiz/Yesenia_scData2026/provenance/yesenia_handoff`

- outputs/: latest same-family copies with clean names, including all distinct
  dataset generations, cohorts and checkpoints. v1 versus v2 datasets remain distinct.
- rendered_notebooks/: available rendered run records moved out of the repository.
- source_at_capture/: exact pre-cleanup scripts/configuration, including selected versions.
- output_manifest.tsv: original/current path, bytes, mtime, SHA256 where requested,
  and a reason when not hashed. This is the authoritative rename crosswalk.
- result_cleanup_plan.json and source_cleanup_plan.json: selection and rename rules.
- source_at_capture/config/: QC, clustering and annotation decisions at capture.
- git_state/: commit, branch, status and diff at capture.

Captured scripts/configs are the working tree at handoff. They are not asserted to
have generated every historical output. Latest numbered copy means filename version,
not verified successful execution. Report/log/COMPLETE files are retained as evidence;
check them before interpreting a reference. Original file contents, report text and
embedded historical paths are preserved. Directory/file names alone are normalized.

Older same-family versions, copied review bundles and repository-generated artifacts
are quarantined outside the active repository/provenance in `/scratch/dsaiz/cleanup_backups/yesenia_cleanup_20261002T225830Z`. Different
historical input datasets and legacy analysis branches are not assumed redundant.
The reference tree is separate from active results; no writable output aliases are made.

## Checks before committing

Cleanup verifies source presence, name collisions and file-move identity. It does not
run R, validate scientific thresholds, or claim successful analyses. Review the diff,
run the existing configuration test in your R 4.4 environment, and open notebook 00:

```sh
Rscript tests/test_cohort_configuration.R
git diff --check
git status --short
```

Stage only reviewed code/config/docs with explicit paths. The cleanup never stages,
commits or pushes. Provenance/results remain outside Git. Keep the rollback directory
until review succeeds; rollback before generating new outputs or editing cleaned files.
