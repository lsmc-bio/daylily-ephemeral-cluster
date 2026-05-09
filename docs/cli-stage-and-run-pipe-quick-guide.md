# Get Going With `dyec`

This is the short path for staging DayOA test data and launching remote dry-run
analysis checks from an existing Daylily ParallelCluster. It does not create,
export, delete, or tear down AWS resources.

## Setup

Run from the `daylily-ephemeral-cluster` repo root:

```bash
source ./activate

export AWS_PROFILE=daylily-service-lsmc
export REGION=us-west-2
export CLUSTER_NAME=mk-gotime3
export REF_BUCKET=s3://lsmc-dayoa-omics-analysis-us-west-2
export RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
```

Optional sanity checks:

```bash
dyec version
dyec runtime status
dyec headnode jobs --profile "$AWS_PROFILE" --region "$REGION" --cluster "$CLUSTER_NAME"
```

## Test Data Cases

| Case | Manifest | Catalog command |
| --- | --- | --- |
| ONT solo | `examples/staging/ont_solo/analysis_samples_manifest.tsv` | `ont_snv_alignstats` |
| Ultima solo | `examples/staging/ultima_solo/analysis_samples_manifest.tsv` | `ultima_snv_alignstats` |
| ILMN solo | `examples/staging/ilmn_solo/analysis_samples_manifest.tsv` | `illumina_snv_alignstats` |
| ILMN+ONT hybrid | `examples/staging/hybrid_ilmn_ont/analysis_samples_manifest.tsv` | `hybrid_ilmn_ont_snv` |

Use a unique destination and session for every run. The remote launcher clones
into `/fsx/analysis_results/ubuntu/<destination>/...` and fails if that
destination already exists.

## One-Step Stage And Dry-Run

This is the easiest path. It stages the manifest, validates compatibility with
the catalog command, launches the remote DayOA workflow dry-run, and writes a
receipt next to the generated config files.

Set one case:

```bash
export CASE=ont_solo
export MANIFEST=examples/staging/ont_solo/analysis_samples_manifest.tsv
export COMMAND_ID=ont_snv_alignstats
export DESTINATION="stg-ex-${CASE}-${RUN_STAMP}"
export SESSION="$DESTINATION"
export CFG_DIR="$PWD/tmp-stage-config/get-going/${RUN_STAMP}/${CASE}"
```

Then run:

```bash
dyec samples run "$MANIFEST" \
  --command-id "$COMMAND_ID" \
  --profile "$AWS_PROFILE" \
  --region "$REGION" \
  --cluster "$CLUSTER_NAME" \
  --reference-bucket "$REF_BUCKET" \
  --config-dir "$CFG_DIR" \
  --destination "$DESTINATION" \
  --session-name "$SESSION" \
  --git-tag 0.7.726 \
  --dry-run
```

Repeat with the table values for the other cases. For full execution, use a new
destination/session and remove `--dry-run`.

## Two-Step Stage Then Launch

Use this path when you want to inspect generated `samples.tsv` and `units.tsv`
before launching.

Stage:

```bash
dyec samples stage "$MANIFEST" \
  --profile "$AWS_PROFILE" \
  --region "$REGION" \
  --reference-bucket "$REF_BUCKET" \
  --config-dir "$CFG_DIR"
```

Copy the printed `Remote FSx stage directory` into `STAGE_DIR`:

```bash
export STAGE_DIR=/fsx/data/staged_sample_data/remote_stage_<timestamp>
```

Pick the dry-run DayOA command for the case:

```bash
# ONT solo
export DY_COMMAND="bin/day_run produce_alignstats produce_sentmm2ont_align_sort produce_sentdont_vcf produce_snv_concordances --config dedupers=['na'] -p -j 5 -k -n"

# Ultima solo
export DY_COMMAND="bin/day_run produce_alignstats produce_sentdug_vcf produce_snv_concordances --config dppl=['na'] -p -j 20 -k -n"

# ILMN solo
export DY_COMMAND="bin/day_run produce_snv_concordances produce_alignstats --config aligners=['sent'] dedupers=['dppl'] snv_callers=['sentd'] -p -k -j 20 -n"

# ILMN+ONT hybrid
export DY_COMMAND="bin/day_run produce_snv_concordances produce_sentdhiom_sv produce_sentdhiom_vcf -p -j 100 -k -n"
```

Launch:

```bash
dyec workflow launch \
  --profile "$AWS_PROFILE" \
  --region "$REGION" \
  --cluster "$CLUSTER_NAME" \
  --stage-dir "$STAGE_DIR" \
  --destination "$DESTINATION" \
  --session-name "$SESSION" \
  --git-tag 0.7.726 \
  --genome hg38_broad \
  --dy-command "$DY_COMMAND" \
  --dry-run
```

For full execution, use a new destination/session, remove `--dry-run`, and use
the same `DY_COMMAND` without the trailing `-n`.

## Check Status And Logs

```bash
dyec --json workflow status \
  --profile "$AWS_PROFILE" \
  --region "$REGION" \
  --cluster "$CLUSTER_NAME" \
  --session "$SESSION"

dyec workflow logs \
  --profile "$AWS_PROFILE" \
  --region "$REGION" \
  --cluster "$CLUSTER_NAME" \
  --session "$SESSION" \
  --lines 100
```

Successful dry-runs should eventually report `exit_code: 0` in workflow status.
