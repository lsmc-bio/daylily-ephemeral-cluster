# 23andMe ILMN Per-Lane Command Log Summary

Run: `23ame-ilmn-1-bylane`

DayOA checkout:

`/fsx/analysis_results/ubuntu/23ame-ilmn-1-bylane/daylily-omics-analysis`

This note summarizes where the command provenance lives for the per-lane 23andMe Illumina run, and pulls out representative `sentieon bwa mem` and `doppelmark` invocations from the exported Snakemake logs.

## Log Links

The presigned URLs generated during the live handoff were valid for 7 days and
are intentionally not persisted here. Use the durable S3 object paths below to
regenerate temporary access links when needed.

| Log | What it contains | Durable S3 object |
| --- | --- | --- |
| `day_cmd.log` | High-level DayOA launch history, including dry-runs, broad production launches, the alignstats priority run, and the later remaining-rule resumes. | `s3://lsmc-dayoa-omics-analysis-us-west-2/FSxLustre20260509T204833Z/analysis_results/ubuntu/23ame-ilmn-1-bylane/daylily-omics-analysis/day_cmd.log` |
| Full production Snakemake log | Initial broad production DAG and command provenance. This is the source used below for the example `sentieon_bwa_sort` and `doppelmark_dups` command bodies. | `s3://lsmc-dayoa-omics-analysis-us-west-2/FSxLustre20260509T204833Z/analysis_results/ubuntu/23ame-ilmn-1-bylane/daylily-omics-analysis/.snakemake/log/2026-05-11T112232.515957.snakemake.log` |
| Remaining `-j 400` Snakemake log | Later high-concurrency resume for the remaining workflow after alignstats prioritization, including downstream rule progress such as SNV concordance, TIDDIT, VEP, and final reporting. | `s3://lsmc-dayoa-omics-analysis-us-west-2/FSxLustre20260509T204833Z/analysis_results/ubuntu/23ame-ilmn-1-bylane/daylily-omics-analysis/.snakemake/log/2026-05-11T163408.869558.snakemake.log` |

Durable S3 object paths:

```text
s3://lsmc-dayoa-omics-analysis-us-west-2/FSxLustre20260509T204833Z/analysis_results/ubuntu/23ame-ilmn-1-bylane/daylily-omics-analysis/day_cmd.log
s3://lsmc-dayoa-omics-analysis-us-west-2/FSxLustre20260509T204833Z/analysis_results/ubuntu/23ame-ilmn-1-bylane/daylily-omics-analysis/.snakemake/log/2026-05-11T112232.515957.snakemake.log
s3://lsmc-dayoa-omics-analysis-us-west-2/FSxLustre20260509T204833Z/analysis_results/ubuntu/23ame-ilmn-1-bylane/daylily-omics-analysis/.snakemake/log/2026-05-11T163408.869558.snakemake.log
```

## DAG Evidence

The full production Snakemake log includes 768 per-lane alignment jobs and 768 per-lane duplicate-marking jobs:

```text
sentieon_bwa_sort    768    192    192
doppelmark_dups      768    192    192
```

The sample/unit names include `23ame-ilmn-1-bylane-<lane>` and the outputs are written under per-unit result directories. That is the command-level evidence that these results were produced lane-by-lane rather than by merging all lanes for a barcode before alignment.

## Representative Sentieon BWA Command

Rule: `sentieon_bwa_sort`

Example lane unit:

`20260507-LH01106-0005-A23K3JVLT4-s77716-23ame-ilmn-1-bylane-8-D0-PCR-FREE-ILMN-NOVASEQ`

Inputs:

```text
results/day/hg38/20260507-LH01106-0005-A23K3JVLT4-s77716-23ame-ilmn-1-bylane-8-D0-PCR-FREE-ILMN-NOVASEQ/20260507-LH01106-0005-A23K3JVLT4-s77716-23ame-ilmn-1-bylane-8-D0-PCR-FREE-ILMN-NOVASEQ.R1.fastq.gz
results/day/hg38/20260507-LH01106-0005-A23K3JVLT4-s77716-23ame-ilmn-1-bylane-8-D0-PCR-FREE-ILMN-NOVASEQ/20260507-LH01106-0005-A23K3JVLT4-s77716-23ame-ilmn-1-bylane-8-D0-PCR-FREE-ILMN-NOVASEQ.R2.fastq.gz
```

Output:

```text
results/day/hg38/20260507-LH01106-0005-A23K3JVLT4-s77716-23ame-ilmn-1-bylane-8-D0-PCR-FREE-ILMN-NOVASEQ/align/sent/20260507-LH01106-0005-A23K3JVLT4-s77716-23ame-ilmn-1-bylane-8-D0-PCR-FREE-ILMN-NOVASEQ.sent.sort.bam
```

Representative command:

```bash
LD_PRELOAD=$LD_PRELOAD /fsx/data/cached_envs/sentieon-genomics-202503.02/bin/sentieon bwa mem \
  -t 96 -k 19 -Y -M -K 10000000 \
  -x /fsx/data/cached_envs/sentieon-genomics-202503.02/bundles/SentieonIlluminaWGS2.2.bundle/bwa.model \
  -R "@RG\tID:20260507-LH01106-0005-A23K3JVLT4-s77716-23ame-ilmn-1-bylane-8-D0-PCR-FREE-ILMN-NOVASEQ-$epocsec\tSM:20260507-LH01106-0005-A23K3JVLT4-s77716-23ame-ilmn-1-bylane-8-D0-PCR-FREE-ILMN-NOVASEQ\tLB:20260507-LH01106-0005-A23K3JVLT4-s77716-23ame-ilmn-1-bylane-8-D0-PCR-FREE-ILMN-NOVASEQ-LB-1\tPL:ILLUMINA" \
  /fsx/data/genomic_data/organism_references/H_sapiens/hg38/fasta_fai_minalt/GRCh38_no_alt_analysis_set.fasta \
  <(igzip -cd -T 32 -q results/day/hg38/20260507-LH01106-0005-A23K3JVLT4-s77716-23ame-ilmn-1-bylane-8-D0-PCR-FREE-ILMN-NOVASEQ/20260507-LH01106-0005-A23K3JVLT4-s77716-23ame-ilmn-1-bylane-8-D0-PCR-FREE-ILMN-NOVASEQ.R1.fastq.gz) \
  <(igzip -cd -T 32 -q results/day/hg38/20260507-LH01106-0005-A23K3JVLT4-s77716-23ame-ilmn-1-bylane-8-D0-PCR-FREE-ILMN-NOVASEQ/20260507-LH01106-0005-A23K3JVLT4-s77716-23ame-ilmn-1-bylane-8-D0-PCR-FREE-ILMN-NOVASEQ.R2.fastq.gz) \
  | mbuffer -m 128G -q -s 2M \
  | /fsx/data/cached_envs/sentieon-genomics-202503.02/bin/sentieon util sort \
      -t 96 \
      --reference /fsx/data/genomic_data/organism_references/H_sapiens/hg38/fasta_fai_minalt/GRCh38_no_alt_analysis_set.fasta \
      --cram_write_options version=3.0,compressor=rans \
      --sortblock_thread_count 96 \
      --bam_compression 1 \
      --temp_dir "$sort_tmp" \
      --intermediate_compress_level 1 \
      --block_size 4G \
      --sam2bam \
      -o results/day/hg38/20260507-LH01106-0005-A23K3JVLT4-s77716-23ame-ilmn-1-bylane-8-D0-PCR-FREE-ILMN-NOVASEQ/align/sent/20260507-LH01106-0005-A23K3JVLT4-s77716-23ame-ilmn-1-bylane-8-D0-PCR-FREE-ILMN-NOVASEQ.sent.sort.bam -
```

Key interpretation:

- `bwa mem` was run against the GRCh38 no-alt reference.
- The read group `SM` is the per-lane DayOA unit name, not just the 96-sample barcode name.
- The command emits a lane-level coordinate-sorted BAM for downstream duplicate marking and QC.
- The Snakemake resource block for this example requested 192 threads and 300000 MB.

## Representative Doppelmark Command

Rule: `doppelmark_dups`

Example lane unit:

`20260507-LH01106-0005-A23K3JVLT4-s77700-23ame-ilmn-1-bylane-5-D0-PCR-FREE-ILMN-NOVASEQ`

Input:

```text
results/day/hg38/20260507-LH01106-0005-A23K3JVLT4-s77700-23ame-ilmn-1-bylane-5-D0-PCR-FREE-ILMN-NOVASEQ/align/sent/20260507-LH01106-0005-A23K3JVLT4-s77700-23ame-ilmn-1-bylane-5-D0-PCR-FREE-ILMN-NOVASEQ.sent.sort.bam
```

Output:

```text
results/day/hg38/20260507-LH01106-0005-A23K3JVLT4-s77700-23ame-ilmn-1-bylane-5-D0-PCR-FREE-ILMN-NOVASEQ/align/sent/dmd/20260507-LH01106-0005-A23K3JVLT4-s77700-23ame-ilmn-1-bylane-5-D0-PCR-FREE-ILMN-NOVASEQ.sent.dmd.cram
```

Representative command:

```bash
OMP_NUM_THREADS=192 OMP_PROC_BIND=close OMP_PLACES=threads OMP_DYNAMIC=true \
OMP_MAX_ACTIVE_LEVELS=1 OMP_SCHEDULE=dynamic OMP_WAIT_POLICY=ACTIVE \
  resources/DOPPLEMARK/doppelmark \
    -parallelism 192 \
    -bam results/day/hg38/20260507-LH01106-0005-A23K3JVLT4-s77700-23ame-ilmn-1-bylane-5-D0-PCR-FREE-ILMN-NOVASEQ/align/sent/20260507-LH01106-0005-A23K3JVLT4-s77700-23ame-ilmn-1-bylane-5-D0-PCR-FREE-ILMN-NOVASEQ.sent.sort.bam \
    -clip-padding 800 \
    -optical-distance 2500 \
    -logtostderr \
    -disk-mate-shards 0 \
    -scratch-dir "$tdir" \
    -min-bases 20000 \
    -queue-length 1250 \
    -shard-size 50000000 \
  | mbuffer -m 1G \
  | samtools view -@ 8 -m 2G --output-fmt-option level=1 -C \
      -T /fsx/data/genomic_data/organism_references/H_sapiens/hg38/fasta_fai_minalt/GRCh38_no_alt_analysis_set.fasta \
      --write-index \
      -o results/day/hg38/20260507-LH01106-0005-A23K3JVLT4-s77700-23ame-ilmn-1-bylane-5-D0-PCR-FREE-ILMN-NOVASEQ/align/sent/dmd/20260507-LH01106-0005-A23K3JVLT4-s77700-23ame-ilmn-1-bylane-5-D0-PCR-FREE-ILMN-NOVASEQ.sent.dmd.cram -
```

Key interpretation:

- Doppelmark consumed the per-lane Sentieon-sorted BAM.
- The deduplicated output is a per-lane CRAM and CRAI in `align/sent/dmd/`.
- The command keeps duplicate marking scoped to each lane-level unit.
- The Snakemake resource block for this example requested 192 threads and 200000 MB.

## How To Read These Logs

- `day_cmd.log` is the top-level operator history. It answers "what DayOA/Snakemake commands were launched and when?"
- The full production Snakemake log is the best provenance source for early heavy compute jobs, including alignment and duplicate marking command bodies.
- The later `-j 400` Snakemake log is the best provenance source for the high-concurrency resume after the alignstats-priority interruption.
- The Snakemake logs contain rule names, wildcard-expanded sample/lane IDs, inputs, outputs, resource requests, and shell commands.
- The lane-level naming convention is the important provenance detail for this run: outputs include `23ame-ilmn-1-bylane-<lane>`, so each FASTQ lane pair was aligned and duplicate-marked independently.
