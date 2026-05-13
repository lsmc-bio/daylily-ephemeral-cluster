#!/usr/bin/env python3
"""Infer sex with one combined mini-reference BWA pass, patch samples.tsv, rerun dry-run."""

from __future__ import annotations

import argparse
import sys
import textwrap

from daylily_ec.aws.ssm import SsmCommandFailedError, run_shell, write_remote_text


REMOTE_SCRIPT = r'''
#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

SESSION = "23andme-run1-ILMN"
RUN_DIR = Path("/home/ubuntu/daylily-runs") / SESSION
REPO = Path("/fsx/analysis_results/ubuntu/23andme-run1-ILMN/daylily-omics-analysis")
STAGE_DIR = Path("/fsx/data/staged_sample_data/remote_stage_20260510T130308Z")
STAGE_SAMPLES = STAGE_DIR / "20260510T130308Z_samples.tsv"
CONFIG_SAMPLES = REPO / "config/samples.tsv"
UNITS = REPO / "config/units.tsv"
WORKDIR = RUN_DIR / "sex_inference"
REF = Path("/fsx/data/genomic_data/organism_references/H_sapiens/hg38/fasta_fai_minalt/GRCh38_no_alt_analysis_set.fasta")
BWA = Path("/fsx/data/cached_envs/sentieon-genomics-202503.02/bin/bwa")
SENTIEON_LICENSE = "/fsx/data/cached_envs/Life_Sciences_Manufacturing_Corporation_eval.lic"
READ_PAIRS = int(os.environ.get("DAYLILY_SEX_READ_PAIRS", "2000"))
THREADS = int(os.environ.get("DAYLILY_SEX_BWA_THREADS", "8"))
S3_BUCKET = "lsmc-dayoa-omics-analysis-us-west-2"
S3_RANGE_BYTES = int(os.environ.get("DAYLILY_SEX_S3_RANGE_BYTES", str(8 * 1024 * 1024)))
CONTIGS = ("chr1", "chrX", "chrY")
DY_COMMAND = "bin/day_run produce_alignstats produce_multiqc_final produce_snv_concordances produce_tiddit --config aligners=['sent'] dedupers=['dppl'] snv_callers=['sentd'] sv_callers=['tiddit'] -p -j 100 -k -n"
RESULT_FIELDS = [
    "sample",
    "read_pairs_sampled",
    "chr1_reads",
    "chrX_reads",
    "chrY_reads",
    "other_mapped_reads",
    "unmapped_reads",
    "x_chr1_ratio",
    "y_chr1_ratio",
    "inferred_biological_sex",
    "confidence",
    "reason",
    "N_X",
    "N_Y",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
    tmp.replace(path)


def backup_file(path: Path) -> Path:
    backup_dir = WORKDIR / "backups"
    backup_dir.mkdir(parents=True, exist_ok=True)
    backup = backup_dir / f"{path.name}.pre_combined_sexfix_{int(time.time())}"
    shutil.copyfile(path, backup)
    return backup


def ensure_mini_reference() -> Path:
    WORKDIR.mkdir(parents=True, exist_ok=True)
    mini = WORKDIR / "hg38.chr1_chrX_chrY.fa"
    index_exts = (".amb", ".ann", ".bwt", ".pac", ".sa")
    if mini.exists() and all((mini.parent / (mini.name + ext)).exists() for ext in index_exts):
        return mini
    fai = REF.with_suffix(REF.suffix + ".fai")
    lengths: dict[str, int] = {}
    with fai.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields and fields[0] in CONTIGS:
                lengths[fields[0]] = int(fields[1])
    missing = sorted(set(CONTIGS) - set(lengths))
    if missing:
        raise RuntimeError(f"Reference fai missing contigs: {missing}")
    with REF.open() as ref, mini.open("w") as out:
        current = None
        copied = 0
        for raw in ref:
            if raw.startswith(">"):
                name = raw[1:].split()[0]
                current = name if name in CONTIGS else None
                copied = 0
                if current:
                    out.write(f">{current}\n")
                continue
            if not current:
                continue
            seq = raw.strip()
            remaining = lengths[current] - copied
            if remaining <= 0:
                current = None
                continue
            if len(seq) > remaining:
                seq = seq[:remaining]
            out.write(seq + "\n")
            copied += len(seq)
    env = os.environ.copy()
    env.setdefault("SENTIEON_LICENSE", SENTIEON_LICENSE)
    subprocess.run([str(BWA), "index", str(mini)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    return mini


def rewrite_header(line: str, sample: str) -> str:
    if not line.startswith("@"):
        return line
    return "@" + sample + "|" + line[1:]


def fsx_path_to_s3_key(path: str) -> str:
    prefix = "/fsx/"
    if not path.startswith(prefix):
        raise ValueError(f"Expected staged /fsx path, got {path}")
    return path[len(prefix):]


def fetch_s3_range(fsx_path: str, dst: Path) -> None:
    key = fsx_path_to_s3_key(fsx_path)
    env = os.environ.copy()
    env.setdefault("AWS_DEFAULT_REGION", "us-west-2")
    env.setdefault("AWS_REGION", "us-west-2")
    cmd = [
        "aws",
        "s3api",
        "get-object",
        "--bucket",
        S3_BUCKET,
        "--key",
        key,
        "--range",
        f"bytes=0-{S3_RANGE_BYTES - 1}",
        str(dst),
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)


def upload_fsx_shadow_to_s3(local_path: Path, fsx_target: Path) -> str:
    s3_uri = f"s3://{S3_BUCKET}/{fsx_path_to_s3_key(str(fsx_target))}"
    env = os.environ.copy()
    env.setdefault("AWS_DEFAULT_REGION", "us-west-2")
    env.setdefault("AWS_REGION", "us-west-2")
    subprocess.run(
        ["aws", "s3", "cp", str(local_path), s3_uri],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
    )
    return s3_uri


def append_reads(src: Path, out, sample: str) -> int:
    max_lines = READ_PAIRS * 4
    line_count = 0
    try:
        with gzip.open(src, "rt", errors="replace") as inp:
            for line in inp:
                if line_count >= max_lines:
                    break
                if line_count % 4 == 0:
                    out.write(rewrite_header(line, sample))
                else:
                    out.write(line)
                line_count += 1
    except EOFError:
        pass
    return line_count // 4


def build_combined_fastqs(unit_rows: list[dict[str, str]]) -> tuple[Path, Path, dict[str, int]]:
    combined_dir = WORKDIR / "combined"
    shutil.rmtree(combined_dir, ignore_errors=True)
    combined_dir.mkdir(parents=True, exist_ok=True)
    r1_out = combined_dir / "all_samples.R1.fq"
    r2_out = combined_dir / "all_samples.R2.fq"
    scratch = combined_dir / "s3_ranges"
    scratch.mkdir(parents=True, exist_ok=True)
    sampled: dict[str, int] = {}
    with r1_out.open("w") as r1_handle, r2_out.open("w") as r2_handle:
        for row in unit_rows:
            sample = row["SAMPLEID"]
            r1_gz = scratch / f"{sample}.R1.fastq.gz"
            r2_gz = scratch / f"{sample}.R2.fastq.gz"
            fetch_s3_range(row["ILMN_R1_PATH"], r1_gz)
            fetch_s3_range(row["ILMN_R2_PATH"], r2_gz)
            n1 = append_reads(r1_gz, r1_handle, sample)
            n2 = append_reads(r2_gz, r2_handle, sample)
            sampled[sample] = min(n1, n2)
            print(f"COMBINED_FASTQ\t{sample}\tread_pairs={sampled[sample]}", flush=True)
            r1_gz.unlink(missing_ok=True)
            r2_gz.unlink(missing_ok=True)
    return r1_out, r2_out, sampled


def run_bwa_counts(mini_ref: Path, r1: Path, r2: Path) -> tuple[dict[str, dict[str, int]], Path]:
    counts: dict[str, dict[str, int]] = {}
    stderr_path = WORKDIR / "combined.bwa.stderr.log"
    env = os.environ.copy()
    env.setdefault("SENTIEON_LICENSE", SENTIEON_LICENSE)
    cmd = [str(BWA), "mem", "-t", str(THREADS), str(mini_ref), str(r1), str(r2)]
    with stderr_path.open("w") as err:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err, text=True, env=env)
        assert proc.stdout is not None
        for line in proc.stdout:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or "|" not in fields[0]:
                continue
            try:
                flag = int(fields[1])
            except ValueError:
                continue
            sample = fields[0].split("|", 1)[0]
            sample_counts = counts.setdefault(
                sample,
                {"chr1": 0, "chrX": 0, "chrY": 0, "other": 0, "unmapped": 0},
            )
            if flag & 4:
                sample_counts["unmapped"] += 1
                continue
            if flag & 0x100 or flag & 0x800:
                continue
            rname = fields[2]
            if rname in ("chr1", "chrX", "chrY"):
                sample_counts[rname] += 1
            elif rname != "*":
                sample_counts["other"] += 1
        rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"combined bwa mem failed rc={rc}; see {stderr_path}")
    return counts, stderr_path


def classify(c1: int, cx: int, cy: int) -> tuple[str, str, str, int, int, float, float]:
    x_ratio = cx / c1 if c1 else 0.0
    y_ratio = cy / c1 if c1 else 0.0
    informative = c1 + cx + cy
    if informative < 50:
        if cy >= 3:
            return "male", "low", "low_depth_y_signal", 1, 1, x_ratio, y_ratio
        return "female", "low", "low_human_content_no_y_signal", 2, 0, x_ratio, y_ratio
    if y_ratio >= 0.025 or cy >= 10:
        return "male", "high", "y_signal", 1, 1, x_ratio, y_ratio
    if y_ratio <= 0.01:
        return "female", "high", "no_y_signal", 2, 0, x_ratio, y_ratio
    return "male", "medium", "weak_y_signal", 1, 1, x_ratio, y_ratio


def patched_sample_rows(path: Path, results: list[dict[str, object]]) -> tuple[list[str], list[dict[str, str]]]:
    result_by_sample = {str(row["sample"]): row for row in results}
    fields, rows = read_tsv(path)
    for row in rows:
        result = result_by_sample.get(row["SAMPLEID"])
        if not result:
            continue
        row["BIOLOGICAL_SEX"] = str(result["inferred_biological_sex"])
        row["N_X"] = str(result["N_X"])
        row["N_Y"] = str(result["N_Y"])
    return fields, rows


def apply_sex_results(results: list[dict[str, object]]) -> None:
    for path in (CONFIG_SAMPLES,):
        fields, rows = patched_sample_rows(path, results)
        backup = backup_file(path)
        write_tsv(path, fields, rows)
        print(f"BACKUP_SAMPLES={backup}", flush=True)
        print(f"UPDATED_CONFIG_SAMPLES={path}", flush=True)

    fields, rows = patched_sample_rows(STAGE_SAMPLES, results)
    backup = backup_file(STAGE_SAMPLES)
    stage_shadow = WORKDIR / STAGE_SAMPLES.name
    write_tsv(stage_shadow, fields, rows)
    print(f"BACKUP_SAMPLES={backup}", flush=True)
    try:
        write_tsv(STAGE_SAMPLES, fields, rows)
        print(f"UPDATED_STAGE_SAMPLES={STAGE_SAMPLES}", flush=True)
    except PermissionError as exc:
        print(f"SKIPPED_STAGE_SAMPLES_FSX_UPDATE={STAGE_SAMPLES}\t{exc}", flush=True)
    try:
        s3_uri = upload_fsx_shadow_to_s3(stage_shadow, STAGE_SAMPLES)
        print(f"UPDATED_STAGE_SAMPLES_S3={s3_uri}", flush=True)
    except subprocess.CalledProcessError as exc:
        print(
            f"SKIPPED_STAGE_SAMPLES_S3_UPDATE={STAGE_SAMPLES}\t"
            f"exit_code={exc.returncode}\tstderr={(exc.stderr or '').strip()}",
            flush=True,
        )


def load_existing_results(path: Path, sample_rows: list[dict[str, str]]) -> list[dict[str, object]] | None:
    if not path.exists():
        return None
    fields, rows = read_tsv(path)
    missing_fields = [field for field in RESULT_FIELDS if field not in fields]
    if missing_fields:
        print(f"IGNORING_EXISTING_SEX_INFERENCE_TSV={path}\tmissing_fields={','.join(missing_fields)}", flush=True)
        return None
    expected_samples = {row["SAMPLEID"] for row in sample_rows}
    observed_samples = {row["sample"] for row in rows}
    if len(rows) != len(expected_samples) or observed_samples != expected_samples:
        print(
            f"IGNORING_EXISTING_SEX_INFERENCE_TSV={path}\t"
            f"rows={len(rows)} expected={len(expected_samples)}",
            flush=True,
        )
        return None
    print(f"REUSING_SEX_INFERENCE_TSV={path}\trows={len(rows)}", flush=True)
    return [{field: row[field] for field in RESULT_FIELDS} for row in rows]


def write_status(exit_code: int | None, started_at: str, completed_at: str | None = None) -> None:
    payload = {
        "session_name": SESSION,
        "repo_path": str(REPO),
        "started_at": started_at,
        "completed_at": completed_at,
        "exit_code": exit_code,
        "command": DY_COMMAND,
        "metadata_note": "rerun after combined sex inference metadata patch",
    }
    (RUN_DIR / "status.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def rerun_dryrun() -> int:
    log = RUN_DIR / "sexfix_tiddit_dryrun.log"
    started_at = utc_now()
    write_status(None, started_at)
    script = f"""
set -u
cd {REPO}
. /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate DAYOA
set +u
. dyoainit --skip-project-check
. bin/day_activate slurm hg38 remote
set -u
{DY_COMMAND}
"""
    with log.open("w") as handle:
        proc = subprocess.run(["bash", "-lc", script], stdout=handle, stderr=subprocess.STDOUT, text=True)
    write_status(proc.returncode, started_at, utc_now())
    with log.open("a") as handle:
        handle.write(f"\n[INFO] combined sexfix dry-run exited with status {proc.returncode}\n")
    return proc.returncode


def main() -> int:
    WORKDIR.mkdir(parents=True, exist_ok=True)
    _sample_fields, sample_rows = read_tsv(CONFIG_SAMPLES)
    _unit_fields, unit_rows = read_tsv(UNITS)
    result_path = WORKDIR / "sex_inference.tsv"
    results = load_existing_results(result_path, sample_rows)
    stderr_path = WORKDIR / "combined.bwa.stderr.log"
    if results is None:
        mini_ref = ensure_mini_reference()
        r1, r2, sampled = build_combined_fastqs(unit_rows)
        counts, stderr_path = run_bwa_counts(mini_ref, r1, r2)
        results = []
        for sample_row in sample_rows:
            sample = sample_row["SAMPLEID"]
            sample_counts = counts.get(sample, {"chr1": 0, "chrX": 0, "chrY": 0, "other": 0, "unmapped": 0})
            sex, confidence, reason, nx, ny, x_ratio, y_ratio = classify(
                sample_counts["chr1"],
                sample_counts["chrX"],
                sample_counts["chrY"],
            )
            result = {
                "sample": sample,
                "read_pairs_sampled": sampled.get(sample, 0),
                "chr1_reads": sample_counts["chr1"],
                "chrX_reads": sample_counts["chrX"],
                "chrY_reads": sample_counts["chrY"],
                "other_mapped_reads": sample_counts["other"],
                "unmapped_reads": sample_counts["unmapped"],
                "x_chr1_ratio": round(x_ratio, 6),
                "y_chr1_ratio": round(y_ratio, 6),
                "inferred_biological_sex": sex,
                "confidence": confidence,
                "reason": reason,
                "N_X": nx,
                "N_Y": ny,
            }
            results.append(result)
            print(
                "SEX_INFER\t{sample}\t{inferred_biological_sex}\t{confidence}\t"
                "chr1={chr1_reads}\tchrX={chrX_reads}\tchrY={chrY_reads}\t"
                "x_chr1={x_chr1_ratio}\ty_chr1={y_chr1_ratio}\t{reason}".format(**result),
                flush=True,
            )
        write_tsv(result_path, RESULT_FIELDS, [{key: str(row[key]) for key in RESULT_FIELDS} for row in results])
    apply_sex_results(results)
    dryrun_rc = rerun_dryrun()
    print(f"SEX_INFERENCE_TSV={result_path}")
    print(f"COMBINED_BWA_STDERR={stderr_path}")
    print(f"UPDATED_CONFIG_SAMPLES={CONFIG_SAMPLES}")
    print(f"STAGE_SAMPLES_TARGET={STAGE_SAMPLES}")
    print(f"SEXFIX_DRYRUN_LOG={RUN_DIR / 'sexfix_tiddit_dryrun.log'}")
    print(f"SEXFIX_DRYRUN_EXIT_CODE={dryrun_rc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
'''


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", default="daylily-service-lsmc")
    parser.add_argument("--region", default="us-west-2")
    parser.add_argument("--instance-id", default="i-047b238434ea745b2")
    parser.add_argument(
        "--remote-script",
        default="/home/ubuntu/daylily-runs/23andme-run1-ILMN/infer_sex_combined_and_rerun.py",
    )
    parser.add_argument("--timeout", type=int, default=7200)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    remote_script = textwrap.dedent(REMOTE_SCRIPT).lstrip()
    write_remote_text(
        args.instance_id,
        args.region,
        args.remote_script,
        remote_script,
        profile=args.profile,
    )
    try:
        result = run_shell(
            args.instance_id,
            args.region,
            f"python3 {args.remote_script}",
            profile=args.profile,
            timeout=args.timeout,
            comment="Infer sex combined and rerun 23andme dry-run",
        )
    except SsmCommandFailedError as exc:
        result = exc.result
        print(result.stdout, end="")
        print(result.stderr, file=sys.stderr, end="")
        return result.response_code or 1
    print(result.stdout, end="")
    print(result.stderr, file=sys.stderr, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
