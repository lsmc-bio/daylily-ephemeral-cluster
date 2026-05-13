#!/usr/bin/env python3
"""Infer biological sex from staged FASTQs on the headnode and rerun dry-run.

The remote script samples a small number of read pairs per sample, maps them
against a chr1/chrX/chrY mini-reference, updates the active samples.tsv, and
runs the requested DayOA dry-run command again.
"""

from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

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
import tempfile
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
READ_PAIRS = int(os.environ.get("DAYLILY_SEX_READ_PAIRS", "10000"))
THREADS = int(os.environ.get("DAYLILY_SEX_BWA_THREADS", "2"))
CONTIGS = ("chr1", "chrX", "chrY")
DY_COMMAND = "bin/day_run produce_alignstats produce_multiqc_final produce_snv_concordances produce_tiddit --config aligners=['sent'] dedupers=['dppl'] snv_callers=['sentd'] sv_callers=['tiddit'] -p -j 100 -k -n"


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


def ensure_mini_reference() -> Path:
    WORKDIR.mkdir(parents=True, exist_ok=True)
    mini = WORKDIR / "hg38.chr1_chrX_chrY.fa"
    if mini.exists() and all((mini.parent / (mini.name + ext)).exists() for ext in (".amb", ".ann", ".bwt", ".pac", ".sa")):
        return mini

    fai = REF.with_suffix(REF.suffix + ".fai")
    if not REF.exists() or not fai.exists():
        raise FileNotFoundError(f"Missing reference or fai: {REF}")

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
        wanted = set(CONTIGS)
        current = None
        copied = 0
        for raw in ref:
            if raw.startswith(">"):
                name = raw[1:].split()[0]
                current = name if name in wanted else None
                copied = 0
                if current:
                    out.write(f">{current}\n")
                continue
            if current:
                seq = raw.strip()
                if not seq:
                    continue
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


def extract_reads(src: str, dst: Path, read_pairs: int) -> int:
    lines_needed = read_pairs * 4
    written_lines = 0
    with gzip.open(src, "rt", errors="replace") as inp, dst.open("w") as out:
        for line in inp:
            if written_lines >= lines_needed:
                break
            out.write(line)
            written_lines += 1
    return written_lines // 4


def classify(c1: int, cx: int, cy: int) -> tuple[str, str, str, int, int, float, float]:
    x_ratio = cx / c1 if c1 else 0.0
    y_ratio = cy / c1 if c1 else 0.0
    informative = c1 + cx + cy
    if informative < 100:
        if cy >= 5:
            return "male", "low", "low_depth_y_signal", 1, 1, x_ratio, y_ratio
        return "female", "low", "low_human_content_no_y_signal", 2, 0, x_ratio, y_ratio
    if y_ratio >= 0.025 or (cy >= 25 and x_ratio < 0.45):
        return "male", "high", "y_signal", 1, 1, x_ratio, y_ratio
    if y_ratio <= 0.01 and x_ratio >= 0.45:
        return "female", "high", "no_y_and_x_diploid", 2, 0, x_ratio, y_ratio
    if x_ratio < 0.45:
        return "male", "medium", "x_ratio_male_range", 1, 1, x_ratio, y_ratio
    return "female", "medium", "x_ratio_female_range", 2, 0, x_ratio, y_ratio


def infer_one(sample: str, r1: str, r2: str, mini_ref: Path) -> dict[str, object]:
    sample_dir = WORKDIR / "tmp" / sample
    shutil.rmtree(sample_dir, ignore_errors=True)
    sample_dir.mkdir(parents=True, exist_ok=True)
    r1_tmp = sample_dir / "r1.fq"
    r2_tmp = sample_dir / "r2.fq"
    n1 = extract_reads(r1, r1_tmp, READ_PAIRS)
    n2 = extract_reads(r2, r2_tmp, READ_PAIRS)
    count_pairs = min(n1, n2)
    cmd = [str(BWA), "mem", "-t", str(THREADS), str(mini_ref), str(r1_tmp), str(r2_tmp)]
    env = os.environ.copy()
    env.setdefault("SENTIEON_LICENSE", SENTIEON_LICENSE)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    c1 = cx = cy = other = unmapped = 0
    assert proc.stdout is not None
    for line in proc.stdout:
        if not line or line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            continue
        try:
            flag = int(fields[1])
        except ValueError:
            continue
        if flag & 4:
            unmapped += 1
            continue
        rname = fields[2]
        if rname == "chr1":
            c1 += 1
        elif rname == "chrX":
            cx += 1
        elif rname == "chrY":
            cy += 1
        elif rname != "*":
            other += 1
    stderr = proc.stderr.read() if proc.stderr is not None else ""
    rc = proc.wait()
    (WORKDIR / f"{sample}.bwa.stderr.log").write_text(stderr, encoding="utf-8")
    shutil.rmtree(sample_dir, ignore_errors=True)
    if rc != 0:
        raise RuntimeError(f"bwa mem failed for {sample} with rc={rc}; see {WORKDIR}/{sample}.bwa.stderr.log")
    sex, confidence, reason, nx, ny, x_ratio, y_ratio = classify(c1, cx, cy)
    return {
        "sample": sample,
        "read_pairs_sampled": count_pairs,
        "chr1_reads": c1,
        "chrX_reads": cx,
        "chrY_reads": cy,
        "other_mapped_reads": other,
        "unmapped_reads": unmapped,
        "x_chr1_ratio": round(x_ratio, 6),
        "y_chr1_ratio": round(y_ratio, 6),
        "inferred_biological_sex": sex,
        "confidence": confidence,
        "reason": reason,
        "N_X": nx,
        "N_Y": ny,
    }


def write_status(exit_code: int | None, started_at: str, completed_at: str | None = None) -> None:
    payload = {
        "session_name": SESSION,
        "repo_path": str(REPO),
        "started_at": started_at,
        "completed_at": completed_at,
        "exit_code": exit_code,
        "command": DY_COMMAND,
        "metadata_note": "rerun after sex inference metadata patch",
    }
    (RUN_DIR / "status.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def apply_sex_results(results: list[dict[str, object]]) -> None:
    result_by_sample = {str(row["sample"]): row for row in results}
    for path in (CONFIG_SAMPLES, STAGE_SAMPLES):
        fields, rows = read_tsv(path)
        backup = path.with_name(path.name + f".pre_sexfix_{int(time.time())}")
        shutil.copy2(path, backup)
        for row in rows:
            sample = row["SAMPLEID"]
            result = result_by_sample.get(sample)
            if not result:
                continue
            row["BIOLOGICAL_SEX"] = str(result["inferred_biological_sex"])
            row["N_X"] = str(result["N_X"])
            row["N_Y"] = str(result["N_Y"])
        write_tsv(path, fields, rows)


def rerun_dryrun() -> int:
    log = RUN_DIR / "sexfix_dryrun.log"
    started_at = utc_now()
    write_status(None, started_at)
    script = f"""
set -u
cd {REPO}
. /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate DAYOA
. dyoainit --skip-project-check
set +u
. bin/day_activate slurm hg38 remote
set -u
{DY_COMMAND}
"""
    with log.open("w") as handle:
        proc = subprocess.run(["bash", "-lc", script], stdout=handle, stderr=subprocess.STDOUT, text=True)
    completed_at = utc_now()
    write_status(proc.returncode, started_at, completed_at)
    with log.open("a") as handle:
        handle.write(f"\n[INFO] sexfix dry-run exited with status {proc.returncode}\n")
    return proc.returncode


def main() -> int:
    WORKDIR.mkdir(parents=True, exist_ok=True)
    if not BWA.exists():
        raise FileNotFoundError(f"Missing bwa executable: {BWA}")
    mini_ref = ensure_mini_reference()
    _sample_fields, sample_rows = read_tsv(CONFIG_SAMPLES)
    _unit_fields, unit_rows = read_tsv(UNITS)
    units_by_sample = {row["SAMPLEID"]: row for row in unit_rows}
    results: list[dict[str, object]] = []
    for row in sample_rows:
        sample = row["SAMPLEID"]
        unit = units_by_sample.get(sample)
        if not unit:
            raise RuntimeError(f"No unit row for sample {sample}")
        result = infer_one(sample, unit["ILMN_R1_PATH"], unit["ILMN_R2_PATH"], mini_ref)
        results.append(result)
        print(
            "SEX_INFER\t{sample}\t{inferred_biological_sex}\t{confidence}\t"
            "chr1={chr1_reads}\tchrX={chrX_reads}\tchrY={chrY_reads}\t"
            "x_chr1={x_chr1_ratio}\ty_chr1={y_chr1_ratio}\t{reason}".format(**result),
            flush=True,
        )
    result_fields = [
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
    result_path = WORKDIR / "sex_inference.tsv"
    write_tsv(result_path, result_fields, [{key: str(row[key]) for key in result_fields} for row in results])
    apply_sex_results(results)
    dryrun_rc = rerun_dryrun()
    print(f"SEX_INFERENCE_TSV={result_path}")
    print(f"UPDATED_CONFIG_SAMPLES={CONFIG_SAMPLES}")
    print(f"UPDATED_STAGE_SAMPLES={STAGE_SAMPLES}")
    print(f"SEXFIX_DRYRUN_LOG={RUN_DIR / 'sexfix_dryrun.log'}")
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
        default="/home/ubuntu/daylily-runs/23andme-run1-ILMN/infer_sex_and_rerun.py",
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
    command = f"python3 {args.remote_script}"
    try:
        result = run_shell(
            args.instance_id,
            args.region,
            command,
            profile=args.profile,
            timeout=args.timeout,
            comment="Infer sex and rerun 23andme dry-run",
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
