#!/usr/bin/env python3
"""Generate 23andMe ILMN per-lane manifests and stage per-lane FASTQs."""

from __future__ import annotations

import argparse
import csv
import os
import sys
import threading
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import boto3
from boto3.s3.transfer import TransferConfig
from botocore.exceptions import ClientError


BUCKET = "lsmc-dayoa-omics-analysis-us-west-2"
RUN_ID_RAW = "20260507_LH01106_0005_A23K3JVLT4"
RUN_ID = "20260507-LH01106-0005-A23K3JVLT4"
EXPERIMENT_ID = "23ame-ilmn-1-bylane"
SMOKE_SAMPLES = {"LSMC_HG003", "LSMC_NTC", "s74478"}
GIAB_HG003 = (
    "/fsx/data/genomic_data/organism_annotations/H_sapiens/hg38/controls/giab/"
    "snv/v4.2.1/HG003/"
)

SAMPLES_HEADER = [
    "SAMPLEID",
    "SAMPLESOURCE",
    "SAMPLECLASS",
    "BIOLOGICAL_SEX",
    "CONCORDANCE_CONTROL_PATH",
    "IS_POSITIVE_CONTROL",
    "IS_NEGATIVE_CONTROL",
    "SAMPLE_TYPE",
    "TUM_NRM_SAMPLEID_MATCH",
    "EXTERNAL_SAMPLE_ID",
    "N_X",
    "N_Y",
    "TRUTH_DATA_DIR",
]

UNITS_HEADER = [
    "RUNID",
    "SAMPLEID",
    "EXPERIMENTID",
    "LANEID",
    "BARCODEID",
    "LIBPREP",
    "SEQ_VENDOR",
    "SEQ_PLATFORM",
    "ILMN_R1_PATH",
    "ILMN_R2_PATH",
    "PACBIO_R1_PATH",
    "PACBIO_R2_PATH",
    "ONT_R1_PATH",
    "ONT_R2_PATH",
    "UG_R1_PATH",
    "UG_R2_PATH",
    "SUBSAMPLE_PCT",
    "ILMN_TRIM_READ_LENGTH",
    "SAMPLEUSE",
    "BWA_KMER",
    "DEEP_MODEL",
    "ULTIMA_CRAM",
    "ULTIMA_CRAM_ALIGNER",
    "ULTIMA_CRAM_SNV_CALLER",
    "ONT_CRAM",
    "ONT_CRAM_ALIGNER",
    "ONT_CRAM_SNV_CALLER",
    "PB_BAM",
    "PB_BAM_ALIGNER",
    "PB_BAM_SNV_CALLER",
    "ROCHE_BAM",
    "ROCHE_BAM_ALIGNER",
    "ROCHE_BAM_SNV_CALLER",
    "ROCHE_DOWNSAMPLE_RATIO",
    "LONGREADTRIM_READ_LENGTH",
    "LONGREADTRIM_MODE",
    "ULTIMA_SUBSAMPLE_PCT",
    "ONT_SUBSAMPLE_PCT",
    "ONT_BAM",
    "ONT_BAM_ALIGNER",
    "ONT_BAM_SNV_CALLER",
]


def parse_s3(uri: str) -> tuple[str, str]:
    if not uri.startswith("s3://"):
        raise ValueError(f"not an S3 URI: {uri}")
    rest = uri[5:]
    bucket, key = rest.split("/", 1)
    return bucket, key


def normalize_identifier(value: str) -> str:
    return value.replace("_", "-").replace(".", "-")


def fsx_to_s3(path: str) -> str:
    if not path.startswith("/fsx/data/"):
        raise ValueError(f"expected /fsx/data path: {path}")
    return f"s3://{BUCKET}/data/{path[len('/fsx/data/'):]}"


def sample_id_for(row: dict[str, str]) -> str:
    return normalize_identifier(row["SAMPLE_ID"])


def lane_dir(row: dict[str, str], stage_fsx: str) -> str:
    lane = int(row["LANE"])
    sample = sample_id_for(row)
    prefix = (
        f"{RUN_ID}_{sample}-NOVASEQ-PCR-FREE-gdna-"
        f"{EXPERIMENT_ID}_D0_L{lane:03d}"
    )
    return f"{stage_fsx}/{prefix}"


def lane_fastq_path(row: dict[str, str], stage_fsx: str, read: int) -> str:
    lane = int(row["LANE"])
    sample = sample_id_for(row)
    prefix = (
        f"{RUN_ID}_{sample}-NOVASEQ-PCR-FREE-gdna-"
        f"{EXPERIMENT_ID}_D0_L{lane:03d}"
    )
    return f"{lane_dir(row, stage_fsx)}/{prefix}_R{read}.fastq.gz"


def sample_row(row: dict[str, str]) -> dict[str, str]:
    raw_sample = row["SAMPLE_ID"]
    sample = sample_id_for(row)
    is_hg003 = raw_sample == "LSMC_HG003"
    is_ntc = raw_sample == "LSMC_NTC"
    return {
        "SAMPLEID": sample,
        "SAMPLESOURCE": row.get("SAMPLESOURCE") or "blood",
        "SAMPLECLASS": row.get("SAMPLECLASS") or "research",
        "BIOLOGICAL_SEX": "male" if is_hg003 else ("female" if is_ntc else "na"),
        "CONCORDANCE_CONTROL_PATH": GIAB_HG003 if is_hg003 else "na",
        "IS_POSITIVE_CONTROL": "true" if is_hg003 else "false",
        "IS_NEGATIVE_CONTROL": "true" if is_ntc else "false",
        "SAMPLE_TYPE": "gdna",
        "TUM_NRM_SAMPLEID_MATCH": row.get("TUM_NRM_SAMPLEID_MATCH") or "na",
        "EXTERNAL_SAMPLE_ID": "HG003" if is_hg003 else ("LSMC_NTC" if is_ntc else raw_sample),
        "N_X": "1" if is_hg003 else ("2" if is_ntc else "na"),
        "N_Y": "1" if is_hg003 else ("0" if is_ntc else "na"),
        "TRUTH_DATA_DIR": GIAB_HG003 if is_hg003 else "na",
    }


def unit_row(row: dict[str, str], stage_fsx: str) -> dict[str, str]:
    unit = {column: "" for column in UNITS_HEADER}
    unit.update(
        {
            "RUNID": RUN_ID,
            "SAMPLEID": sample_id_for(row),
            "EXPERIMENTID": EXPERIMENT_ID,
            "LANEID": str(int(row["LANE"])),
            "BARCODEID": normalize_identifier(row["SEQBC_ID"]),
            "LIBPREP": "PCR-FREE",
            "SEQ_VENDOR": "ILMN",
            "SEQ_PLATFORM": "NOVASEQ",
            "ILMN_R1_PATH": lane_fastq_path(row, stage_fsx, 1),
            "ILMN_R2_PATH": lane_fastq_path(row, stage_fsx, 2),
            "SUBSAMPLE_PCT": row.get("SUBSAMPLE_PCT") or "na",
            "ILMN_TRIM_READ_LENGTH": row.get("ILMN_TRIM_READ_LENGTH", ""),
            "SAMPLEUSE": row["SAMPLEUSE"],
            "BWA_KMER": row.get("BWA_KMER") or "19",
            "DEEP_MODEL": row.get("DEEP_MODEL") or "WGS",
        }
    )
    return unit


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 768:
        raise SystemExit(f"Expected 768 input rows, saw {len(rows)}")
    if any("Undetermined" in row["SAMPLE_ID"] for row in rows):
        raise SystemExit("Unexpected Undetermined row in input manifest")
    return rows


def validate_input(rows: list[dict[str, str]]) -> None:
    by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_sample[row["SAMPLE_ID"]].append(row)
    if len(by_sample) != 96:
        raise SystemExit(f"Expected 96 samples, saw {len(by_sample)}")
    for sample, sample_rows in sorted(by_sample.items()):
        lanes = sorted(int(row["LANE"]) for row in sample_rows)
        if lanes != list(range(1, 9)):
            raise SystemExit(f"{sample} lanes were {lanes}, expected 1..8")
        for row in sample_rows:
            if not row["ILMN_R1_FQ"].startswith("s3://") or not row["ILMN_R2_FQ"].startswith("s3://"):
                raise SystemExit(f"{sample} lane {row['LANE']} has missing FASTQ URI")


def write_tsv(path: Path, rows: list[dict[str, str]], header: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_outputs(
    rows: list[dict[str, str]],
    *,
    smoke_stage_fsx: str,
    full_stage_fsx: str,
) -> dict[str, list[dict[str, str]]]:
    by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_sample[row["SAMPLE_ID"]].append(row)

    sample_rows = [sample_row(sorted(by_sample[sample], key=lambda item: int(item["LANE"]))[0]) for sample in sorted(by_sample)]
    full_units = [unit_row(row, full_stage_fsx) for row in sorted(rows, key=lambda item: (sample_id_for(item), int(item["LANE"])))]
    smoke_source_rows = [
        row
        for row in rows
        if row["SAMPLE_ID"] in SMOKE_SAMPLES
    ]
    smoke_samples = [
        sample_row(sorted(by_sample[sample], key=lambda item: int(item["LANE"]))[0])
        for sample in sorted(SMOKE_SAMPLES)
    ]
    smoke_units = [
        unit_row(row, smoke_stage_fsx)
        for row in sorted(smoke_source_rows, key=lambda item: (sample_id_for(item), int(item["LANE"])))
    ]
    return {
        "full_samples": sample_rows,
        "full_units": full_units,
        "smoke_samples": smoke_samples,
        "smoke_units": smoke_units,
        "smoke_source_rows": smoke_source_rows,
    }


class S3Stager:
    def __init__(self, workers: int) -> None:
        self.s3 = boto3.client("s3", region_name="us-west-2")
        self.transfer_config = TransferConfig(
            multipart_threshold=64 * 1024 * 1024,
            multipart_chunksize=128 * 1024 * 1024,
            max_concurrency=8,
            use_threads=True,
        )
        self.workers = workers
        self.lock = threading.Lock()
        self.done = 0
        self.bytes_copied_or_present = 0

    def object_size(self, uri: str) -> int:
        bucket, key = parse_s3(uri)
        return int(self.s3.head_object(Bucket=bucket, Key=key)["ContentLength"])

    def destination_complete(self, uri: str, expected_size: int) -> bool:
        bucket, key = parse_s3(uri)
        try:
            size = int(self.s3.head_object(Bucket=bucket, Key=key)["ContentLength"])
        except ClientError as exc:
            code = exc.response.get("Error", {}).get("Code")
            if code in {"404", "NoSuchKey", "NotFound"}:
                return False
            raise
        return size == expected_size

    def copy_one(self, source: str, dest: str, total: int) -> str:
        expected = self.object_size(source)
        if self.destination_complete(dest, expected):
            copied = False
        else:
            src_bucket, src_key = parse_s3(source)
            dst_bucket, dst_key = parse_s3(dest)
            self.s3.copy(
                {"Bucket": src_bucket, "Key": src_key},
                dst_bucket,
                dst_key,
                Config=self.transfer_config,
            )
            copied = True
            if not self.destination_complete(dest, expected):
                raise RuntimeError(f"destination size mismatch after copy: {dest}")
        with self.lock:
            self.done += 1
            self.bytes_copied_or_present += expected
            done = self.done
        action = "copied" if copied else "skip"
        return f"[{done}/{total}] {action} bytes={expected} {dest}"

    def upload_file(self, path: Path, dest: str) -> None:
        bucket, key = parse_s3(dest)
        self.s3.upload_file(str(path), bucket, key)

    def stage(self, jobs: list[tuple[str, str]]) -> None:
        total = len(jobs)
        with ThreadPoolExecutor(max_workers=self.workers) as pool:
            futures = [pool.submit(self.copy_one, source, dest, total) for source, dest in jobs]
            for future in as_completed(futures):
                print(future.result(), flush=True)


def validate_outputs(outputs: dict[str, list[dict[str, str]]]) -> None:
    if len(outputs["full_samples"]) != 96:
        raise SystemExit("full_samples row count mismatch")
    if len(outputs["full_units"]) != 768:
        raise SystemExit("full_units row count mismatch")
    if len(outputs["smoke_samples"]) != 3:
        raise SystemExit("smoke_samples row count mismatch")
    if len(outputs["smoke_units"]) != 24:
        raise SystemExit("smoke_units row count mismatch")
    for name in ("full_units", "smoke_units"):
        units = outputs[name]
        if any(row["LANEID"] == "0" for row in units):
            raise SystemExit(f"{name} contains LANEID=0")
        if any("Undetermined" in row["SAMPLEID"] for row in units):
            raise SystemExit(f"{name} contains Undetermined")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--timestamp", default=time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()))
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument(
        "--stage-scope",
        choices=("smoke", "full", "both", "none"),
        default="smoke",
        help="Which FASTQ set to stage after writing manifests.",
    )
    parser.add_argument("--no-stage", action="store_true", help="Alias for --stage-scope none.")
    args = parser.parse_args()

    os.environ.setdefault("AWS_PROFILE", "daylily-service-lsmc")
    os.environ.setdefault("AWS_REGION", "us-west-2")
    os.environ.setdefault("AWS_DEFAULT_REGION", "us-west-2")

    smoke_stage = f"remote_stage_{args.timestamp}_23ame_ilmn_1_bylane_smoke3"
    full_stage = f"remote_stage_{args.timestamp}_23ame_ilmn_1_bylane_full"
    smoke_stage_fsx = f"/fsx/data/staged_sample_data/{smoke_stage}"
    full_stage_fsx = f"/fsx/data/staged_sample_data/{full_stage}"
    smoke_stage_s3 = f"s3://{BUCKET}/data/staged_sample_data/{smoke_stage}"
    full_stage_s3 = f"s3://{BUCKET}/data/staged_sample_data/{full_stage}"

    rows = read_manifest(args.input)
    validate_input(rows)
    outputs = build_outputs(rows, smoke_stage_fsx=smoke_stage_fsx, full_stage_fsx=full_stage_fsx)
    validate_outputs(outputs)

    paths = {
        "full_samples": args.out_dir / "full_samples.tsv",
        "full_units": args.out_dir / "full_units.tsv",
        "smoke_samples": args.out_dir / "smoke3_samples.tsv",
        "smoke_units": args.out_dir / "smoke3_units.tsv",
    }
    write_tsv(paths["full_samples"], outputs["full_samples"], SAMPLES_HEADER)
    write_tsv(paths["full_units"], outputs["full_units"], UNITS_HEADER)
    write_tsv(paths["smoke_samples"], outputs["smoke_samples"], SAMPLES_HEADER)
    write_tsv(paths["smoke_units"], outputs["smoke_units"], UNITS_HEADER)

    smoke_jobs: list[tuple[str, str]] = []
    for row in outputs["smoke_source_rows"]:
        smoke_jobs.append((row["ILMN_R1_FQ"], fsx_to_s3(lane_fastq_path(row, smoke_stage_fsx, 1))))
        smoke_jobs.append((row["ILMN_R2_FQ"], fsx_to_s3(lane_fastq_path(row, smoke_stage_fsx, 2))))

    full_jobs: list[tuple[str, str]] = []
    for row in rows:
        full_jobs.append((row["ILMN_R1_FQ"], fsx_to_s3(lane_fastq_path(row, full_stage_fsx, 1))))
        full_jobs.append((row["ILMN_R2_FQ"], fsx_to_s3(lane_fastq_path(row, full_stage_fsx, 2))))

    stage_scope = "none" if args.no_stage else args.stage_scope
    stager = S3Stager(args.workers)

    smoke_bytes = 0
    for source, _dest in smoke_jobs:
        smoke_bytes += stager.object_size(source)

    full_bytes = 0
    if stage_scope in {"full", "both"}:
        for source, _dest in full_jobs:
            full_bytes += stager.object_size(source)

    if stage_scope in {"smoke", "both"}:
        stager.stage(smoke_jobs)
        for key, path in paths.items():
            stager.upload_file(path, f"{smoke_stage_s3}/manifests/{path.name}")
            print(f"uploaded_manifest {key} {smoke_stage_s3}/manifests/{path.name}", flush=True)

    if stage_scope in {"full", "both"}:
        stager.stage(full_jobs)
        for key in ("full_samples", "full_units"):
            path = paths[key]
            stager.upload_file(path, f"{full_stage_s3}/manifests/{path.name}")
            print(f"uploaded_manifest {key} {full_stage_s3}/manifests/{path.name}", flush=True)

    print(f"timestamp={args.timestamp}")
    print(f"smoke_stage_s3={smoke_stage_s3}/")
    print(f"smoke_stage_fsx={smoke_stage_fsx}/")
    print(f"full_stage_s3={full_stage_s3}/")
    print(f"full_stage_fsx={full_stage_fsx}/")
    print(f"full_samples={paths['full_samples']}")
    print(f"full_units={paths['full_units']}")
    print(f"smoke_samples={paths['smoke_samples']}")
    print(f"smoke_units={paths['smoke_units']}")
    print(f"full_sample_rows={len(outputs['full_samples'])}")
    print(f"full_unit_rows={len(outputs['full_units'])}")
    print(f"smoke_sample_rows={len(outputs['smoke_samples'])}")
    print(f"smoke_unit_rows={len(outputs['smoke_units'])}")
    print(f"stage_scope={stage_scope}")
    print(f"smoke_fastq_objects={len(smoke_jobs)}")
    print(f"smoke_fastq_bytes={smoke_bytes}")
    print(f"smoke_fastq_gib={smoke_bytes / 1024 / 1024 / 1024:.6f}")
    print(f"full_fastq_objects={len(full_jobs)}")
    print(f"full_fastq_bytes={full_bytes}")
    print(f"full_fastq_gib={full_bytes / 1024 / 1024 / 1024:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
