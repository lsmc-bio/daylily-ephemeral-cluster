#!/usr/bin/env python3
"""Parallel S3-side staging for the 23andMe ILMN merged FASTQs."""

from __future__ import annotations

import argparse
import csv
import os
import sys
import tempfile
import threading
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable

import boto3
from botocore.exceptions import ClientError


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

RAW_TO_UNIT = {"1": "ILMN_R1_PATH", "2": "ILMN_R2_PATH"}
S3_MULTIPART_MIN_PART_SIZE = 5 * 1024 * 1024


def parse_s3(uri: str) -> tuple[str, str]:
    if not uri.startswith("s3://"):
        raise ValueError(f"not an S3 URI: {uri}")
    rest = uri[5:]
    bucket, key = rest.split("/", 1)
    return bucket, key


def normalise_identifier(value: str) -> str:
    return value.replace("_", "-")


def normalise_run_id(value: str) -> str:
    return normalise_identifier(value).replace(".", "-")


def canonical_libprep(value: str, vendor: str) -> str:
    normalized = normalise_identifier(value)
    if vendor == "ILMN" and normalized.upper() in {
        "NOAMPWGS",
        "PCRFREE",
        "PCR-FREE",
        "TRUSEQPF",
        "TRUSEQ-PF",
    }:
        return "PCR-FREE"
    return normalized


def canonical_platform(value: str, vendor: str) -> str:
    normalized = normalise_identifier(value)
    if vendor == "ILMN" and normalized.upper() in {"NOVASEQ", "NOVASEQX"}:
        return "NOVASEQ"
    return normalized


def headnode_visible(path: str) -> str:
    if path == "/data":
        return "/fsx/data"
    if path.startswith("/data/"):
        return f"/fsx{path}"
    return path


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def grouping_key(row: dict[str, str]) -> tuple[str, ...]:
    vendor = normalise_identifier(row["SEQ_VENDOR"]).upper()
    return (
        normalise_run_id(row["RUN_ID"]),
        normalise_identifier(row["SAMPLE_ID"]),
        normalise_identifier(row["EXPERIMENTID"]),
        normalise_identifier(row["SAMPLE_TYPE"]),
        canonical_libprep(row["LIB_PREP"], vendor),
        vendor,
        canonical_platform(row["SEQ_PLATFORM"], vendor),
    )


def sample_prefix(row: dict[str, str]) -> str:
    run_id, sample_id, experiment_id, sample_type, libprep, vendor, platform = grouping_key(row)
    seqbc = normalise_identifier(row["SEQBC_ID"])
    return f"{run_id}_{sample_id}-{platform}-{libprep}-{sample_type}-{experiment_id}_{seqbc}_0"


def build_outputs(rows: list[dict[str, str]], remote_stage: str) -> tuple[list[dict[str, object]], list[dict[str, str]], list[dict[str, str]]]:
    grouped: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[grouping_key(row)].append(row)

    copy_jobs: list[dict[str, object]] = []
    samples: dict[str, dict[str, str]] = {}
    units: list[dict[str, str]] = []

    for key in sorted(grouped):
        entries = sorted(grouped[key], key=lambda item: int(item["LANE"]))
        if [int(item["LANE"]) for item in entries] != list(range(1, 9)):
            raise SystemExit(f"{key[1]} does not have lanes 1-8")
        first = entries[0]
        prefix = sample_prefix(first)
        dest_fsx_dir = f"{remote_stage}/{prefix}"
        dest_s3_dir = dest_fsx_dir.replace("/data/", "s3://lsmc-dayoa-omics-analysis-us-west-2/data/", 1)
        source_values: dict[str, str] = {}
        for read_no, unit_field in RAW_TO_UNIT.items():
            sources = [entry[f"ILMN_R{read_no}_FQ"] for entry in entries]
            merged_name = f"{prefix}_merged_R{read_no}.fastq.gz"
            dest_s3 = f"{dest_s3_dir}/{merged_name}"
            dest_fsx = f"{dest_fsx_dir}/{merged_name}"
            copy_jobs.append({"sample": key[1], "read": read_no, "sources": sources, "dest_s3": dest_s3})
            source_values[unit_field] = headnode_visible(dest_fsx)

        sample_id = key[1]
        concordance = first["CONCORDANCE_CONTROL_PATH"].strip() or "na"
        samples[sample_id] = {
            "SAMPLEID": sample_id,
            "SAMPLESOURCE": first["SAMPLESOURCE"] or first["SAMPLE_TYPE"],
            "SAMPLECLASS": first["SAMPLECLASS"] or "research",
            "BIOLOGICAL_SEX": first["BIOLOGICAL_SEX"] or "na",
            "CONCORDANCE_CONTROL_PATH": concordance,
            "IS_POSITIVE_CONTROL": first["IS_POS_CTRL"] or "false",
            "IS_NEGATIVE_CONTROL": first["IS_NEG_CTRL"] or "false",
            "SAMPLE_TYPE": normalise_identifier(first["SAMPLE_TYPE"]),
            "TUM_NRM_SAMPLEID_MATCH": first["TUM_NRM_SAMPLEID_MATCH"] or "na",
            "EXTERNAL_SAMPLE_ID": first["EXTERNAL_SAMPLE_ID"] or "na",
            "N_X": first["N_X"] or "na",
            "N_Y": first["N_Y"] or "na",
            "TRUTH_DATA_DIR": concordance,
        }

        unit = {column: "" for column in UNITS_HEADER}
        unit.update(
            {
                "RUNID": key[0],
                "SAMPLEID": sample_id,
                "EXPERIMENTID": key[2],
                "LANEID": "0",
                "BARCODEID": normalise_identifier(first["SEQBC_ID"]),
                "LIBPREP": key[4],
                "SEQ_VENDOR": key[5],
                "SEQ_PLATFORM": key[6],
                "SUBSAMPLE_PCT": first["SUBSAMPLE_PCT"] or "na",
                "ILMN_TRIM_READ_LENGTH": first.get("ILMN_TRIM_READ_LENGTH", ""),
                "SAMPLEUSE": first["SAMPLEUSE"] or ("posControl" if first["IS_POS_CTRL"] == "true" else "sample"),
                "BWA_KMER": first["BWA_KMER"] or "19",
                "DEEP_MODEL": first["DEEP_MODEL"] or "WGS",
                **source_values,
            }
        )
        units.append(unit)

    return copy_jobs, [samples[sample_id] for sample_id in sorted(samples)], units


class FastStager:
    def __init__(self, workers: int) -> None:
        self.s3 = boto3.client("s3", region_name="us-west-2")
        self.workers = workers
        self.meta_cache: dict[str, int] = {}
        self.lock = threading.Lock()
        self.completed = 0

    def object_size(self, uri: str) -> int:
        with self.lock:
            cached = self.meta_cache.get(uri)
        if cached is not None:
            return cached
        bucket, key = parse_s3(uri)
        size = int(self.s3.head_object(Bucket=bucket, Key=key)["ContentLength"])
        with self.lock:
            self.meta_cache[uri] = size
        return size

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

    def copy_job(self, job: dict[str, object], total: int) -> str:
        sources = list(job["sources"])
        dest_s3 = str(job["dest_s3"])
        source_sizes = [self.object_size(source) for source in sources]
        expected_size = sum(source_sizes)
        if self.destination_complete(dest_s3, expected_size):
            with self.lock:
                self.completed += 1
                done = self.completed
            return f"[{done}/{total}] skip {dest_s3}"

        if any(size < S3_MULTIPART_MIN_PART_SIZE for size in source_sizes[:-1]):
            self.local_concat_upload(sources, dest_s3)
            with self.lock:
                self.completed += 1
                done = self.completed
            return f"[{done}/{total}] copied-small {dest_s3} bytes={expected_size}"

        bucket, key = parse_s3(dest_s3)
        upload_id = self.s3.create_multipart_upload(Bucket=bucket, Key=key)["UploadId"]
        parts: list[dict[str, object]] = []
        try:
            for idx, source in enumerate(sources, start=1):
                src_bucket, src_key = parse_s3(source)
                result = self.s3.upload_part_copy(
                    Bucket=bucket,
                    Key=key,
                    PartNumber=idx,
                    UploadId=upload_id,
                    CopySource={"Bucket": src_bucket, "Key": src_key},
                )
                parts.append(
                    {
                        "ETag": result["CopyPartResult"]["ETag"],
                        "PartNumber": idx,
                    }
                )
            self.s3.complete_multipart_upload(
                Bucket=bucket,
                Key=key,
                UploadId=upload_id,
                MultipartUpload={"Parts": parts},
            )
        except Exception:
            self.s3.abort_multipart_upload(Bucket=bucket, Key=key, UploadId=upload_id)
            raise

        with self.lock:
            self.completed += 1
            done = self.completed
        return f"[{done}/{total}] copied {dest_s3} bytes={expected_size}"

    def local_concat_upload(self, sources: list[str], dest_s3: str) -> None:
        bucket, key = parse_s3(dest_s3)
        fd, tmp_name = tempfile.mkstemp(prefix="daylily-small-merge-", suffix=".fastq.gz")
        os.close(fd)
        try:
            with open(tmp_name, "wb") as out_handle:
                for source in sources:
                    src_bucket, src_key = parse_s3(source)
                    self.s3.download_fileobj(src_bucket, src_key, out_handle)
            self.s3.upload_file(tmp_name, bucket, key)
        finally:
            try:
                os.unlink(tmp_name)
            except FileNotFoundError:
                pass

    def run(self, jobs: Iterable[dict[str, object]]) -> None:
        job_list = list(jobs)
        total = len(job_list)
        started = time.time()
        with ThreadPoolExecutor(max_workers=self.workers) as pool:
            futures = [pool.submit(self.copy_job, job, total) for job in job_list]
            for future in as_completed(futures):
                print(future.result(), flush=True)
        elapsed = time.time() - started
        print(f"copy_jobs={total}", flush=True)
        print(f"copy_elapsed_seconds={elapsed:.1f}", flush=True)


def write_tsv(path: Path, header: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def upload_file(local_path: Path, s3_uri: str) -> None:
    bucket, key = parse_s3(s3_uri)
    boto3.client("s3", region_name="us-west-2").upload_file(str(local_path), bucket, key)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--remote-stage", required=True)
    parser.add_argument("--config-dir", required=True)
    parser.add_argument("--workers", type=int, default=24)
    args = parser.parse_args()

    manifest = Path(args.manifest).resolve()
    config_dir = Path(args.config_dir).resolve()
    remote_stage = args.remote_stage.rstrip("/")
    if not remote_stage.startswith("/data/staged_sample_data/"):
        raise SystemExit(f"remote stage must start with /data/staged_sample_data/: {remote_stage}")

    rows = read_manifest(manifest)
    copy_jobs, samples_rows, units_rows = build_outputs(rows, remote_stage)
    if len(copy_jobs) != 192:
        raise SystemExit(f"expected 192 merged FASTQ copy jobs, saw {len(copy_jobs)}")
    if len(samples_rows) != 96:
        raise SystemExit(f"expected 96 samples rows, saw {len(samples_rows)}")
    if len(units_rows) != 96:
        raise SystemExit(f"expected 96 unit rows, saw {len(units_rows)}")

    print(f"remote_fsx_stage={headnode_visible(remote_stage)}", flush=True)
    print(f"merged_fastq_jobs={len(copy_jobs)}", flush=True)
    print(f"samples_rows={len(samples_rows)}", flush=True)
    print(f"units_rows={len(units_rows)}", flush=True)

    FastStager(workers=args.workers).run(copy_jobs)

    timestamp = remote_stage.rsplit("remote_stage_", 1)[-1]
    samples_path = config_dir / f"{timestamp}_samples.tsv"
    units_path = config_dir / f"{timestamp}_units.tsv"
    write_tsv(samples_path, SAMPLES_HEADER, samples_rows)
    write_tsv(units_path, UNITS_HEADER, units_rows)

    remote_s3_stage = remote_stage.replace(
        "/data/",
        "s3://lsmc-dayoa-omics-analysis-us-west-2/data/",
        1,
    )
    upload_file(samples_path, f"{remote_s3_stage}/{samples_path.name}")
    upload_file(units_path, f"{remote_s3_stage}/{units_path.name}")

    print("Remote staging completed successfully.", flush=True)
    print(f"Remote FSx stage directory: {headnode_visible(remote_stage)}", flush=True)
    print("Generated configuration files:", flush=True)
    print(f"  samples.tsv -> {remote_s3_stage}/{samples_path.name}", flush=True)
    print(f"  units.tsv   -> {remote_s3_stage}/{units_path.name}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
