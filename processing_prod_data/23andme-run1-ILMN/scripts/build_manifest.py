#!/usr/bin/env python3
"""Build the 23andMe run1 ILMN dyec staging manifest from S3."""

from __future__ import annotations

import csv
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


RUN_ID = "20260507_LH01106_0005_A23K3JVLT4"
EXPERIMENT_ID = "23andme-run1-ILMN"
FASTQ_PREFIX = (
    "s3://lsmc-ssf-sequencing-data/basecalls/lsmc/ssf-hq/LH01106/2026/"
    "20260507_LH01106_0005_A23K3JVLT4/Analysis/1/Data/BCLConvert/fastq/"
)
SAMPLE_SHEET = (
    "s3://lsmc-ssf-sequencing-data/basecalls/lsmc/ssf-hq/LH01106/2026/"
    "20260507_LH01106_0005_A23K3JVLT4/SampleSheet.csv"
)
GIAB_HG003 = (
    "/fsx/data/genomic_data/organism_annotations/H_sapiens/hg38/controls/giab/"
    "snv/v4.2.1/HG003/"
)

HEADER = [
    "RUN_ID",
    "SAMPLE_ID",
    "EXPERIMENTID",
    "SAMPLESOURCE",
    "SAMPLECLASS",
    "BIOLOGICAL_SEX",
    "SAMPLE_TYPE",
    "LIB_PREP",
    "SEQ_VENDOR",
    "SEQ_PLATFORM",
    "LANE",
    "SEQBC_ID",
    "PATH_TO_CONCORDANCE_DATA_DIR",
    "CONCORDANCE_CONTROL_PATH",
    "TRUTH_DATA_DIR",
    "ILMN_R1_FQ",
    "ILMN_R2_FQ",
    "STAGE_DIRECTIVE",
    "STAGE_TARGET",
    "SUBSAMPLE_PCT",
    "SAMPLEUSE",
    "BWA_KMER",
    "DEEP_MODEL",
    "IS_POS_CTRL",
    "IS_NEG_CTRL",
    "TUM_NRM_SAMPLEID_MATCH",
    "N_X",
    "N_Y",
    "EXTERNAL_SAMPLE_ID",
]

FASTQ_RE = re.compile(
    r"^(?P<sample>.+)_S(?P<sample_no>[0-9]+)_L(?P<lane>[0-9]{3})_"
    r"R(?P<read>[12])_001\.fastq\.gz$"
)


def run_aws(args: list[str]) -> str:
    env = dict(os.environ)
    env.setdefault("AWS_PROFILE", "daylily-service-lsmc")
    env.setdefault("AWS_REGION", "us-west-2")
    env.setdefault("AWS_DEFAULT_REGION", "us-west-2")
    proc = subprocess.run(
        ["aws", *args],
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise SystemExit(proc.returncode)
    return proc.stdout


def sample_sheet_ids() -> list[str]:
    text = run_aws(["s3", "cp", SAMPLE_SHEET, "-"])
    rows: list[str] = []
    in_data = False
    reader = csv.reader(text.splitlines())
    for row in reader:
        if not row:
            continue
        if in_data and row[0].startswith("["):
            break
        if row[0] == "[BCLConvert_Data]":
            in_data = True
            continue
        if not in_data or row[0] == "Sample_ID":
            continue
        rows.append(row[0])
    seen: set[str] = set()
    unique: list[str] = []
    for sample_id in rows:
        if sample_id not in seen:
            unique.append(sample_id)
            seen.add(sample_id)
    return unique


def fastq_pairs() -> dict[tuple[str, str], dict[str, str]]:
    listing = run_aws(["s3", "ls", FASTQ_PREFIX, "--recursive"])
    pairs: dict[tuple[str, str], dict[str, str]] = defaultdict(dict)
    total_fastq = 0
    for line in listing.splitlines():
        fields = line.split(maxsplit=3)
        if len(fields) != 4:
            continue
        key = fields[3]
        name = key.rsplit("/", 1)[-1]
        if not name.endswith(".fastq.gz"):
            continue
        total_fastq += 1
        match = FASTQ_RE.match(name)
        if not match:
            raise SystemExit(f"Could not parse FASTQ name: {name}")
        sample = match.group("sample")
        lane = match.group("lane")
        read = match.group("read")
        pairs[(sample, lane)][read] = f"s3://lsmc-ssf-sequencing-data/{key}"

    missing = sorted((sample, lane, reads) for (sample, lane), reads in pairs.items() if set(reads) != {"1", "2"})
    if missing:
        raise SystemExit(f"Missing mates: {missing[:10]}")

    undetermined = sum(1 for sample, _lane in pairs if sample == "Undetermined")
    named = len(pairs) - undetermined
    print(f"observed_fastq_objects={total_fastq}")
    print(f"observed_pairs={len(pairs)}")
    print(f"observed_named_pairs={named}")
    print(f"observed_undetermined_pairs={undetermined}")
    return pairs


def row_for(sample_id: str, lane: str, reads: dict[str, str]) -> dict[str, str]:
    is_hg003 = sample_id == "LSMC_HG003"
    is_ntc = sample_id == "LSMC_NTC"
    concordance = GIAB_HG003 if is_hg003 else "na"
    return {
        "RUN_ID": RUN_ID,
        "SAMPLE_ID": sample_id,
        "EXPERIMENTID": EXPERIMENT_ID,
        "SAMPLESOURCE": "blood",
        "SAMPLECLASS": "research",
        "BIOLOGICAL_SEX": "male" if is_hg003 else "na",
        "SAMPLE_TYPE": "gdna",
        "LIB_PREP": "NOAMPWGS",
        "SEQ_VENDOR": "ILMN",
        "SEQ_PLATFORM": "NOVASEQX",
        "LANE": str(int(lane)),
        "SEQBC_ID": "D0",
        "PATH_TO_CONCORDANCE_DATA_DIR": concordance,
        "CONCORDANCE_CONTROL_PATH": concordance,
        "TRUTH_DATA_DIR": concordance,
        "ILMN_R1_FQ": reads["1"],
        "ILMN_R2_FQ": reads["2"],
        "STAGE_DIRECTIVE": "stage_data",
        "STAGE_TARGET": "/data/staged_sample_data",
        "SUBSAMPLE_PCT": "na",
        "SAMPLEUSE": "posControl" if is_hg003 else ("negControl" if is_ntc else "sample"),
        "BWA_KMER": "19",
        "DEEP_MODEL": "WGS",
        "IS_POS_CTRL": "true" if is_hg003 else "false",
        "IS_NEG_CTRL": "true" if is_ntc else "false",
        "TUM_NRM_SAMPLEID_MATCH": "na",
        "N_X": "1" if is_hg003 else "na",
        "N_Y": "1" if is_hg003 else "na",
        "EXTERNAL_SAMPLE_ID": "HG003" if is_hg003 else sample_id,
    }


def main() -> int:
    out_path = Path(sys.argv[1]).resolve()
    samples = sample_sheet_ids()
    pairs = fastq_pairs()
    named_samples = [sample for sample in samples if sample != "Undetermined"]

    if len(named_samples) != 96:
        raise SystemExit(f"Expected 96 named samples, saw {len(named_samples)}")
    if len(pairs) != 776:
        raise SystemExit(f"Expected 776 total pairs, saw {len(pairs)}")
    if sum(1 for sample, _lane in pairs if sample != "Undetermined") != 768:
        raise SystemExit("Expected 768 named pairs")
    if sum(1 for sample, _lane in pairs if sample == "Undetermined") != 8:
        raise SystemExit("Expected 8 Undetermined pairs")

    missing_samples = sorted(set(named_samples) - {sample for sample, _lane in pairs})
    if missing_samples:
        raise SystemExit(f"SampleSheet samples missing FASTQs: {missing_samples}")

    rows: list[dict[str, str]] = []
    for sample_id in named_samples:
        lanes = sorted(lane for sample, lane in pairs if sample == sample_id)
        if lanes != [f"{idx:03d}" for idx in range(1, 9)]:
            raise SystemExit(f"{sample_id} lanes were {lanes}, expected L001-L008")
        for lane in lanes:
            rows.append(row_for(sample_id, lane, pairs[(sample_id, lane)]))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"sample_sheet_named_samples={len(named_samples)}")
    print(f"manifest_rows={len(rows)}")
    print(f"manifest_path={out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
