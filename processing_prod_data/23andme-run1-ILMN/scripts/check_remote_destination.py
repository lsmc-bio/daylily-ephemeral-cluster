#!/usr/bin/env python3
"""Fail if the target DayOA destination/session already exists on the headnode."""

from __future__ import annotations

import sys

from daylily_ec.aws.ssm import resolve_headnode_instance_id, run_shell, wait_for_ssm_online


PROFILE = "daylily-service-lsmc"
REGION = "us-west-2"
CLUSTER = "fk-260509-use"
DESTINATION = "23andme-run1-ILMN"


SCRIPT = f"""
set -euo pipefail
if [[ "$(id -un)" != "ubuntu" ]]; then
  echo "__DAYLILY_ERROR__=wrong_user"
  exit 64
fi
destination={DESTINATION!r}
analysis_path="/fsx/analysis_results/ubuntu/${{destination}}"
run_dir="/home/ubuntu/daylily-runs/${{destination}}"
if [[ -e "$analysis_path" ]]; then
  echo "__DAYLILY_ERROR__=destination_exists"
  echo "$analysis_path"
  exit 20
fi
if [[ -e "$run_dir" ]]; then
  echo "__DAYLILY_ERROR__=run_dir_exists"
  echo "$run_dir"
  exit 21
fi
if command -v tmux >/dev/null 2>&1 && tmux has-session -t "$destination" 2>/dev/null; then
  echo "__DAYLILY_ERROR__=tmux_session_exists"
  echo "$destination"
  exit 22
fi
echo "__DAYLILY_DESTINATION_OK__=$destination"
"""


def main() -> int:
    target = resolve_headnode_instance_id(CLUSTER, REGION, profile=PROFILE)
    wait_for_ssm_online(target.instance_id, REGION, profile=PROFILE, timeout=120)
    result = run_shell(
        target.instance_id,
        REGION,
        SCRIPT,
        profile=PROFILE,
        timeout=120,
        comment="Check 23andMe ILMN DayOA destination",
    )
    if result.stdout:
        print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
