#!/usr/bin/env python3
"""Launch the 23andMe per-lane smoke run inside a remote tmux session."""

from __future__ import annotations

import argparse
import shlex
from pathlib import Path

from daylily_ec.aws.ssm import run_shell, write_remote_text


SESSION = "23ame-ilmn-1-bylane"
DESTINATION = "23ame-ilmn-1-bylane"
SMOKE_STAGE_FSX = (
    "/fsx/data/staged_sample_data/"
    "remote_stage_20260511T084859Z_23ame_ilmn_1_bylane_smoke3"
)
FULL_STAGE_FSX = (
    "/fsx/data/staged_sample_data/"
    "remote_stage_20260511T084859Z_23ame_ilmn_1_bylane_full"
)
RUN_DIR = f"/home/ubuntu/daylily-runs/{SESSION}"
REMOTE_SCRIPT = f"{RUN_DIR}/run_smoke_controller.sh"
REPO_PATH = f"/fsx/analysis_results/ubuntu/{DESTINATION}/daylily-omics-analysis"


def remote_controller(run_mode: str) -> str:
    if run_mode == "smoke":
        stage_fsx = SMOKE_STAGE_FSX
        sample_manifest = "smoke3_samples.tsv"
        unit_manifest = "smoke3_units.tsv"
        status_prefix = ""
        expected_samples = 3
        expected_units = 24
        exact_rule_count_check = False
    elif run_mode == "full":
        stage_fsx = FULL_STAGE_FSX
        sample_manifest = "full_samples.tsv"
        unit_manifest = "full_units.tsv"
        status_prefix = "full_"
        expected_samples = 96
        expected_units = 768
        exact_rule_count_check = False
    else:
        raise ValueError(f"unsupported run mode: {run_mode}")

    edit_yaml = """from pathlib import Path
p = Path('config/day_profiles/slurm/rule_config.yaml')
lines = p.read_text().splitlines()
out = []
in_multiqc = False
seen_include = False
seen_enable = False
seen_disable = False
disabled_tools = ['expansionhunter', 'peddy', 'site_mix', 'verifybamid2']
disabled_line = '    disable_tools: ' + repr(disabled_tools)
for idx, line in enumerate(lines):
    if line == 'multiqc_qc:':
        in_multiqc = True
        out.append(line)
        continue
    if in_multiqc and line and not line.startswith(' '):
        if not seen_include:
            out.append('    include_no_dedup_alignment_qc: false')
        if not seen_enable:
            out.append(\"    enable_tools: ['vep']\")
        if not seen_disable:
            out.append(disabled_line)
        in_multiqc = False
    if in_multiqc and line.strip().startswith('include_no_dedup_alignment_qc:'):
        out.append('    include_no_dedup_alignment_qc: false')
        seen_include = True
        continue
    if in_multiqc and line.strip().startswith('enable_tools:'):
        out.append(\"    enable_tools: ['vep']\")
        seen_enable = True
        continue
    if in_multiqc and line.strip().startswith('disable_tools:'):
        out.append(disabled_line)
        seen_disable = True
        continue
    out.append(line)
if in_multiqc:
    if not seen_include:
        out.append('    include_no_dedup_alignment_qc: false')
    if not seen_enable:
        out.append(\"    enable_tools: ['vep']\")
    if not seen_disable:
        out.append(disabled_line)
text = '\\n'.join(out) + '\\n'
if 'include_no_dedup_alignment_qc: false' not in text:
    raise SystemExit('multiqc_qc.include_no_dedup_alignment_qc was not written')
if \"enable_tools: ['vep']\" not in text:
    raise SystemExit('multiqc_qc.enable_tools vep was not written')
if disabled_line.strip() not in text:
    raise SystemExit('multiqc_qc.disable_tools was not written')
p.write_text('\\n'.join(out) + '\\n')
"""
    patch_multiqc_intro = """from pathlib import Path
p = Path('workflow/scripts/build_multiqc_intro.py')
text = p.read_text()
old = \"\"\"    if not raw:
        raise ValueError(f\\\"Benchmark row {line_number} has empty task_cost: {path}\\\")
\"\"\"
new = \"\"\"    if not raw or raw.upper() in {\\\"NA\\\", \\\"N/A\\\", \\\"NONE\\\"}:
        return Decimal(\\\"0\\\")
\"\"\"
if old in text:
    p.write_text(text.replace(old, new))
elif new not in text:
    raise SystemExit('build_multiqc_intro.py task_cost patch did not match')
"""
    validate_manifests = f"""import csv
from pathlib import Path
samples = list(csv.DictReader(Path('config/samples.tsv').open(), delimiter='\\t'))
units = list(csv.DictReader(Path('config/units.tsv').open(), delimiter='\\t'))
expected_samples = {expected_samples}
expected_units = {expected_units}
if len(samples) != expected_samples:
    raise SystemExit(f'samples.tsv rows={{len(samples)}} expected={{expected_samples}}')
if len(units) != expected_units:
    raise SystemExit(f'units.tsv rows={{len(units)}} expected={{expected_units}}')
if any(row.get('LANEID') == '0' for row in units):
    raise SystemExit('units.tsv contains LANEID=0')
if any('Undetermined' in row.get('SAMPLEID', '') for row in units):
    raise SystemExit('units.tsv contains Undetermined')
"""
    return f"""#!/usr/bin/env bash
set -euo pipefail
SESSION={shlex.quote(SESSION)}
DESTINATION={shlex.quote(DESTINATION)}
RUN_MODE={shlex.quote(run_mode)}
RUN_DIR={shlex.quote(RUN_DIR)}
REPO_PATH={shlex.quote(REPO_PATH)}
STAGE_FSX={shlex.quote(stage_fsx)}
SAMPLE_MANIFEST={shlex.quote(sample_manifest)}
UNIT_MANIFEST={shlex.quote(unit_manifest)}
EXPECTED_SAMPLES={expected_samples}
EXPECTED_UNITS={expected_units}
EXACT_RULE_COUNT_CHECK={str(exact_rule_count_check).lower()}
STATUS_FILE="$RUN_DIR/{status_prefix}status.env"
DRYRUN_LOG="$RUN_DIR/{status_prefix}dryrun.log"
PROD_LOG="$RUN_DIR/{status_prefix}production.log"
COUNTS_FILE="$RUN_DIR/{status_prefix}dryrun_counts.env"
FORBIDDEN_FILE="$RUN_DIR/{status_prefix}dryrun_forbidden_hits.txt"
DRYRUN_EXIT_FILE="$RUN_DIR/{status_prefix}dryrun_exit.env"
PROD_EXIT_FILE="$RUN_DIR/{status_prefix}production_exit.env"
mkdir -p "$RUN_DIR"
write_status() {{
  printf 'stage=%s\\nrun_mode=%s\\nupdated_at=%s\\nrepo_path=%s\\n' "$1" "$RUN_MODE" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$REPO_PATH" > "$STATUS_FILE"
}}
keep_shell() {{
  write_status "$1"
  exec bash -il
}}
write_status starting
if [[ "$(id -un)" != "ubuntu" ]]; then
  echo "ERROR: expected ubuntu user, got $(id -un)"
  keep_shell wrong_user
fi
if [[ ! -f "$STAGE_FSX/manifests/$SAMPLE_MANIFEST" || ! -f "$STAGE_FSX/manifests/$UNIT_MANIFEST" ]]; then
  echo "ERROR: manifests missing under $STAGE_FSX/manifests"
  keep_shell missing_manifests
fi
write_status cloning
if [[ -e "$REPO_PATH" ]]; then
  echo "INFO: reusing existing cloned repo path: $REPO_PATH"
else
  day-clone -t main -d "$DESTINATION"
fi
cd "$REPO_PATH"
mkdir -p config
cp "$STAGE_FSX/manifests/$SAMPLE_MANIFEST" config/samples.tsv
cp "$STAGE_FSX/manifests/$UNIT_MANIFEST" config/units.tsv
python -c {shlex.quote(validate_manifests)}
write_status activating
shopt -s expand_aliases
if [[ ! -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
  echo "ERROR: missing conda profile script"
  keep_shell missing_conda_profile
fi
source "$HOME/miniconda3/etc/profile.d/conda.sh"
set +e
set +u
source dyoainit --project da-us-west-2d-fk-260509-use
dyoainit_rc=$?
if [[ "$dyoainit_rc" != "0" ]]; then
  set -u
  set -e
  echo "ERROR: dyoainit failed with status $dyoainit_rc"
  keep_shell dyoainit_failed
fi
dy-a slurm hg38
dya_rc=$?
set -u
set -e
if [[ "$dya_rc" != "0" ]]; then
  echo "ERROR: dy-a failed with status $dya_rc"
  keep_shell dya_failed
fi
python -c {shlex.quote(edit_yaml)}
python -c {shlex.quote(patch_multiqc_intro)}
grep -A6 '^multiqc_qc:' config/day_profiles/slurm/rule_config.yaml | tee "$RUN_DIR/active_multiqc_qc.txt"
if ! grep -q 'include_no_dedup_alignment_qc: false' config/day_profiles/slurm/rule_config.yaml; then
  echo "ERROR: active slurm profile did not disable no-dedup alignment QC"
  keep_shell profile_edit_failed
fi
CMD=(
  bin/day_run
  produce_snv_concordances
  produce_alignstats
  produce_tiddit
  produce_vep
  produce_relatedness
  produce_multiqc_final
  -j 120
  -k
  -p
  --config
  "aligners=['sent']"
  "dedupers=['dppl']"
  "snv_callers=['sentd']"
  "sv_callers=['tiddit']"
)
printf '%q ' "${{CMD[@]}}" > "$RUN_DIR/{status_prefix}production.command"
printf '%q ' "${{CMD[@]}}" > "$RUN_DIR/{status_prefix}dryrun.command"
printf ' -n\\n' >> "$RUN_DIR/{status_prefix}dryrun.command"
write_status dryrun_running
set +e
"${{CMD[@]}}" -n > "$DRYRUN_LOG" 2>&1
dryrun_rc=$?
set -e
echo "DRYRUN_EXIT=$dryrun_rc" | tee "$DRYRUN_EXIT_FILE"
if [[ "$dryrun_rc" != "0" ]]; then
  tail -200 "$DRYRUN_LOG"
  keep_shell dryrun_failed
fi
sentieon_count=$(grep -c '^rule sentieon_bwa_sort:' "$DRYRUN_LOG" || true)
doppel_count=$(grep -c '^rule doppelmark_dups:' "$DRYRUN_LOG" || true)
tiddit_count=$(grep -c '^rule tiddit:' "$DRYRUN_LOG" || true)
tiddit_sort_count=$(grep -c '^rule tiddit_sort_index:' "$DRYRUN_LOG" || true)
vep_rule_count=$(grep -Ec '^rule .*vep' "$DRYRUN_LOG" || true)
expansionhunter_rule_count=$(grep -Ec '^rule .*expansionhunter' "$DRYRUN_LOG" || true)
verifybamid2_rule_count=$(grep -Ec '^rule .*verifybamid2' "$DRYRUN_LOG" || true)
site_mix_rule_count=$(grep -Ec '^rule .*site_mix' "$DRYRUN_LOG" || true)
peddy_rule_count=$(grep -Ec '^rule .*peddy' "$DRYRUN_LOG" || true)
relatedness_extract_count=$(grep -c '^rule relatedness_batch_somalier_extract:' "$DRYRUN_LOG" || true)
relatedness_relate_count=$(grep -c '^rule relatedness_batch_somalier_relate:' "$DRYRUN_LOG" || true)
relatedness_report_count=$(grep -c '^rule relatedness_batch_report:' "$DRYRUN_LOG" || true)
relatedness_gather_count=$(grep -c '^rule relatedness_batch_gather:' "$DRYRUN_LOG" || true)
multiqc_final_wgs_count=$(grep -c '^rule multiqc_final_wgs:' "$DRYRUN_LOG" || true)
printf 'sentieon_bwa_sort=%s\\ndoppelmark_dups=%s\\ntiddit=%s\\ntiddit_sort_index=%s\\nvep_rules=%s\\nrelatedness_batch_somalier_extract=%s\\nrelatedness_batch_somalier_relate=%s\\nrelatedness_batch_report=%s\\nrelatedness_batch_gather=%s\\nmultiqc_final_wgs=%s\\nexpansionhunter_rules=%s\\nverifybamid2_rules=%s\\nsite_mix_rules=%s\\npeddy_rules=%s\\n' \\
  "$sentieon_count" "$doppel_count" "$tiddit_count" "$tiddit_sort_count" "$vep_rule_count" "$relatedness_extract_count" "$relatedness_relate_count" "$relatedness_report_count" "$relatedness_gather_count" "$multiqc_final_wgs_count" "$expansionhunter_rule_count" "$verifybamid2_rule_count" "$site_mix_rule_count" "$peddy_rule_count" | tee "$COUNTS_FILE"
grep -E '/na/|[.]na[.]|duphold|manta|dysgu' "$DRYRUN_LOG" > "$FORBIDDEN_FILE" || true
if [[ -s "$FORBIDDEN_FILE" ]]; then
  echo "ERROR: forbidden dry-run hits found"
  cat "$FORBIDDEN_FILE"
  keep_shell dryrun_forbidden_hits
fi
if [[ "$EXACT_RULE_COUNT_CHECK" == "true" && ( "$sentieon_count" != "24" || "$doppel_count" != "24" || "$tiddit_count" != "24" ) ]]; then
  echo "ERROR: unexpected dry-run counts"
  cat "$COUNTS_FILE"
  keep_shell dryrun_unexpected_counts
fi
if [[ "$vep_rule_count" == "0" ]]; then
  echo "ERROR: no VEP dry-run rules found"
  cat "$COUNTS_FILE"
  keep_shell dryrun_missing_vep
fi
if [[ "$expansionhunter_rule_count" != "0" ]]; then
  echo "ERROR: ExpansionHunter dry-run rules were present despite disabled sample-sex blocker"
  cat "$COUNTS_FILE"
  keep_shell dryrun_unexpected_expansionhunter
fi
if [[ "$verifybamid2_rule_count" != "0" || "$site_mix_rule_count" != "0" || "$peddy_rule_count" != "0" ]]; then
  echo "ERROR: disabled optional QC rules were present in dry-run"
  cat "$COUNTS_FILE"
  keep_shell dryrun_unexpected_optional_qc
fi
if [[ ! -f results/day/hg38/other_reports/relatedness_mqc.tsv ]]; then
  if [[ "$relatedness_extract_count" != "$EXPECTED_UNITS" || "$relatedness_relate_count" != "1" || "$relatedness_report_count" != "1" || "$relatedness_gather_count" != "1" ]]; then
    echo "ERROR: unexpected relatedness dry-run counts"
    cat "$COUNTS_FILE"
    keep_shell dryrun_unexpected_relatedness
  fi
fi
write_status production_running
set +e
"${{CMD[@]}}" > "$PROD_LOG" 2>&1
prod_rc=$?
set -e
echo "PRODUCTION_EXIT=$prod_rc" | tee "$PROD_EXIT_FILE"
if [[ "$prod_rc" == "0" ]]; then
  keep_shell production_completed
fi
tail -200 "$PROD_LOG"
keep_shell production_failed
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--instance-id", required=True)
    parser.add_argument("--region", default="us-west-2")
    parser.add_argument("--profile", default="daylily-service-lsmc")
    parser.add_argument("--run-mode", choices=("smoke", "full"), default="smoke")
    args = parser.parse_args()

    write_remote_text(
        args.instance_id,
        args.region,
        REMOTE_SCRIPT,
        remote_controller(args.run_mode),
        profile=args.profile,
    )
    launch_script = "\n".join(
        [
            "set -euo pipefail",
            f"SESSION={shlex.quote(SESSION)}",
            f"RUN_DIR={shlex.quote(RUN_DIR)}",
            f"REMOTE_SCRIPT={shlex.quote(REMOTE_SCRIPT)}",
            'mkdir -p "$RUN_DIR"',
            'chmod 700 "$REMOTE_SCRIPT"',
            'if tmux has-session -t "$SESSION" 2>/dev/null; then',
            '  windows=$(tmux list-windows -t "$SESSION" | wc -l | tr -d " ")',
            '  panes=$(tmux list-panes -t "$SESSION" | wc -l | tr -d " ")',
            '  if [[ "$windows" != "1" || "$panes" != "1" ]]; then echo "ERROR: tmux session has windows=$windows panes=$panes"; exit 9; fi',
            '  tmux send-keys -t "$SESSION" "bash -ilc \'$REMOTE_SCRIPT >>$RUN_DIR/tmux.log 2>&1\'" C-m',
            'else',
            '  tmux new-session -d -s "$SESSION" "bash -ilc \'$REMOTE_SCRIPT >>$RUN_DIR/tmux.log 2>&1\'"',
            'fi',
            'sleep 2',
            'tmux has-session -t "$SESSION"',
            'printf "session=%s\\nrun_dir=%s\\nremote_script=%s\\n" "$SESSION" "$RUN_DIR" "$REMOTE_SCRIPT"',
        ]
    )
    result = run_shell(
        args.instance_id,
        args.region,
        launch_script,
        profile=args.profile,
        timeout=120,
        comment="Launch 23andMe per-lane smoke tmux",
    )
    print(result.stdout, end="")
    print(result.stderr, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
