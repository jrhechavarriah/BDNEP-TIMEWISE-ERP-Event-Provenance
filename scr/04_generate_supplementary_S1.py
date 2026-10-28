#!/usr/bin/env python3
"""Rebuild Supplementary Figure S1 from the verified P14 independent replay.

Scientific inputs: P14 `recomputed_participants.csv` and `replay_verification_report.json`.
No EEG reprocessing or artifact rejection is performed. Output is an editorial
visualization of behaviorally eligible trials, not a new scientific analysis.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
from bdnep_inputs import verified_artifact

EXPECTED_IDS = [f"sub-{i:03d}" for i in range(4, 41)]
P14_PASS = "FULL_INDEPENDENT_REPLAY_MATCHES_FROZEN_NUMERICAL_REFERENCES"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def check_p14_file(replay_dir: Path, name: str) -> Path:
    return verified_artifact(replay_dir, name)

def collect_trials(replay_dir: Path):
    report_path = check_p14_file(replay_dir, "replay_verification_report.json")
    rows_path = check_p14_file(replay_dir, "recomputed_participants.csv")
    report = json.loads(report_path.read_text(encoding="utf-8-sig"))
    if report.get("status") != P14_PASS or report.get("n_subjects") != 37:
        raise RuntimeError("P14 full replay status not verified")
    with rows_path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        required = {"subject_id", "n_A", "n_B", "n_B_face", "n_B_car",
                    "n_B_scrambled_face", "n_B_scrambled_car"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise RuntimeError("P14 participant table is missing required columns")
        rows = list(reader)
    if [r["subject_id"] for r in rows] != EXPECTED_IDS:
        raise RuntimeError("P14 participant order or cohort mismatch")
    y = []
    for r in rows:
        n_a, n_b = int(r["n_A"]), int(r["n_B"])
        condition_sum = sum(int(r[c]) for c in (
            "n_B_face", "n_B_car", "n_B_scrambled_face", "n_B_scrambled_car"))
        if n_a != 320 or not (0 <= n_b <= 320) or condition_sum != n_b:
            raise RuntimeError(f"Invalid trial accounting: {r['subject_id']}")
        y.append(n_b)
    if sum(y) != 10810 or sum(int(r["n_A"]) for r in rows) != 11840:
        raise RuntimeError("Arm B or Arm A aggregate count diverged from frozen evidence")
    return EXPECTED_IDS, y, {"replay_report_sha256": sha256(report_path),
                             "participant_csv_sha256": sha256(rows_path),
                             "arm_a_stimuli": 11840,
                             "arm_b_eligible": 10810,
                             "arm_b_ineligible": 1030}


def draw_s1(ids, eligible, output: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    x = np.arange(len(ids))
    fig, ax = plt.subplots(figsize=(12.6, 5.2))
    ax.bar(x, eligible, width=0.75)
    ax.axhline(320, linestyle="--", linewidth=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels([s.removeprefix("sub-") for s in ids], rotation=90, fontsize=8)
    ax.set_xlabel("Participant")
    ax.set_ylabel("Behaviorally eligible N170 trials")
    ax.set_title("Arm B behavioral trial availability")
    ax.set_ylim(0, 340)
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    inputs = ap.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--repo", dest="replay_dir", type=Path)
    inputs.add_argument("--replay-dir", dest="replay_dir", type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()
    ids, eligible, provenance = collect_trials(args.replay_dir)
    if args.out_dir.exists():
        raise RuntimeError("Refusing to overwrite output folder: " + str(args.out_dir))
    args.out_dir.mkdir(parents=True, exist_ok=False)
    filename = "Supplementary_Figure_S1_Behavioral_Trial_Availability.png"
    target = args.out_dir / filename
    draw_s1(ids, eligible, target)
    receipt = {"status": "S1_REGENERATED_FROM_VERIFIED_P14",
               "figure": filename, "figure_sha256": sha256(target),
               "n_subjects": 37, "event_selection": "Arm B correct response, RT 200-1000 ms",
               "artifact_rejection": "NOT_PERFORMED", **provenance,
               }
    (args.out_dir / "Supplementary_S1_provenance.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print("STATUS: S1_REGENERATED_FROM_VERIFIED_P14")
    print("ELIGIBLE_ARM_B: 10810/11840")
    print("TITLE: Arm B behavioral trial availability")
    print("FIGURE:", target)
    print("SHA256:", receipt["figure_sha256"])


if __name__ == "__main__":
    main()
