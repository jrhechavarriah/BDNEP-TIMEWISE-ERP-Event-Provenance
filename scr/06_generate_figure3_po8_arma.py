#!/usr/bin/env python3
"""Generate publication Figure 3 from a checksum-verified P14 FULL replay.

Source: P14 recomputed_participants.csv, not archived/frozen original results.
No EEG reprocessing, model fitting, OSF downloads or publication actions.
Only Arm A (prespecified primary outcome) is plotted; Arm B is descriptive
and intentionally absent. This is an internal-staging, IP-Gate-HELD asset.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from bdnep_inputs import verified_artifact

EXPECTED_CONTRACT_SHA256 = (
    "f75505d2b30f57d4b6ce88a5efcd60c2cf27398e7021dbe2545f34dab4d4827c"
)
EXPECTED_IDS = [f"sub-{n:03d}" for n in range(4, 41)]
EXPECTED_MEAN_UV = -2.136471428450268
MEAN_TOL_UV = 1e-6


def require(ok: bool, why: str) -> None:
    if not ok:
        raise RuntimeError(why)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verified_source(replay: Path, relative: str, sums: dict[str, str]) -> Path:
    return verified_artifact(replay, relative)

def parse_sums(path: Path) -> dict[str, str]:
    require(path.is_file(), "Missing P14 SHA256SUMS.txt")
    result = {}
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        if not line.strip():
            continue
        bits = line.split(maxsplit=1)
        require(len(bits) == 2, "Malformed SHA256SUMS line")
        h, name = bits
        name = name.strip().lstrip("*").replace("\\", "/")
        require(len(h) == 64 and all(c in "0123456789abcdefABCDEF" for c in h),
                "Malformed SHA256SUMS digest")
        require(name not in result and not name.startswith("/") and ".." not in Path(name).parts,
                "Duplicate/unsafe SHA256SUMS entry")
        result[name] = h.lower()
    return result


def load_and_verify(replay: Path):
    sums = parse_sums(replay / "SHA256SUMS.txt")
    report_path = verified_source(replay, "replay_verification_report.json", sums)
    participants_path = verified_source(replay, "recomputed_participants.csv", sums)
    report = json.loads(report_path.read_text(encoding="utf-8-sig"))
    require(report.get("status") == "FULL_INDEPENDENT_REPLAY_MATCHES_FROZEN_NUMERICAL_REFERENCES",
            "P14 FULL replay not verified")
    require(report.get("source_contract_sha256") == EXPECTED_CONTRACT_SHA256,
            "Unexpected P14 source contract")
    require(report.get("subjects") == EXPECTED_IDS and report.get("n_subjects") == 37,
            "Unexpected P14 participants")
    require(report.get("group", {}).get("n_A_stimuli") == 11840,
            "Unexpected Arm A event count")
    with participants_path.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        needed = {
            "subject_id", "n_A", "PO8_A_original_uV", "PO8_A_corrected_uV",
            "PO8_A_delta_uV", "PO8_A_sign_flip",
        }
        require(needed.issubset(set(reader.fieldnames or [])), "P14 columns incomplete")
        rows = list(reader)
    ids = [r["subject_id"] for r in rows]
    require(ids == EXPECTED_IDS, "P14 CSV participant ordering or cohort mismatch")
    values = []
    flips = 0
    for r in rows:
        require(int(r["n_A"]) == 320, "Unexpected Arm A events: " + r["subject_id"])
        original = float(r["PO8_A_original_uV"])
        corrected = float(r["PO8_A_corrected_uV"])
        delta = float(r["PO8_A_delta_uV"])
        require(all(math.isfinite(x) for x in (original, corrected, delta)),
                "Nonfinite ERP measure: " + r["subject_id"])
        require(abs((corrected - original) - delta) < 1e-10,
                "Unreconciled participant difference: " + r["subject_id"])
        sign_flip = (original * corrected < 0)
        require(r["PO8_A_sign_flip"].strip().lower() == str(sign_flip).lower(),
                "Sign-flip flag disagreement: " + r["subject_id"])
        flips += sign_flip
        values.append(delta)
    mean = math.fsum(values) / 37
    require(abs(mean - EXPECTED_MEAN_UV) <= MEAN_TOL_UV,
            "Primary Arm A mean differs from frozen reference")
    require(abs(mean - float(report["group"]["PO8_A_mean_delta_uV"])) <= MEAN_TOL_UV,
            "P14 report mean differs from recalculated mean")
    require(flips == 7 and flips == report["group"]["PO8_A_sign_flips"],
            "Unexpected primary sign reversals")
    require(sum(v < 0 for v in values) == 31 and sum(v > 0 for v in values) == 6,
            "Unexpected directions of paired differences")
    return report_path, participants_path, values, mean, flips


def generate(replay: Path, out: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    report_path, participants_path, values, mean, flips = load_and_verify(replay)
    require(not out.exists(), "Output folder already exists: " + str(out))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=False, exist_ok=False)

    fig, ax = plt.subplots(figsize=(10.6, 5.0))
    x = list(range(37))
    ax.scatter(x, values, s=29, zorder=4, label="Participant-level difference")
    ax.axhline(0, linewidth=0.9, color="0.40", label="Zero difference")
    ax.axhline(mean, linewidth=1.1, linestyle="--", color="0.25",
               label=f"Arm A group mean ({mean:.3f} µV)")
    ax.set_xticks(x)
    ax.set_xticklabels(EXPECTED_IDS, rotation=90, fontsize=7.3)
    ax.set_xlabel("Participant")
    ax.set_ylabel("Corrected minus original PO8 face-minus-car contrast (µV)")
    ax.set_title("Participant-level PO8 Arm A digital-reference sensitivity")
    ax.grid(axis="y", alpha=0.14)
    ax.legend(loc="lower left", frameon=False, fontsize=8)
    fig.tight_layout()

    png = out / "Figure_3_PO8_ArmA_Participant_Deltas.png"
    svg = out / "Figure_3_PO8_ArmA_Participant_Deltas.svg"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)

    metadata = {
        "stage": "P19_FIGURE3_FROM_VERIFIED_P14_REPLAY",
        "status": "FIGURE3_PO8_ARMA_GENERATED_VERIFIED",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_replay_status": "FULL_INDEPENDENT_REPLAY_MATCHES_FROZEN_NUMERICAL_REFERENCES",
        "source_contract_sha256": EXPECTED_CONTRACT_SHA256,
        "recomputed_participants_sha256": sha256(participants_path),
        "replay_verification_report_sha256": sha256(report_path),
        "n_participants": 37,
        "displayed_analysis": "PO8 Arm A, primary only",
        "arm_b_displayed": False,
        "po8_arm_a_mean_delta_uV": mean,
        "n_negative_differences": 31,
        "n_positive_differences": 6,
        "face_minus_car_sign_flips": flips,
        "no_eeg_reprocessing": True,
    }
    meta_path = out / "figure3_generation_report.json"
    meta_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    sums_path = out / "SHA256SUMS.txt"
    with sums_path.open("x", encoding="utf-8") as handle:
        for path in (png, svg, meta_path):
            handle.write(f"{sha256(path)}  {path.name}\n")
    print("STATUS: FIGURE3_PO8_ARMA_GENERATED_VERIFIED")
    print("PARTICIPANTS: 37/37")
    print("PO8_ARMA_GROUP_MEAN_uV:", f"{mean:.15f}")
    print("PAIRED_DIFFERENCES: 31 negative, 6 positive")
    print("FACE_MINUS_CAR_SIGN_REVERSALS:", flips)
    print("PNG:", png)
    print("SVG:", svg)
    print("SHA256:", sums_path)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sources = p.add_mutually_exclusive_group(required=True)
    sources.add_argument("--repo", dest="replay_dir", type=Path)
    sources.add_argument("--replay-dir", dest="replay_dir", type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    args = p.parse_args()
    generate(args.replay_dir.resolve(), args.out_dir.resolve())


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[FAIL] FIGURE3: {type(exc).__name__}: {exc}", file=sys.stderr)
        sys.exit(2)
