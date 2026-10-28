#!/usr/bin/env python3
"""Generate Supplementary Figure S2 from checksum-verified BDNEP P14 replay counts.

Shows the two conditions that form the face-minus-car ERP contrast under Arm B;
this is a descriptive trial-availability visualization, not an alternative analysis.

Usage:
  python scr/05_generate_supplementary_S2_face_car.py --repo <BDNEP path> --check
  python scr/05_generate_supplementary_S2_face_car.py --repo <BDNEP path> --run

The --run mode verifies existing SHA256SUMS entries before creating new files,
refuses to overwrite output, and then atomically updates the root manifest.
No EEG is loaded, processed, modified or released.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np

EXPECTED_SUBJECTS = [f"sub-{i:03d}" for i in range(4, 41)]
P14_SUCCESS = "FULL_INDEPENDENT_REPLAY_MATCHES_FROZEN_NUMERICAL_REFERENCES"
EXPECTED_COUNTS = {
    "n_B_face": 2720,
    "n_B_car": 2691,
    "n_B_scrambled_face": 2685,
    "n_B_scrambled_car": 2714,
}
SCRIPT_NAME = "05_generate_supplementary_S2_face_car.py"
FIGURE_NAME = "Supplementary_Figure_S2_Face_Car_Behavioral_Eligibility"


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_checksums(repo: Path) -> dict[str, str]:
    path = repo / "SHA256SUMS.txt"
    require(path.is_file() and not path.is_symlink(), "Missing root SHA256SUMS.txt")
    records: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        if not line.strip():
            continue
        parts = line.split(maxsplit=1)
        require(len(parts) == 2, "Malformed SHA256SUMS line")
        value, rel = parts
        rel = rel.strip().lstrip("*").replace("\\", "/")
        require(re.fullmatch(r"[0-9a-fA-F]{64}", value) is not None,
                "Malformed checksum: " + rel)
        parts_path = Path(rel)
        require(not parts_path.is_absolute() and ".." not in parts_path.parts and ":" not in rel,
                "Unsafe manifest path: " + rel)
        require(rel not in records, "Duplicate manifest entry: " + rel)
        records[rel] = value.lower()
    require(len(records) >= 20, "Unexpectedly small BDNEP repository manifest")
    # Do not accept a repository that has already drifted from its root manifest.
    for rel, expected in records.items():
        target = repo / Path(rel)
        require(target.is_file() and not target.is_symlink(), "Manifest file absent: " + rel)
        require(sha256(target) == expected, "Manifest checksum mismatch: " + rel)
    return records


def validated_counts(repo: Path) -> tuple[list[str], np.ndarray, np.ndarray, dict]:
    report_path = repo / "results" / "replay_verification_report.json"
    counts_path = repo / "results" / "recomputed_participants.csv"
    require(report_path.is_file() and counts_path.is_file(), "Missing P14 replay outputs")
    report = json.loads(report_path.read_text(encoding="utf-8-sig"))
    require(report.get("status") == P14_SUCCESS and report.get("n_subjects") == 37,
            "Full independent replay has not been validated")

    with counts_path.open(newline="", encoding="utf-8-sig") as source:
        reader = csv.DictReader(source)
        required = {"subject_id", "n_A", "n_B", *EXPECTED_COUNTS}
        require(required.issubset(set(reader.fieldnames or [])), "Missing trial-count columns")
        rows = list(reader)

    subjects = [r["subject_id"] for r in rows]
    require(subjects == EXPECTED_SUBJECTS, "Expected exactly sub-004 to sub-040, in order")

    totals = {name: 0 for name in EXPECTED_COUNTS}
    face: list[int] = []
    car: list[int] = []
    arm_b_sum = 0
    for row in rows:
        sid = row["subject_id"]
        a = int(row["n_A"])
        b = int(row["n_B"])
        require(a == 320 and 0 <= b <= 320, "Invalid trial total: " + sid)
        parts: dict[str, int] = {}
        for name in EXPECTED_COUNTS:
            v = int(row[name])
            require(0 <= v <= 80, "Condition eligibility out of bounds: " + sid + "/" + name)
            parts[name] = v
            totals[name] += v
        require(sum(parts.values()) == b, "Arm B condition mismatch: " + sid)
        face.append(parts["n_B_face"])
        car.append(parts["n_B_car"])
        arm_b_sum += b

    require(totals == EXPECTED_COUNTS, "Frozen Arm B condition totals disagree")
    require(arm_b_sum == 10810, "Frozen Arm B total disagrees")
    idx = subjects.index("sub-016")
    require(face[idx] == 2 and car[idx] == 77,
            "Expected sub-016 count imbalance not verified")

    metadata = {
        "status": "S2_COUNTS_VERIFIED_FROM_P14",
        "participants": 37,
        "arm_a_stimuli": 11840,
        "arm_b_eligible_all_conditions": arm_b_sum,
        "face_eligible": totals["n_B_face"],
        "car_eligible": totals["n_B_car"],
        "scrambled_face_eligible_not_plotted": totals["n_B_scrambled_face"],
        "scrambled_car_eligible_not_plotted": totals["n_B_scrambled_car"],
        "sub_016_face": 2,
        "sub_016_car": 77,
        "participant_table_sha256": sha256(counts_path),
        "replay_report_sha256": sha256(report_path),
        "interpretation": "Descriptive Arm B trial availability; no post-hoc exclusion or artifact rejection",
    }
    return subjects, np.asarray(face), np.asarray(car), metadata


def plot_eligibility(subjects: list[str], face: np.ndarray, car: np.ndarray,
                     png: Path, svg: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xs = np.arange(len(subjects), dtype=float)
    fig, ax = plt.subplots(figsize=(13.2, 5.7))
    ax.scatter(xs - 0.16, face, s=30, marker="o", label="Face", zorder=3)
    ax.scatter(xs + 0.16, car, s=34, marker="s", label="Car", zorder=3)
    ax.axhline(80, color="0.48", linewidth=1, linestyle="--",
               label="Maximum per condition (80)", zorder=1)
    ax.axhline(0, color="0.78", linewidth=0.6)
    idx = subjects.index("sub-016")
    ax.annotate("sub-016: 2 eligible face trials", xy=(xs[idx]-0.16, 2),
                xytext=(xs[idx]+3.3, 20),
                arrowprops=dict(arrowstyle="->", color="0.25", linewidth=0.9),
                fontsize=9, ha="left")
    ax.set_title("Arm B face and car trial availability by participant")
    ax.set_xlabel("Participant")
    ax.set_ylabel("Behaviorally eligible trials per condition (of 80)")
    ax.set_xticks(xs)
    ax.set_xticklabels(subjects, rotation=90, fontsize=7.0)
    ax.set_xlim(-0.8, len(subjects)-0.2)
    ax.set_ylim(-2, 88)
    ax.grid(axis="y", alpha=0.16)
    ax.legend(frameon=False, loc="upper right", ncol=3, fontsize=9)
    fig.tight_layout()
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)


def register_new_artifacts(repo: Path, prior: dict[str, str], additions: list[Path]) -> None:
    manifest_path = repo / "SHA256SUMS.txt"
    updated = dict(prior)
    for path in additions:
        rel = path.relative_to(repo).as_posix()
        require(rel not in prior, "Manifest entry already exists: " + rel)
        updated[rel] = sha256(path)
    # Atomic manifest replacement, without deleting existing verified entries.
    with NamedTemporaryFile("w", encoding="utf-8", dir=repo, prefix=".SHA256SUMS_",
                            suffix=".tmp", delete=False, newline="\n") as tmp:
        temp_path = Path(tmp.name)
        for rel, digest in sorted(updated.items()):
            tmp.write(f"{digest}  {rel}\n")
    os.replace(temp_path, manifest_path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--check", action="store_true", help="Verify input only; do not create files")
    action.add_argument("--run", action="store_true", help="Generate S2 and register hashes")
    args = parser.parse_args()
    repo = args.repo.resolve()
    require(repo.is_dir(), "Repository path does not exist")
    require((repo / "scr" / SCRIPT_NAME).is_file(), "Install S2 generator in scr first")
    require((repo / "scr" / SCRIPT_NAME).resolve() == Path(__file__).resolve(),
            "Run the canonical script installed under repository scr/")
    records = load_checksums(repo)
    subjects, face, car, metadata = validated_counts(repo)

    print("STATUS: S2_SOURCE_COUNTS_VERIFIED")
    print("PARTICIPANTS: 37/37")
    print("FACE_ELIGIBLE:", metadata["face_eligible"])
    print("CAR_ELIGIBLE:", metadata["car_eligible"])
    print("SUB016_FACE/CAR: 2/77")
    if args.check:
        print("PLOT_GENERATION: NOT_EXECUTED")
        return

    outdir = repo / "supplementary" / "figures"
    require(outdir.is_dir(), "Missing canonical supplementary/figures directory")
    png = outdir / f"{FIGURE_NAME}.png"
    svg = outdir / f"{FIGURE_NAME}.svg"
    receipt = outdir / f"{FIGURE_NAME}_provenance.json"
    for path in (png, svg, receipt):
        require(not path.exists(), "Refusing to overwrite existing artifact: " + str(path))

    plot_eligibility(subjects, face, car, png, svg)
    metadata["figure_png_sha256"] = sha256(png)
    metadata["figure_svg_sha256"] = sha256(svg)
    metadata["status"] = "S2_GENERATED_FROM_VERIFIED_P14"
    receipt.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n",
                       encoding="utf-8")
    register_new_artifacts(repo, records, [repo / "scr" / SCRIPT_NAME, png, svg, receipt])
    print("STATUS: S2_GENERATED_FROM_VERIFIED_P14")
    print("PNG:", png)
    print("SVG:", svg)
    print("ROOT_MANIFEST: SHA256SUMS.txt UPDATED")
    print("GITHUB_PUSH: NOT_EXECUTED")


if __name__ == "__main__":
    main()
