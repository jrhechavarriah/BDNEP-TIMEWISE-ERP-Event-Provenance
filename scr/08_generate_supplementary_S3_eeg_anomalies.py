#!/usr/bin/env python3
"""Generate descriptive Supplementary Figure S3 from original NEMAR EEG bytes.

This QC visualization is NOT an artifact removal, participant exclusion rule,
source EEG identity assessment, or independent physiological timing calibration.
It plots selected unfiltered segments of sub-029 and sub-035 in microvolts.
Only verified BIDS .set + .fdt source files are used. Source files are read-only.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

EXPECTED_CONTRACT_SHA256 = "f75505d2b30f57d4b6ce88a5efcd60c2cf27398e7021dbe2545f34dab4d4827c"
# Descriptive segments guided by the original effect-agnostic QC images, NOT the N170 effect.
WINDOWS_S = {"sub-029": (151.0, 162.5), "sub-035": (159.5, 175.0)}
CHANNELS = {"sub-029": ("Oz", "PO7", "PO8", "O1", "O2"),
            "sub-035": ("Oz", "PO7", "PO8", "HEOG_left", "VEOG_lower")}
DESCRIPTIVE_THRESHOLD_UV = 250000.0  # QC illustration only; never an exclusion threshold.


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verified_inputs(contract_path: Path, nemar_dir: Path):
    if sha256(contract_path) != EXPECTED_CONTRACT_SHA256:
        raise RuntimeError("Version-locked source contract SHA256 mismatch")
    contract = json.loads(contract_path.read_text(encoding="utf-8-sig"))
    if contract.get("file_count") != 261:
        raise RuntimeError("Incomplete source contract")
    index = {r["path"]: r for r in contract["files"]}
    if len(index) != 261:
        raise RuntimeError("Duplicated source file paths")
    sources = {}
    for subject in WINDOWS_S:
        files = {}
        for ext in ("set", "fdt"):
            rel = f"{subject}/eeg/{subject}_task-N170_eeg.{ext}"
            entry = index.get("nemar/" + rel)
            if entry is None:
                raise RuntimeError("Source missing from contract: " + rel)
            file = nemar_dir / rel
            if not file.is_file() or sha256(file) != entry["sha256"]:
                raise RuntimeError("Official NEMAR source integrity failure: " + rel)
            files[ext] = file
        sources[subject] = files
    return sources


def load_segment(subject: str, set_file: Path):
    import mne
    import numpy as np
    # preload=False: only requested recording segments are loaded as arrays.
    raw = mne.io.read_raw_eeglab(str(set_file), preload=False, verbose="ERROR")
    try:
        sfreq = float(raw.info["sfreq"])
        if sfreq != 1024.0 or len(raw.ch_names) != 33:
            raise RuntimeError(f"Unexpected EEG sampling/channel contract for {subject}")
        for channel in CHANNELS[subject]:
            if channel not in raw.ch_names:
                raise RuntimeError(f"Missing {subject} channel: {channel}")
        start_s, stop_s = WINDOWS_S[subject]
        start = round(start_s * sfreq)
        stop = round(stop_s * sfreq) + 1
        if start < 0 or stop > raw.n_times or stop <= start:
            raise RuntimeError(f"QC window outside {subject} recording")
        # MNE signal values are volts. Multiply exactly once to obtain µV.
        samples = raw.get_data(start=start, stop=stop) * 1e6
        if not np.isfinite(samples).all():
            raise RuntimeError(f"Nonfinite unfiltered source values: {subject}")
        t = np.arange(start, stop, dtype=float) / sfreq
        ch = {name: i for i, name in enumerate(raw.ch_names)}
        selected = {name: samples[ch[name], :] for name in CHANNELS[subject]}
        median = np.median(samples, axis=0)
        dispersion = np.std(samples, axis=0)
        exceed = (np.abs(samples) >= DESCRIPTIVE_THRESHOLD_UV).sum(axis=0)
        return {"t": t, "selected": selected,
                "median": median, "sd": dispersion, "exceed": exceed,
                "n_recorded_channels": len(raw.ch_names), "fs": sfreq,
                "samples": len(t), "exceed_max_channels": int(exceed.max()),
                "max_selected_abs_uV": float(max(np.max(np.abs(v)) for v in selected.values()))}
    finally:
        raw.close()


def render_s3(segments: dict, target: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(15.4, 11.3),
                             sharex="col", gridspec_kw={"height_ratios": [1.25, 1.0, 0.55]})
    for col, subject in enumerate(("sub-029", "sub-035")):
        s = segments[subject]
        ax = axes[0, col]
        for label, series in s["selected"].items():
            ax.plot(s["t"], series, linewidth=0.8, label=label)
        ax.set_title(f"{'A' if col == 0 else 'B'}. {subject}: selected unfiltered channels")
        ax.set_ylabel("Unfiltered signal (µV)")
        ax.legend(loc="best", fontsize=8, ncol=3, frameon=False)
        ax = axes[1, col]
        ax.plot(s["t"], s["median"], linewidth=0.8, label="Across-channel median")
        ax.plot(s["t"], s["sd"], linewidth=0.8, label="Across-channel SD")
        ax.set_ylabel("Across-channel summary (µV)")
        ax.legend(loc="best", fontsize=8, frameon=False)
        ax = axes[2, col]
        ax.plot(s["t"], s["exceed"], linewidth=0.8)
        ax.set_ylabel("Channels above\ndescriptive threshold")
        ax.set_ylim(-0.8, 34)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
        ax.set_xlabel("Recording time (s)")
        for row in range(3):
            axes[row, col].set_xlim(WINDOWS_S[subject])
            axes[row, col].grid(axis="y", alpha=0.15)
    fig.suptitle("Supplementary Figure S3. Representative unfiltered EEG signal anomalies",
                 fontsize=14)
    fig.text(0.5, 0.008,
             "Descriptive QC only; 250,000 µV is not an artifact-rejection or exclusion criterion.",
             ha="center", va="bottom", fontsize=9)
    fig.tight_layout(rect=[0, 0.032, 1, 0.965])
    fig.savefig(target, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--contract", required=True, type=Path)
    ap.add_argument("--nemar-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()
    files = verified_inputs(args.contract, args.nemar_dir)
    segments = {s: load_segment(s, files[s]["set"]) for s in WINDOWS_S}
    if args.out_dir.exists():
        raise RuntimeError("Refusing to overwrite output folder: " + str(args.out_dir))
    args.out_dir.mkdir(parents=True, exist_ok=False)
    figure = args.out_dir / "Supplementary_Figure_S3_Unfiltered_EEG_Anomalies.png"
    render_s3(segments, figure)
    summary = {
        "status": "S3_RENDERED_FROM_VERIFIED_NEMAR_UNFILTERED_EEG",
        "figure": figure.name, "figure_sha256": sha256(figure),
        "source_dataset": "nm000132", "dataset_version": "v1.1.1",
        "source_contract_sha256": EXPECTED_CONTRACT_SHA256,
        "window_s": {s: list(w) for s, w in WINDOWS_S.items()},
        "source_sha256": {s: {ext: sha256(p) for ext, p in fs.items()} for s, fs in files.items()},
        "diagnostics": {s: {k: v for k, v in data.items() if k not in ("t", "selected", "median", "sd", "exceed")}
                        for s, data in segments.items()},
        "channels_shown": {s: list(ch) for s, ch in CHANNELS.items()},
        "all_33_channels_used_for_summary": True,
        "descriptive_threshold_uV": DESCRIPTIVE_THRESHOLD_UV,
        "artifact_rejection_or_exclusion": False,
        "not_an_ERP_result": True,
        "no_cross_source_byte_identity_claim": True,
    }
    (args.out_dir / "Supplementary_S3_provenance.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("STATUS: S3_RENDERED_FROM_VERIFIED_NEMAR_UNFILTERED_EEG")
    print("SOURCES: sub-029 and sub-035; original unfiltered NEMAR bytes")
    for subject, d in segments.items():
        print(f"{subject}: selected channels={list(d['selected'])}; "
              f"max descriptive threshold channels={d['exceed_max_channels']}")
    print("FIGURE:", figure)
    print("SHA256:", summary["figure_sha256"])


if __name__ == "__main__":
    main()
