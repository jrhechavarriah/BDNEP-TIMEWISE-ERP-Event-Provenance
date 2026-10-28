#!/usr/bin/env python3
"""IOSA D5-K10-7-P14, INTERNAL independent replay from verified public source bytes.

Frozen sources: NEMAR ERP CORE v1.1.1 (sub-004..040), OSF original/shifted
EEGLAB event headers, P11 manifest.  This does NOT establish EEG byte identity
between OSF and NEMAR, true physical visual onset, or permission to publish.

This deliberately does NOT import the original D5K104 numerical implementation
nor read archived contrasts until after new source-derived estimates are computed.
"""
from __future__ import annotations
import argparse
import bisect
import csv
import hashlib
import json
import math
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path

import mne
import numpy as np
import scipy
from scipy.io import loadmat

EXPECTED = {
    "contract": "8448aa61d793bebadbf193ed81a46de6e14be2d7eb9192de62e47c2a4b99d161",
    "protocol": "38db77378341b403fa2edc78bf5cd8480b073f15cff9bdb92db26ffa4ffadcb1",
    "cohort": "84d372485eb83f334f2674a0554d3d31434cd719431f1a2e1d840e477e392ade",
    "reference_master": "db678ef566d38b0cdd843f96ad5a42a262f416a3aa2172eac2f8a54cebaa0a58",
    "reference_summary": "e1c8f406c1c38bcf3c8e85e34d77b6f049f47cf27a524276d7036fc647dffab8",
}
IDS = [f"sub-{i:03d}" for i in range(4, 41)]
CONDS = ("face", "car", "scrambled_face", "scrambled_car")
CHANNELS = ("PO8", "PO7")
EOG = ("HEOG_left", "HEOG_right", "VEOG_lower")
FS = 1024.0
TOL_UV = 1e-6  # Frozen P8-R2 cross-check tolerance, not a new statistical threshold.
CHUNK = 1024 * 1024


def require(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def as_int(value):
    arr = np.asarray(value)
    require(arr.size == 1, "NON_SCALAR_EVENT_FIELD")
    try:
        dec = Decimal(str(arr.item()).strip())
    except (InvalidOperation, ValueError) as exc:
        raise RuntimeError("NON_NUMERIC_EVENT_FIELD") from exc
    require(dec.is_finite() and dec == dec.to_integral_value(),
            "NON_INTEGER_EVENT_FIELD")
    return int(dec)


def parse_bids_events(path):
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        require({"sample", "value", "trial_type", "event_type"}.issubset(
            set(reader.fieldnames or ())), "BIDS_EVENT_COLUMNS_MISSING")
        rows = list(reader)
    events = []
    for row in rows:
        typ = row["trial_type"].strip().lower()
        label = row["event_type"].strip()
        require(typ in ("stimulus", "response"), "INVALID_BIDS_EVENT_TYPE")
        sample1 = as_int(row["sample"])
        require(sample1 >= 1, "INVALID_BIDS_SAMPLE")
        code = as_int(row["value"])
        if typ == "stimulus":
            bounds = {"face": (1, 40), "car": (41, 80),
                      "scrambled_face": (101, 140),
                      "scrambled_car": (141, 180)}
            require(label in bounds and bounds[label][0] <= code <= bounds[label][1],
                    "STIMULUS_CODE_MAPPING_MISMATCH")
        else:
            require((label, code) in (("correct_response", 201),
                                      ("incorrect_response", 202)),
                    "RESPONSE_CODE_MAPPING_MISMATCH")
        events.append((code, sample1, typ, label))
    require(all(a[1] <= b[1] for a, b in zip(events, events[1:])),
            "EVENT_ORDER_NOT_MONOTONIC")
    require(len([e for e in events if e[2] == "stimulus"]) == 320,
            "STIMULUS_COUNT_NOT_320")
    require(Counter(e[3] for e in events if e[2] == "stimulus") ==
            Counter(dict.fromkeys(CONDS, 80)), "CONDITION_COUNTS_NOT_80")
    return events


def parse_osf_events(path):
    eeg = loadmat(str(path), variable_names=["EEG"], struct_as_record=False,
                  squeeze_me=True).get("EEG")
    require(eeg is not None and hasattr(eeg, "event") and hasattr(eeg, "srate"),
            "OSF_EEGLAB_HEADER_INVALID")
    srate = float(np.asarray(eeg.srate).item())
    require(srate == FS, "OSF_SAMPLE_RATE_NOT_1024")
    return [(as_int(x.type), as_int(x.latency))
            for x in np.asarray(eeg.event, dtype=object).reshape(-1)]


def validate_event_identity(bids, original, shifted, subject):
    require(len(bids) == len(original) == len(shifted),
            f"EVENT_COUNT_MISMATCH_{subject}")
    stim, response = 0, 0
    for (code, sample, typ, _), (oc, osample), (sc, ssample) in zip(
            bids, original, shifted):
        require(code == oc == sc and sample == osample,
                f"OSF_BIDS_DIGITAL_IDENTITY_FAILED_{subject}")
        if typ == "stimulus":
            require(ssample - osample == 26,
                    f"STIMULUS_SHIFT_FAILED_{subject}")
            stim += 1
        else:
            require(ssample == osample, f"RESPONSE_SHIFT_FAILED_{subject}")
            response += 1
    require(stim == 320 and stim + response == len(bids),
            f"EVENT_COUNTS_FAILED_{subject}")
    return {"subject_id": subject, "all_events": len(bids),
            "stimuli_shifted_26": stim, "responses_unchanged": response}


def classify_b(bids, sid):
    stimuli = [r for r in bids if r[2] == "stimulus"]
    responses = [r for r in bids if r[2] == "response"]
    centers = [r[1] - 1 for r in stimuli]
    require(len(set(centers)) == 320, f"DUPLICATE_STIMULI_{sid}")
    linked = [[] for _ in stimuli]
    orphans = 0
    for response in responses:
        pos = bisect.bisect_left(centers, response[1] - 1) - 1
        if pos < 0:
            orphans += 1
        else:
            linked[pos].append(response)
    flags = []
    classifications = Counter()
    for stim, group in zip(stimuli, linked):
        if not group:
            state = "NO_RESPONSE"
        elif len(group) > 1:
            state = "MULTIPLE_RESPONSES"
        elif group[0][3] == "incorrect_response":
            state = "INCORRECT_RESPONSE"
        else:
            rt = ((group[0][1] - 1) - (stim[1] - 1)) * 1000.0 / FS
            if rt < 200 or rt > 1000:
                state = "CORRECT_BUT_RT_OUTSIDE_200_1000_MS"
            else:
                state = "BEHAVIORALLY_ELIGIBLE_PRE_ARTIFACT"
        flags.append(state == "BEHAVIORALLY_ELIGIBLE_PRE_ARTIFACT")
        classifications[state] += 1
    mask = np.asarray(flags, dtype=bool)
    require(len(mask) == 320, "ARM_B_MASK_LENGTH_MISMATCH")
    return stimuli, mask, orphans, classifications


def source_index(manifest, root, subject):
    mapping = {item["path"]: item for item in manifest["files"]}
    files = {}
    nemar_eeg = root["nemar"] / subject / "eeg"
    stem = f"{subject}_task-N170"
    for name, extension in (("bids_set", "_eeg.set"),
                            ("bids_fdt", "_eeg.fdt"),
                            ("bids_json", "_eeg.json"),
                            ("bids_channels", "_channels.tsv"),
                            ("bids_events", "_events.tsv")):
        path = nemar_eeg / (stem + extension)
        rel = "nemar/" + path.relative_to(root["nemar"]).as_posix()
        files[name] = (rel, path)
    number = int(subject[4:])
    for kind, suffix in (("original", ""), ("shifted", "_shifted")):
        path = root["osf"] / subject / f"{number}_N170{suffix}.set"
        rel = "osf/" + path.relative_to(root["osf"]).as_posix()
        files[kind] = (rel, path)
    for name, (rel, path) in files.items():
        require(rel in mapping, f"NOT_IN_PINNED_CONTRACT_{subject}_{name}")
        require(path.is_file(), f"SOURCE_MISSING_{subject}_{name}: {path}")
        require(sha(path) == mapping[rel]["sha256"],
                f"SOURCE_SHA256_CHANGED_{subject}_{name}")
    return {name: path for name, (rel, path) in files.items()}


def load_reference(reference):
    required = ("source_manifest_v1.json", "IOSA_D5K10_Protocol_v0_1.json",
                "extended_cohort_lock.json", "table_01_participant_master.csv",
                "table_02_results_summary.csv")
    for name in required:
        require((reference / name).is_file(), "MISSING_REFERENCE_" + name)
    for key, name in (("contract", "source_manifest_v1.json"),
                      ("protocol", "IOSA_D5K10_Protocol_v0_1.json"),
                      ("cohort", "extended_cohort_lock.json")):
        require(sha(reference / name) == EXPECTED[key],
                "FROZEN_REFERENCE_HASH_MISMATCH_" + key)
    for key, name in (("reference_master", "table_01_participant_master.csv"),
                      ("reference_summary", "table_02_results_summary.csv")):
        require(sha(reference / name) == EXPECTED[key],
                "FROZEN_TABLE_HASH_MISMATCH_" + key)
    contract = json.loads((reference / "source_manifest_v1.json").read_text(encoding="utf-8-sig"))
    protocol = json.loads((reference / "IOSA_D5K10_Protocol_v0_1.json").read_text(encoding="utf-8-sig"))
    cohort = json.loads((reference / "extended_cohort_lock.json").read_text(encoding="utf-8-sig"))
    require(contract.get("file_count") == 261 and len(contract["files"]) == 261,
            "PORTABLE_SOURCE_CONTRACT_INCOMPLETE")
    require(cohort.get("status") == "LOCKED_COMPUTATIONALLY" and
            cohort["cohort"]["included_ids"] == IDS,
            "COHORT_LOCK_INVALID")
    a = protocol["frozen_analysis"]
    require(a["epoch_ms"] == [-200, 800] and a["baseline_ms"] == [-200, 0] and
            a["primary_window_ms"] == [110, 150] and
            a["bootstrap_seed"] == 20261019 and
            a["bootstrap_replicates"] == 10000 and
            a["no_ICA"] and a["no_artifact_rejection"] and
            a["no_bad_channel_interpolation"], "FROZEN_PROTOCOL_DRIFT")
    with (reference / "table_01_participant_master.csv").open(
            encoding="utf-8-sig", newline="") as stream:
        master = list(csv.DictReader(stream))
    with (reference / "table_02_results_summary.csv").open(
            encoding="utf-8-sig", newline="") as stream:
        summary = list(csv.DictReader(stream))
    require([r["subject_id"] for r in master] == IDS and len(summary) == 2,
            "ARCHIVED_REFERENCE_COHORT_INCOMPLETE")
    require(all(r.get("subject_id") for r in master), "REFERENCE_MISSING_IDS")
    require(all("PO8_A_sign_flip" in row and "PO8_B_sign_flip" in row
                for row in master), "FROZEN_PO8_SIGN_FLIP_COLUMNS_MISSING")
    return contract, cohort, master, summary


def channel_metadata_ok(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    eeg = {r["name"] for r in rows if r["type"].upper() == "EEG"}
    eog = {r["name"] for r in rows if r["type"].upper() == "EOG"}
    require(len(eeg) == 30 and eog == set(EOG) and
            CHANNELS[0] in eeg and CHANNELS[1] in eeg,
            "BIDS_CHANNEL_METADATA_NOT_30EEG_3EOG")


def verify_archived_po8_sign_flips(result, reference_row, subject):
    """Compare only the PO8 sign-flip flags present in the frozen CSV.

    PO7 original/corrected amplitudes are verified separately for both arms.
    This schema fix does not change any EEG, trial, or statistical computation.
    """
    for arm in ("A", "B"):
        column = f"PO8_{arm}_sign_flip"
        flag = reference_row.get(column)
        require(isinstance(flag, str) and flag.strip().lower() in ("true", "false"),
                f"FROZEN_PO8_SIGN_FLIP_FLAG_INVALID_{subject}_{arm}")
        require(result[column] == (flag.strip().lower() == "true"),
                f"ARCHIVED_SIGN_FLIP_MISMATCH_{subject}_PO8_{arm}")


def compute_subject(subject, contract, source_roots, reference_row, frozen_subject):
    paths = source_index(contract, source_roots, subject)
    require(sha(paths["bids_events"]) ==
            frozen_subject["sha256"]["bids_events_recorded"],
            "COHORT_LOCK_BIDS_EVENT_HASH_MISMATCH_" + subject)
    require(sha(paths["original"]) ==
            frozen_subject["sha256"]["original_osf_set"] and
            sha(paths["shifted"]) ==
            frozen_subject["sha256"]["shifted_osf_set"],
            "COHORT_LOCK_OSF_HASH_MISMATCH_" + subject)
    channel_metadata_ok(paths["bids_channels"])
    with paths["bids_json"].open(encoding="utf-8-sig") as stream:
        meta = json.load(stream)
    require(int(meta.get("SamplingFrequency", -1)) == 1024 and
            int(meta.get("EEGChannelCount", -1)) == 30 and
            int(meta.get("EOGChannelCount", -1)) == 3,
            "BIDS_EEG_METADATA_MISMATCH")
    bids = parse_bids_events(paths["bids_events"])
    identity = validate_event_identity(
        bids, parse_osf_events(paths["original"]),
        parse_osf_events(paths["shifted"]), subject)
    require(len(bids) == frozen_subject["events"] and
            identity["stimuli_shifted_26"] == frozen_subject["stimuli_shifted_26_samples"] and
            identity["responses_unchanged"] == frozen_subject["responses_shifted_0_samples"],
            "FROZEN_PARTICIPANT_EVENT_TOTAL_MISMATCH_" + subject)
    stimuli, selected_b, orphans, classifications = classify_b(bids, subject)
    require(len(bids) == 320 + identity["responses_unchanged"],
            "BIDS_EVENT_TALLY_INVALID")
    code_ids = {label: i + 1 for i, label in enumerate(CONDS)}
    samples0 = np.asarray([r[1] - 1 for r in stimuli], dtype=np.int64)
    codes = np.asarray([code_ids[r[3]] for r in stimuli], dtype=np.int64)
    labels = np.asarray([r[3] for r in stimuli])
    require(np.all(samples0 >= 0) and np.all(np.diff(samples0) > 0),
            "NONSTRICT_STIMULUS_SAMPLES")

    raw = mne.io.read_raw_eeglab(str(paths["bids_set"]), preload=True, verbose="ERROR")
    require(raw.info["sfreq"] == FS, "MNE_SAMPLING_RATE_MISMATCH")
    require(set(EOG).issubset(raw.ch_names), "MNE_EOG_CHANNELS_MISSING")
    raw.set_channel_types({c: "eog" for c in EOG}, verbose="ERROR")
    eeg_picks = mne.pick_types(raw.info, eeg=True, eog=False, exclude=[])
    require(len(eeg_picks) == 30, "REFERENCE_MUST_USE_EXACTLY_30_EEG")
    raw.filter(l_freq=0.1, h_freq=30.0, picks="eeg", method="iir",
               iir_params={"order": 2, "ftype": "butter", "output": "sos"},
               phase="zero", verbose="ERROR")
    raw.set_eeg_reference(ref_channels="average", projection=False,
                          verbose="ERROR")

    records = {}
    waveforms = {}
    epoch_time = None
    selected_condition_counts = Counter(labels[selected_b])
    for branch, shift in (("original", 0), ("corrected", 26)):
        event_matrix = np.column_stack([samples0 + shift,
                                        np.zeros(len(samples0), dtype=np.int64),
                                        codes])
        epochs = mne.Epochs(raw, event_matrix, event_id=code_ids,
                            tmin=-0.2, tmax=0.8,
                            baseline=(-0.2, 0.0),
                            picks=list(CHANNELS), preload=True,
                            reject=None, flat=None, reject_by_annotation=False,
                            detrend=None, verbose="ERROR")
        require(len(epochs) == 320 and
                np.array_equal(epochs.events, event_matrix),
                f"EPOCH_IDENTITY_DROPPED_{subject}_{branch}")
        window = (epochs.times >= 0.110) & (epochs.times <= 0.150)
        require(window.sum() == 41, "N170_WINDOW_NOT_41_SAMPLES")
        if epoch_time is None:
            epoch_time = epochs.times.copy()
        else:
            require(np.array_equal(epoch_time, epochs.times),
                    "BRANCH_EPOCH_AXIS_MISMATCH")
        data = epochs.get_data(copy=False) * 1e6  # microvolts
        require(np.isfinite(data).all(), f"NONFINITE_EPOCHS_{subject}")
        for arm, trial_mask in (("A", np.ones(320, dtype=bool)),
                                ("B", selected_b)):
            for chan_idx, ch in enumerate(CHANNELS):
                face = trial_mask & (labels == "face")
                car = trial_mask & (labels == "car")
                require(face.any() and car.any(),
                        f"ZERO_FACE_OR_CAR_EPOCHS_{subject}_{arm}")
                face_mean = data[face, chan_idx, :].mean(axis=0)
                car_mean = data[car, chan_idx, :].mean(axis=0)
                contrast = float((face_mean - car_mean)[window].mean())
                records[(arm, branch, ch)] = contrast
                if arm == "A" and ch == "PO8":
                    waveforms[branch] = (face_mean - car_mean).copy()
    del raw
    require(epoch_time is not None, "MISSING_EPOCH_TIME_VECTOR")
    result = {"subject_id": subject, "n_A": 320,
              "n_B": int(selected_b.sum()),
              "n_B_face": selected_condition_counts["face"],
              "n_B_car": selected_condition_counts["car"],
              "n_B_scrambled_face": selected_condition_counts["scrambled_face"],
              "n_B_scrambled_car": selected_condition_counts["scrambled_car"],
              "n_orphans": orphans,
              "n_responses": identity["responses_unchanged"]}
    for ch in CHANNELS:
        for arm in ("A", "B"):
            for branch in ("original", "corrected"):
                result[f"{ch}_{arm}_{branch}_uV"] = records[(arm, branch, ch)]
            a = records[(arm, "original", ch)]
            b = records[(arm, "corrected", ch)]
            result[f"{ch}_{arm}_delta_uV"] = b - a
            result[f"{ch}_{arm}_sign_flip"] = bool(np.sign(a) != np.sign(b))
    numerical_checks = 0
    for col in (f"{ch}_{arm}_{branch}_uV" for ch in CHANNELS
                for arm in ("A", "B") for branch in ("original", "corrected")):
        expected = float(reference_row[col])
        actual = result[col]
        require(abs(actual - expected) <= TOL_UV,
                f"ARCHIVED_NUMERICAL_MISMATCH_{subject}_{col}: actual={actual} reference={expected}")
        numerical_checks += 1
    for col in ("n_A", "n_B", "n_B_face", "n_B_car",
                "n_B_scrambled_face", "n_B_scrambled_car"):
        require(int(result[col]) == int(reference_row[col]),
                f"ARCHIVED_BEHAVIOR_COUNTS_MISMATCH_{subject}_{col}")
    # The frozen reference-master schema has PO8 sign-flip flags only.
    # PO7 contrasts remain numerically validated above; its computed flags
    # are retained descriptively, but must not be fetched from absent columns.
    verify_archived_po8_sign_flips(result, reference_row, subject)
    require(result["n_B_face"] >= 0 and result["n_B_car"] >= 0,
            "IMPOSSIBLE_TRIAL_COUNTS")
    print(f"[PASS] {subject} signal+events, {numerical_checks} contrasts, B={result['n_B']}/320",
          flush=True)
    return result, identity, waveforms, epoch_time, classifications


def compare_summary(results, archive):
    a = [r[f"PO8_A_delta_uV"] for r in results]
    b = [r[f"PO8_B_delta_uV"] for r in results]
    rng = np.random.default_rng(20261019)
    idx = rng.integers(0, len(a), size=(10000, len(a)))
    interval = np.percentile(np.asarray(a)[idx].mean(axis=1), [2.5, 97.5])
    require(len(results) == 37, "FULL_COHORT_NOT_37")
    expected = {r["analysis"]: r for r in archive}
    for arm, values, label in (("A", a, "PO8 Arm A primary"),
                               ("B", b, "PO8 Arm B secondary")):
        frozen = expected[label]
        mean = float(np.mean(values))
        require(abs(mean - float(frozen["mean_delta_uV"])) <= TOL_UV,
                f"GROUP_MEAN_CHANGED_{arm}")
        require(sum(r[f"PO8_{arm}_sign_flip"] for r in results) ==
                int(frozen["sign_flips"]), f"GROUP_SIGN_FLIPS_CHANGED_{arm}")
    require(sum(r["n_A"] for r in results) == 11840 and
            sum(r["n_B"] for r in results) == 10810,
            "TOTAL_TRIAL_COUNT_CHANGED")
    require(np.allclose(interval,
                         [float(expected["PO8 Arm A primary"]["ci95_low_uV"]),
                          float(expected["PO8 Arm A primary"]["ci95_high_uV"])],
                         rtol=0, atol=TOL_UV), "BOOTSTRAP_CI_CHANGED")
    return {"n_participants": 37, "n_A_stimuli": 11840,
            "PO7_A_mean_delta_uV_descriptive": float(np.mean(
                [r["PO7_A_delta_uV"] for r in results])),
            "PO7_B_mean_delta_uV_descriptive": float(np.mean(
                [r["PO7_B_delta_uV"] for r in results])),
            "n_B_stimuli": 10810,
            "PO8_A_mean_delta_uV": float(np.mean(a)),
            "PO8_A_bootstrap_ci95_uV": interval.tolist(),
            "PO8_A_sign_flips": int(sum(r["PO8_A_sign_flip"] for r in results)),
            "PO8_B_mean_delta_uV_descriptive": float(np.mean(b)),
            "PO8_B_sign_flips_descriptive": int(sum(r["PO8_B_sign_flip"] for r in results)),
            "primary_bootstrap_replicates": 10000,
            "primary_bootstrap_seed": 20261019,
            "secondary_inference": "DESCRIPTIVE_ONLY"}


def write_csv(path, rows):
    require(bool(rows), "CANNOT_WRITE_EMPTY_CSV")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def figures(out, results, waves, time_s, full):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figdir = out / "figures"
    figdir.mkdir(exist_ok=False)
    if full:
        original = np.mean([waves[r["subject_id"]]["original"] for r in results], axis=0)
        corrected = np.mean([waves[r["subject_id"]]["corrected"] for r in results], axis=0)
        fig, ax = plt.subplots(figsize=(9, 4.6))
        ax.plot(time_s * 1000, original, label="Original digital markers", lw=1.8)
        ax.plot(time_s * 1000, corrected, label="Author-corrected digital markers", lw=1.8)
        ax.axvspan(110, 150, alpha=0.12, color="gray", label="Prespecified ROI")
        ax.axhline(0, color="black", lw=0.7)
        ax.set(xlabel="Time relative to digital marker (ms)",
               ylabel="Face minus car, PO8 (µV)",
               title="Grand-average PO8 ERP contrast (n = 37)")
        ax.legend(loc="best", frameon=False)
        fig.tight_layout()
        fig.savefig(figdir / "figure_01_PO8_group_ERP.png", dpi=220)
        plt.close(fig)
    fig, ax = plt.subplots(figsize=(9.8, 5.0))
    x = np.arange(len(results))
    ax.axhline(0, color="black", lw=0.7)
    ax.scatter(x, [r["PO8_A_delta_uV"] for r in results], s=19)
    ax.set_xticks(x)
    ax.set_xticklabels([r["subject_id"].replace("sub-", "") for r in results],
                       rotation=90, fontsize=7)
    ax.set(xlabel="Participant", ylabel="Corrected minus original contrast (µV)",
           title="Participant-paired PO8 Arm A estimates")
    fig.tight_layout()
    fig.savefig(figdir / "figure_02_participant_deltas.png", dpi=220)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(9.8, 4.5))
    ax.bar(x, [r["n_B"] for r in results], width=0.72)
    ax.axhline(320, linestyle="--", lw=0.9, color="gray")
    ax.set_xticks(x)
    ax.set_xticklabels([r["subject_id"].replace("sub-", "") for r in results],
                       rotation=90, fontsize=7)
    ax.set(xlabel="Participant", ylabel="Behaviorally eligible N170 trials",
           title="Arm B trial availability before artifact rejection")
    fig.tight_layout()
    fig.savefig(figdir / "figure_03_behavioral_trials.png", dpi=220)
    plt.close(fig)


def selftest():
    # Synthetic, non-scientific unit checks only; no frozen contrasts as input.
    assert as_int(np.array(26.0)) == 26
    try:
        as_int(1.5)
    except RuntimeError:
        pass
    else:
        raise RuntimeError("NON_INTEGER_EVENT_NOT_REJECTED")
    events = []
    for index in range(320):
        sample = 1000 + 1500 * index
        events.append((1 if index < 160 else 41, sample,
                       "stimulus", "face" if index < 160 else "car"))
        if index < 319:
            events.append((201, sample + 300, "response", "correct_response"))
        else:
            events.append((202, sample + 300, "response", "incorrect_response"))
    _, mask, orphans, counts = classify_b(events, "synthetic")
    assert len(mask) == 320 and mask.sum() == 319 and orphans == 0
    assert counts["INCORRECT_RESPONSE"] == 1
    result = validate_event_identity(
        events, [(e[0], e[1]) for e in events],
        [(e[0], e[1] + (26 if e[2] == "stimulus" else 0))
         for e in events], "synthetic")
    assert result["stimuli_shifted_26"] == 320
    assert result["responses_unchanged"] == 320
    # Regression: archived master has PO8 flags, but *no* PO7 flags.
    # This must pass, and a changed PO8 flag must still fail closed.
    test_result = {"PO8_A_sign_flip": True, "PO8_B_sign_flip": False,
                   "PO7_A_sign_flip": False, "PO7_B_sign_flip": True}
    test_reference = {"PO8_A_sign_flip": "True", "PO8_B_sign_flip": "False"}
    verify_archived_po8_sign_flips(test_result, test_reference, "synthetic")
    for corrupt_result, corrupt_reference in (
        ({**test_result, "PO8_A_sign_flip": False}, test_reference),
        (test_result, {"PO8_A_sign_flip": "True"}),
    ):
        try:
            verify_archived_po8_sign_flips(corrupt_result, corrupt_reference,
                                          "synthetic")
        except RuntimeError:
            pass
        else:
            raise RuntimeError("PO8_SIGN_FLIP_NEGATIVE_CONTROL_DID_NOT_FAIL")
    print("[PASS] P14-R1 PO8 flag schema regression tests")
    print("[PASS] P14 synthetic unit tests (no raw EEG)")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", type=Path, required=True,
                   help="Reference directory for the reproducible analysis")
    p.add_argument("--nemar", type=Path)
    p.add_argument("--osf", type=Path)
    p.add_argument("--output-root", type=Path)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--self-test", action="store_true")
    mode.add_argument("--verify-contract", action="store_true")
    mode.add_argument("--subject", choices=IDS)
    mode.add_argument("--all", action="store_true")
    args = p.parse_args()

    if args.self_test:
        selftest()
        return
    contract, cohort, master, archive = load_reference(args.reference)
    print("[PASS] Public package reference integrity SHA256 validated")
    if args.verify_contract:
        print("STATUS: PUBLIC_PACKAGE_REFERENCES_VERIFIED_ONLY")
        print("SCIENTIFIC_REPLAY: NOT_EXECUTED")
        return
    require(args.nemar is not None and args.osf is not None and
            args.output_root is not None,
            "NEMAR_OSF_OUTPUT_ROOT_ARGUMENTS_REQUIRED")
    roots = {"nemar": args.nemar, "osf": args.osf}
    participants = IDS if args.all else [args.subject]
    refs = {r["subject_id"]: r for r in master}
    locks = {r["subject_id"]: r for r in cohort["cohort"]["participants"]}
    require(set(locks) == set(IDS), "COHORT_PARTICIPANT_LOCK_MAPPING_INVALID")
    for root_name in ("participants.tsv", "dataset_description.json"):
        entry = next((x for x in contract["files"]
                      if x["path"] == "nemar/" + root_name), None)
        require(entry is not None and
                sha(roots["nemar"] / root_name) == entry["sha256"],
                "NEMAR_ROOT_METADATA_HASH_MISMATCH_" + root_name)
    results, identities, waveforms = [], [], {}
    global_time = None
    total_class = Counter()
    total_orphans = 0

    for sid in participants:
        result, identity, waves, times, classification = compute_subject(
            sid, contract, roots, refs[sid], locks[sid])
        results.append(result)
        identities.append(identity)
        waveforms[sid] = waves
        require(global_time is None or np.array_equal(global_time, times),
                "SUBJECT_TIME_AXES_DIVERGE")
        global_time = times
        total_class.update(classification)
        total_orphans += result["n_orphans"]
    full = len(results) == 37
    group = compare_summary(results, archive) if full else None
    if full:
        require(sum(x["all_events"] for x in identities) == 23654 and
                sum(x["responses_unchanged"] for x in identities) == 11814 and
                sum(x["stimuli_shifted_26"] for x in identities) == 11840,
                "COHORT_DIGITAL_EVENT_TOTALS_CHANGED")
        require(total_orphans == 67, "ORPHAN_RESPONSES_DRIFT")
        expected_counts = {
            "BEHAVIORALLY_ELIGIBLE_PRE_ARTIFACT": 10810,
            "INCORRECT_RESPONSE": 876,
            "NO_RESPONSE": 100,
            "MULTIPLE_RESPONSES": 7,
            "CORRECT_BUT_RT_OUTSIDE_200_1000_MS": 47,
        }
        require(dict(total_class) == expected_counts,
                "ARM_B_CLASSIFICATION_DISTRIBUTION_DRIFT: " + str(dict(total_class)))
        require({x["subject_id"] for x in results} == set(IDS),
                "COHORT_COVERAGE_MISMATCH")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f")
    out = args.output_root / (f"IOSA_D5K107_P14_REPLAY_{'FULL' if full else 'SMOKE'}_{stamp}")
    require(not out.exists(), "OUTPUT_DIRECTORY_ALREADY_EXISTS")
    out.mkdir(parents=True, exist_ok=False)
    write_csv(out / "recomputed_participants.csv", results)
    write_csv(out / "recomputed_digital_event_provenance.csv", identities)
    times_ms = (global_time * 1000.0).tolist()
    if full:
        group_wave_original = np.mean(
            [waveforms[sid]["original"] for sid in participants], axis=0)
        group_wave_corrected = np.mean(
            [waveforms[sid]["corrected"] for sid in participants], axis=0)
        write_csv(out / "figure_01_PO8_waveform_data.csv", [
            {"time_ms": t, "original_face_minus_car_uV": float(a),
             "corrected_face_minus_car_uV": float(b)}
            for t, a, b in zip(times_ms, group_wave_original, group_wave_corrected)
        ])
    np.savez_compressed(
        out / "subject_PO8_contrast_waveforms.npz",
        times_s=global_time,
        subject_ids=np.asarray(participants),
        po8_original_uV=np.asarray([waveforms[s]["original"] for s in participants]),
        po8_corrected_uV=np.asarray([waveforms[s]["corrected"] for s in participants]),
    )
    figures(out, results, waveforms, global_time, full=full)
    report = {
        "stage": "D5-K10-7-P14", "status": (
            "FULL_INDEPENDENT_REPLAY_MATCHES_FROZEN_NUMERICAL_REFERENCES"
            if full else "SUBJECT_INTEGRATION_MATCHES_FROZEN_NUMERICAL_REFERENCES"),
        "subjects": participants, "n_subjects": len(participants),
        "validated_subject_source_files": len(participants) * 7,
        "sample_event_shift": {"stimuli": 26, "responses": 0},
        "source_contract_sha256": EXPECTED["contract"],
        "frozen_protocol_sha256": EXPECTED["protocol"],
        "cohort_lock_sha256": EXPECTED["cohort"],
        "numerical_tolerance_uV": TOL_UV,
        "group": group,
        "arm_b_classification_counts": dict(sorted(total_class.items())),
        "orphan_responses": total_orphans,
        "mne_version": mne.__version__, "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "osf_nemar_signal_byte_identity": "NOT_ESTABLISHED",
        "physical_display_onset": "NOT_VERIFIED",
        "reuse_of_frozen_contrasts_as_input": False,
        "archive_used_only_for_post_computation_comparison": True,
        "no_ica_no_rejection_no_interpolation": True,
    }
    report_path = out / "replay_verification_report.json"
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                           encoding="utf-8")
    paths = sorted(x for x in out.rglob("*") if x.is_file())
    with (out / "SHA256SUMS.txt").open("x", encoding="utf-8") as stream:
        for path in paths:
            stream.write(f"{sha(path)}  {path.relative_to(out).as_posix()}\n")
    print("=" * 70)
    print("STATUS:", report["status"])
    print("SUBJECTS:", len(participants), "/37")
    if group:
        print("PO8_A_DELTA_uV:", group["PO8_A_mean_delta_uV"])
        print("PO8_A_BOOTSTRAP_CI95_uV:", group["PO8_A_bootstrap_ci95_uV"])
        print("PO8_B_DELTA_uV:", group["PO8_B_mean_delta_uV_descriptive"])
    print("REPORT:", report_path)
    print("OUTPUT:", out)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("[FAIL] P14:", type(exc).__name__ + ": " + str(exc),
              file=sys.stderr, flush=True)
        sys.exit(1)
