"""Unit tests for supplementary figure builders; NO synthetic scientific release files."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

HERE = Path(__file__).resolve().parent

def load_script(stem):
    path = HERE / stem
    spec = importlib.util.spec_from_file_location(stem[:-3], path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

S1 = load_script("04_generate_supplementary_S1.py")
S2 = load_script("05_generate_supplementary_S2_face_car.py")
S3 = load_script("08_generate_supplementary_S3_eeg_anomalies.py")


class SupplementaryTests(unittest.TestCase):
    def test_S1_rejects_no_manifest(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(RuntimeError):
                S1.collect_trials(Path(d))

    def test_S1_rejects_fabricated_eligibility(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            rows = ("subject_id,n_A,n_B,n_B_face,n_B_car,n_B_scrambled_face,n_B_scrambled_car\n"
                    "sub-004,320,321,80,80,80,81\n")
            (d / "recomputed_participants.csv").write_text(rows)
            (d / "replay_verification_report.json").write_text(json.dumps({
                "status": S1.P14_PASS, "n_subjects": 37}))
            (d / "SHA256SUMS.txt").write_text("".join(
                S1.sha256(d / name) + "  " + name + "\n" for name in (
                    "recomputed_participants.csv", "replay_verification_report.json")))
            with self.assertRaises(RuntimeError):
                S1.collect_trials(d)


    def test_S2_is_face_car_not_eeg(self):
        self.assertEqual(S2.SCRIPT_NAME, "05_generate_supplementary_S2_face_car.py")
        self.assertEqual(S2.FIGURE_NAME, "Supplementary_Figure_S2_Face_Car_Behavioral_Eligibility")
        self.assertEqual(S2.EXPECTED_COUNTS["n_B_face"], 2720)
        self.assertEqual(S2.EXPECTED_COUNTS["n_B_car"], 2691)

    def test_S3_has_only_expected_participants(self):
        self.assertEqual(set(S3.WINDOWS_S), {"sub-029", "sub-035"})
        self.assertEqual(S3.DESCRIPTIVE_THRESHOLD_UV, 250000.0)

    def test_S3_plotting_synthetic_isolated(self):
        # Synthetic plotting contract only; this is not EEG evidence.
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "mock_S3.png"
            fake = {}
            for subject, window in S3.WINDOWS_S.items():
                t = np.linspace(*window, 150)
                curves = {name: np.sin(t) * 100 for name in S3.CHANNELS[subject]}
                fake[subject] = {"t": t, "selected": curves,
                                 "median": np.sin(t) * 30, "sd": np.abs(np.cos(t)) * 20,
                                 "exceed": np.zeros(len(t), dtype=int)}
            S3.render_s3(fake, out)
            self.assertTrue(out.is_file())
            self.assertGreater(out.stat().st_size, 5000)


if __name__ == "__main__":
    unittest.main()
