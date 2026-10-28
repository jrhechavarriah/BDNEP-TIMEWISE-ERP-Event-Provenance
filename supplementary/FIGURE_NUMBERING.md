# Supplementary figure numbering — manuscript-aligned (P24)

The *final* numbering follows the order in the manuscript. Do not reuse the earlier labels.

| Figure | Scientific content | Script in `scr/` | Paper reference |
|---|---|---|---|
| Supplementary Figure S1 | All-condition behavioral eligibility under Arm B | `04_generate_supplementary_S1.py` | Results §3.3 |
| **Supplementary Figure S2** | Face/car behavioral eligibility; sub-016 = 2 face, 77 car | `05_generate_supplementary_S2_face_car.py` | Results §3.3 |
| **Supplementary Figure S3** | Unfiltered EEG anomalies sub-029 and sub-035 (six panels) | `08_generate_supplementary_S3_eeg_anomalies.py` | Results §3.5 |

S2 counts: 2,720 eligible face and 2,691 eligible car trials; 37 subjects. The 80-trial reference is per stimulus category.

S3: sub-029 151.0–162.5 s, sub-035 159.5–175.0 s. The 250,000 µV threshold is descriptive only, not an exclusion rule. Physical event onset calibration was not performed.

**Reproducibility:** S2 uses the verified participant-level replay CSV; S3 regenerates from checksum-verified NEMAR v1.1.1 raw EEG. Figure S1 and the frozen primary results are unchanged. All software is **internal until IP approval**.

## Final captions

**Supplementary Figure S2. Face and car trial availability under the behaviorally restricted Arm B analysis.** Numbers of eligible face (blue circles) and car (orange squares) trials by participant (n = 37). The dashed line marks 80 trials per condition. Eligibility required one correct response associated with the preceding stimulus and a reaction time of 200–1000 ms. There were 2,720 eligible face and 2,691 eligible car trials. Participant sub-016 had 2 eligible face and 77 eligible car trials. Scrambled conditions are not plotted. No post hoc exclusion or EEG artifact rejection was performed.

**Supplementary Figure S3. Representative unfiltered EEG signal anomalies and across-channel diagnostics.** Left: sub-029 (151.0–162.5 s); right: sub-035 (159.5–175.0 s). Panels display selected unfiltered channels, median and standard deviation across all 33 channels, and the number of channels exceeding the descriptive absolute-amplitude threshold of 250,000 µV. The threshold did not determine trial, channel, or participant exclusion. The figure is a descriptive quality-control visualization, not a physiological latency calibration or an alternative N170 result.
