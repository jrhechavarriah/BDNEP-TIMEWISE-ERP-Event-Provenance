#!/usr/bin/env python3
"""Build Table 1/2/S1 from SHA256-verified full P14 replay, canonical BDNEP or P14 legacy.

--repo <BDNEP ROOT> or --replay-dir <P14 ROOT> are alternatives.
--out-dir must not exist; scientific inputs and references are never overwritten.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import statistics
from pathlib import Path
from bdnep_inputs import verified_artifact


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--repo', dest='replay_dir', type=Path)
    source.add_argument('--replay-dir', dest='replay_dir', type=Path)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    report = json.loads(verified_artifact(args.replay_dir, 'replay_verification_report.json').read_text(encoding='utf-8-sig'))
    if report.get('status') != 'FULL_INDEPENDENT_REPLAY_MATCHES_FROZEN_NUMERICAL_REFERENCES' or report.get('n_subjects') != 37:
        raise RuntimeError('Unverified full replay status')
    with verified_artifact(args.replay_dir, 'recomputed_participants.csv').open(newline='', encoding='utf-8-sig') as handle:
        rows = list(csv.DictReader(handle))
    if [r['subject_id'] for r in rows] != [f'sub-{i:03d}' for i in range(4, 41)]:
        raise RuntimeError('Unexpected cohort size/order')
    group = report['group']
    if group['n_A_stimuli'] != 11840 or group['n_B_stimuli'] != 10810:
        raise RuntimeError('Frozen trial totals differ')

    def mean(ch, arm, term):
        return statistics.mean(float(r[f'{ch}_{arm}_{term}_uV']) for r in rows)

    table2 = []
    for channel in ('PO8', 'PO7'):
        for arm in ('A', 'B'):
            original, corrected, difference = (mean(channel, arm, kind) for kind in ('original', 'corrected', 'delta'))
            if abs((corrected - original) - difference) > 1e-9:
                raise RuntimeError('Computed contrast mismatch')
            ci = group['PO8_A_bootstrap_ci95_uV'] if (channel, arm) == ('PO8', 'A') else None
            flips = sum(str(r.get(f'{channel}_{arm}_sign_flip', '')).lower() == 'true' for r in rows)
            table2.append({
                'Electrode / Arm': f'{channel} / {arm}',
                'Original mean (µV)': f'{original:.6f}',
                'Corrected mean (µV)': f'{corrected:.6f}',
                'Corrected − original (µV)': f'{difference:.6f}',
                '95% bootstrap CI (µV)': f'[{ci[0]:.6f}, {ci[1]:.6f}]' if ci else '—',
                'Contrast sign reversals (n)': str(flips),
                'Analysis role': 'Primary' if (channel, arm) == ('PO8', 'A') else 'Descriptive',
            })
    if abs(float(table2[0]['Corrected − original (µV)']) - group['PO8_A_mean_delta_uV']) > 1e-6:
        raise RuntimeError('Primary PO8 result changed')

    counts = report['arm_b_classification_counts']
    # Do not repeat the Arm B eligible total; it already appears as a named row.
    categories = [(key, value) for key, value in sorted(counts.items())
                  if key != 'BEHAVIORALLY_ELIGIBLE_PRE_ARTIFACT']
    if sum(counts.values()) != 11840 or counts['BEHAVIORALLY_ELIGIBLE_PRE_ARTIFACT'] != 10810:
        raise RuntimeError('Arm B classification accounting changed')
    table1 = [
        ('Participants', 37),
        ('Stimulus events (Arm A)', group['n_A_stimuli']),
        ('Eligible events (Arm B)', group['n_B_stimuli']),
        ('Ineligible stimulus events', group['n_A_stimuli'] - group['n_B_stimuli']),
        ('Orphan response events', report['orphan_responses']),
    ] + categories
    if args.out_dir.exists():
        raise RuntimeError('Refusing to overwrite existing tables output directory')
    args.out_dir.mkdir(parents=True)

    def save(name, content, fields):
        file = args.out_dir / name
        with file.open('x', newline='', encoding='utf-8-sig') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(content)
        print('BUILT:', file, 'SHA256:', digest(file))

    save('Table_1.csv', [{'Measure': key, 'Value': value} for key, value in table1], ['Measure', 'Value'])
    save('Table_2.csv', table2, list(table2[0]))
    save('Supplementary_Table_S1.csv', rows, list(rows[0]))


if __name__ == '__main__':
    main()
